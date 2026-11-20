#include <AMReX.H>
#include <AMReX_Print.H>
#include "Utilities.H"
#include "Table.H"
#include "BlackBoxFunction.H"
#include "Filter.H"
#include <AMReX_MultiFab.H>
#include <AMReX_GpuAsyncArray.H>

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

namespace {

// Smooth field, periodic on the test domain, not a tensor product, and
// different for each component
void
fill_wavy(amrex::MultiFab& mf, const amrex::Box& domain)
{
  const auto& a = mf.arrays();
  const amrex::IntVect len = domain.length();
  const amrex::Real twopi = 2.0 * amrex::Math::pi<amrex::Real>();
  amrex::ParallelFor(
    mf, mf.nGrowVect(), mf.nComp(),
    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
      const amrex::Real x = twopi * i / len[0];
      const amrex::Real y = twopi * j / len[1];
      const amrex::Real z = (AMREX_SPACEDIM > 2) ? twopi * k / len[2] : 0.0;
      a[nbx](i, j, k, n) =
        std::sin(x + 2.0 * y + 3.0 * z + n) + std::cos(2.0 * x - y) * (n + 1);
    });
  amrex::Gpu::synchronize();
}

// Linear field with a different offset and slope for each component
void
fill_linear(amrex::MultiFab& mf)
{
  const auto& a = mf.arrays();
  amrex::ParallelFor(
    mf, mf.nGrowVect(), mf.nComp(),
    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
      a[nbx](i, j, k, n) = 10.0 * (n + 1) + (n + 1) * i - 2.0 * j + 3.0 * k;
    });
  amrex::Gpu::synchronize();
}

// Direct multi-dimensional convolution with the filter's 1D weights
void
reference_filter(Filter& f, const amrex::MultiFab& in, amrex::MultiFab& out)
{
  const auto& wv = f.get_weights();
  amrex::Gpu::AsyncArray<amrex::Real> weights(wv.data(), wv.size());
  const amrex::Real* w = weights.data();
  const int ng = f.get_filter_ngrow();
  const auto& q = in.const_arrays();
  const auto& qh = out.arrays();
  out.setVal(0.0, 0, out.nComp(), out.nGrow());
  amrex::ParallelFor(
    out, out.nGrowVect(), out.nComp(),
    [=] AMREX_GPU_DEVICE(int b, int i, int j, int k, int n) noexcept {
      run_filter(i, j, k, n, ng, 0, w, q[b], qh[b]);
    });
  amrex::Gpu::synchronize();
}

amrex::Real
max_diff(
  const amrex::MultiFab& a,
  const amrex::MultiFab& b,
  int comp,
  int n,
  int nghost = -1)
{
  if (nghost < 0) {
    nghost = a.nGrow();
  }
  amrex::MultiFab d(a.boxArray(), a.DistributionMap(), n, nghost);
  amrex::MultiFab::Copy(d, a, comp, 0, n, nghost);
  amrex::MultiFab::Subtract(d, b, comp, 0, n, nghost);
  amrex::Real m = 0.0;
  for (int c = 0; c < n; ++c) {
    m = amrex::max(m, d.norm0(c, nghost));
  }
  return m;
}

void
test_filters()
{
  const int ncomp = 3;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(31, 23, 15)));
  amrex::BoxArray ba(domain);
  ba.maxSize(8);
  const amrex::DistributionMapping dm(ba);
  const int out_ng = 1;

  for (int type = no_filter; type < num_filter_types; ++type) {
    for (int fgr : {2, 4, 8}) {
      Filter f(type, fgr);
      const int ng = f.get_filter_ngrow();
      // Round-off grows with the sum of |weight| over the stencil, which is
      // large for the polynomial approximations at large fgr
      amrex::Real wsum = 0.0;
      for (const auto x : f.get_weights()) {
        wsum += std::abs(x);
      }
      const amrex::Real tol = 1e-12 * std::pow(wsum, AMREX_SPACEDIM);
      amrex::MultiFab in(ba, dm, ncomp, out_ng + ng);
      amrex::MultiFab out(ba, dm, ncomp, out_ng);
      amrex::MultiFab ref(ba, dm, ncomp, out_ng);

      // Weights are symmetric and sum to one: a linear field is unchanged
      fill_linear(in);
      f.apply_filter(in, out);
      AMREX_ALWAYS_ASSERT(max_diff(out, in, 0, ncomp) < tol);

      // Same as the direct multi-dimensional convolution, ghost cells too
      fill_wavy(in, domain);
      f.apply_filter(in, out);
      reference_filter(f, in, ref);
      AMREX_ALWAYS_ASSERT(max_diff(out, ref, 0, ncomp) < tol);

      // A component sub-range gives the same result as filtering everything
      // and leaves the other components alone
      {
        amrex::MultiFab part(ba, dm, ncomp, out_ng);
        part.setVal(-123.0);
        f.apply_filter(in, part, 1, 2);
        AMREX_ALWAYS_ASSERT(max_diff(part, ref, 1, 2) < tol);
        AMREX_ALWAYS_ASSERT(part.min(0, out_ng) == -123.0);
        AMREX_ALWAYS_ASSERT(part.max(0, out_ng) == -123.0);
      }

      // FArrayBox interface, on a sub-box and a component sub-range
      for (amrex::MFIter mfi(in); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), out_ng);
        amrex::FArrayBox fout(bx, ncomp);
        fout.setVal<amrex::RunOn::Device>(-123.0);
        f.apply_filter(bx, in[mfi], fout, 1, 2);
        amrex::Gpu::streamSynchronize();
        amrex::FArrayBox diff(bx, 2);
        diff.copy<amrex::RunOn::Device>(fout, bx, 1, bx, 0, 2);
        diff.minus<amrex::RunOn::Device>(ref[mfi], bx, bx, 1, 0, 2);
        AMREX_ALWAYS_ASSERT(diff.norm<amrex::RunOn::Device>(bx, 0, 0, 2) < tol);
        AMREX_ALWAYS_ASSERT(fout.max<amrex::RunOn::Device>(bx, 0) == -123.0);
      }
    }
  }
}

} // namespace

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  // Test table interpolation
  pele::physics::TabulatedFunctionParams tabfunc_par;
  tabfunc_par.host_only_parm().parm_parse_prefix = "test";
  tabfunc_par.initialize();
  pele::physics::TabulatedFunction tabfunc(
    static_cast<const pele::physics::TabulatedFunctionData*>(
      &tabfunc_par.host_parm()));
  amrex::Real outval;
  const amrex::GpuArray<amrex::Real, 2> invals{1.0, 1.0};
  const std::string varname = "T";
  const int idxT =
    pele::physics::get_var_index(varname.c_str(), &tabfunc_par.host_parm());
  tabfunc.get_value(idxT, invals.data(), outval);
  AMREX_ALWAYS_ASSERT(amrex::almostEqual(outval, 3.60));
  tabfunc_par.deallocate();
  amrex::Print() << "TabulatedFunction Tests passed" << std::endl;

  // Test Unit Conversions
  AMREX_ALWAYS_ASSERT(m2c::Length(1.0) == 100.0);
  AMREX_ALWAYS_ASSERT(c2m::Length(1.0, 2) == 1.0e-4);
  AMREX_ALWAYS_ASSERT(m2c::Mass(c2m::Mass(1.0)) == 1.0);
  AMREX_ALWAYS_ASSERT(m2c::Mu(1.0) == m2c::Mass(1.0) / m2c::Length(1.0));
  amrex::Print() << "Unit Conversion Tests Passed" << std::endl;

  // Test locate function
  amrex::Array<amrex::Real, 5> sorted_data{0.0, 0.2, 0.5, 0.9, 1.0};
  int length = sorted_data.size();
  int idxlo = -1;
  amrex::Real x = -50.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 0);
  x = 1.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == length - 1);
  x = 0.25;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 1);
  amrex::Print() << "Locate Tests Passed" << std::endl;

  // Test rectangle circle intersection area
  // Quarter circle intersection
  amrex::Real x0 = 0.0;
  amrex::Real x1 = 1.0;
  amrex::Real y0 = 0.0;
  amrex::Real y1 = 1.0;
  amrex::Real cx = 1.0;
  amrex::Real cy = 1.0;
  amrex::Real r = 1.0;
  amrex::Real A1 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A1 - 0.25 * M_PI) < 1e-14);
  // No intersection
  cx = 3.0;
  cy = 3.0;
  amrex::Real A2 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A2) < 1e-14);
  // Square covered by circle
  r = 20.0;
  amrex::Real A3 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A3 - 1.0) < 1e-14);
  amrex::Print() << "Intersection Area Tests Passed" << std::endl;

  test_filters();
  amrex::Print() << "Filter Tests Passed" << std::endl;

  amrex::Finalize();
  return 0;
}
