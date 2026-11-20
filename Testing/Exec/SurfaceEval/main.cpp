#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Print.H>

#include "mechanism.H"
#include "SurfaceChemistry.H"
#include "reference.H"

#if NUM_SURFACE_SPECIES == 0
#error "SurfaceEval needs a heterogeneous mechanism (Chemistry_Model=ptcombust)"
#endif

namespace {

// Max relative error, ignoring entries negligible against the largest one.
amrex::Real
maxRelErr(const amrex::Real* actual, const amrex::Real* expected, const int n)
{
  amrex::Real scale = 1.0e-300;
  for (int i = 0; i < n; ++i) {
    scale = amrex::max(scale, std::abs(expected[i]));
  }
  scale *= 1.0e-12;
  amrex::Real err = 0.0;
  for (int i = 0; i < n; ++i) {
    err = amrex::max(
      err, std::abs(actual[i] - expected[i]) / (std::abs(expected[i]) + scale));
  }
  return err;
}

amrex::Real
maxAbsErr(const amrex::Real* actual, const amrex::Real* expected, const int n)
{
  amrex::Real err = 0.0;
  for (int i = 0; i < n; ++i) {
    err = amrex::max(err, std::abs(actual[i] - expected[i]));
  }
  return err;
}

} // namespace

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    constexpr int ng = NUM_GAS_SPECIES;
    constexpr int ns = NUM_SURFACE_SPECIES;
    static_assert(reference::n_gas == ng && reference::n_surf == ns);

    // Copy the reference inputs into arrays the device lambda can capture.
    amrex::GpuArray<amrex::Real, ng> Y;
    amrex::GpuArray<amrex::Real, ns> theta, theta0;
    for (int k = 0; k < ng; ++k) {
      Y[k] = reference::Y[k];
    }
    for (int k = 0; k < ns; ++k) {
      theta[k] = reference::theta[k];
      theta0[k] = reference::theta0[k];
    }
    const amrex::Real rho = reference::rho;
    const amrex::Real T = reference::T;

    // Output layout: flux, dthetadt, theta_ss from a bare surface, flux_ss,
    // theta_ss from uniform coverages, converged flags
    constexpr int i_flux = 0;
    constexpr int i_dthetadt = i_flux + ng;
    constexpr int i_theta_ss = i_dthetadt + ns;
    constexpr int i_flux_ss = i_theta_ss + ns;
    constexpr int i_theta_ss_u = i_flux_ss + ng;
    constexpr int i_ok = i_theta_ss_u + ns;
    amrex::Gpu::DeviceVector<amrex::Real> d_out(i_ok + 1);
    amrex::Real* out = d_out.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*i*/) noexcept {
      namespace surf = pele::physics::surface;
      surf::surfaceRates(
        rho, T, Y.data(), theta.data(), out + i_flux, out + i_dthetadt);

      amrex::Real th[ns], dthetadt_ss[ns];
      for (int k = 0; k < ns; ++k) {
        th[k] = theta0[k];
      }
      const bool ok = surf::steadyCoverages(rho, T, Y.data(), th);
      surf::surfaceRates(rho, T, Y.data(), th, out + i_flux_ss, dthetadt_ss);
      for (int k = 0; k < ns; ++k) {
        out[i_theta_ss + k] = th[k];
      }

      // a poor initial guess must reach the same steady state
      for (int k = 0; k < ns; ++k) {
        th[k] = 1.0 / ns;
      }
      const bool ok_u = surf::steadyCoverages(rho, T, Y.data(), th);
      for (int k = 0; k < ns; ++k) {
        out[i_theta_ss_u + k] = th[k];
      }
      out[i_ok] = (ok && ok_u) ? 1.0 : 0.0;
    });

    amrex::Vector<amrex::Real> h_out(d_out.size());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, d_out.begin(), d_out.end(), h_out.begin());
    const amrex::Real* res = h_out.data();

    const amrex::Real err_flux = maxRelErr(res + i_flux, reference::flux, ng);
    const amrex::Real err_dthetadt =
      maxRelErr(res + i_dthetadt, reference::dthetadt, ns);
    const amrex::Real err_theta_ss =
      maxAbsErr(res + i_theta_ss, reference::theta_ss, ns);
    const amrex::Real err_flux_ss =
      maxRelErr(res + i_flux_ss, reference::flux_ss, ng);
    const amrex::Real err_theta_ss_u =
      maxAbsErr(res + i_theta_ss_u, reference::theta_ss, ns);

    amrex::Print() << "flux      max rel err: " << err_flux << "\n"
                   << "dthetadt  max rel err: " << err_dthetadt << "\n"
                   << "steady coverages converged: " << (res[i_ok] > 0.0)
                   << "\n"
                   << "theta_ss  max abs err: " << err_theta_ss << "\n"
                   << "flux_ss   max rel err: " << err_flux_ss << "\n"
                   << "theta_ss (uniform start) max abs err: " << err_theta_ss_u
                   << "\n";

    const bool pass = err_flux < 1.0e-10 && err_dthetadt < 1.0e-10 &&
                      res[i_ok] > 0.0 && err_theta_ss < 1.0e-8 &&
                      err_flux_ss < 1.0e-6 && err_theta_ss_u < 1.0e-8;
    if (!pass) {
      amrex::Abort("SurfaceEval: surface rates do not match Cantera");
    }
    amrex::Print() << "SurfaceEval: PASS\n";
  }
  amrex::Finalize();
  return 0;
}
