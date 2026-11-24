#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>
#include <TestCheck.H>

namespace {
// The checks that apply to the EOS, with a tolerance on the largest relative
// error over the (T, Y, p) grid: 10x or more the measured maximum (gnu,
// double precision, in the comments). A check tagged with the fix it needs
// is an expected failure.
const char* const eos_check_names[eoscheck::NUM] = {
  "T -> e -> T (RTY2E, REY2T)",
  "T -> h -> T (RTY2H, RHY2T)",
  "p,T -> rho -> p (PYT2R, RTY2P)",
  "rho,p -> T (RYP2T)",
  "rho,e -> T (RYET2P)",
  "rho,e -> p (RYET2P)",
  "h - e = p/rho",
  "cv = de/dT at constant rho",
  "cp = dh/dT at constant p",
  "cp - cv = R/W"};

struct EosCheckSpec
{
  int comp;
  amrex::Real tol;
  const char* tag = "";
};

// Iterative inversions converge to an absolute T tolerance (Fuego: Newton
// step < 1e-6 K); the finite differences use dT = 1e-4 T.
const EosCheckSpec eos_check_specs[] = {
#if defined(USE_MANIFOLD_EOS)
  // p(rho(Y)) = nominal p (T is not used)
  {eoscheck::PRP, 1.e-13}, // 1.2e-16 (table and network)
#elif defined(USE_SRK_EOS)
  // No h inversion (RHY2T aborts), and cp - cv != R/W for a real gas
  {eoscheck::TET, 1.e-12},  // 1.2e-14
  {eoscheck::PRP, 1.e-13},  // 5.9e-16
  {eoscheck::RPT, 2.e-10},  // 1.3e-11
  {eoscheck::REPT, 1.e-12}, // 1.2e-14
  {eoscheck::REPP, 1.e-12}, // 1.2e-14
  {eoscheck::HE, 1.e-13},   // 2.8e-15
  {eoscheck::CVFD, 2.e-8},  // 1.6e-9
  {eoscheck::CPFD, 2.e-8},  // 1.2e-9
#elif defined(USE_GAMMALAW_EOS)
  {eoscheck::TET, 1.e-13}, // 2.7e-16
  {eoscheck::THT, 1.e-13}, // 1.9e-16
  {eoscheck::PRP, 1.e-13}, // 2.9e-16
  {eoscheck::RPT, 1.e-13}, // 2.6e-16
  {eoscheck::REPT, 1.e-13,
   "bug: GammaLaw::RYET2P takes T by value and does not return it"},
  {eoscheck::REPP, 1.e-13}, // 2.9e-16
  {eoscheck::HE, 1.e-13},   // 4.9e-16
  {eoscheck::CVFD, 3.e-11}, // 1.5e-12
  {eoscheck::CPFD, 3.e-11}, // 2.1e-12
  {eoscheck::CPCV, 1.e-13}, // 1.6e-16
#else
  {eoscheck::TET, 2.e-8},   // 1.6e-9
  {eoscheck::THT, 2.e-8},   // 1.4e-9
  {eoscheck::PRP, 1.e-13},  // 3.7e-16
  {eoscheck::RPT, 1.e-13},  // 5.0e-16
  {eoscheck::REPT, 2.e-8},  // 1.6e-9
  {eoscheck::REPP, 2.e-8},  // 1.6e-9
  {eoscheck::HE, 1.e-12},   // 1.5e-14
  {eoscheck::CVFD, 2.e-8},  // 1.6e-9
  {eoscheck::CPFD, 2.e-8},  // 1.2e-9
  {eoscheck::CPCV, 1.e-13}, // 3.7e-15
#endif
};
} // namespace

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  int status = 0;

  {
    pele::physics::PeleParams<
      pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>>
      eos_parms;
    amrex::Print() << " Initialization of EOS (CPP)... \n";
    eos_parms.initialize();
    auto const* leosparm = eos_parms.device_parm();

    amrex::ParmParse pp;

    int gridsize = 128;
    pp.query("gridsize", gridsize);

    // Define geometry
    amrex::Array<int, AMREX_SPACEDIM> npts{AMREX_D_DECL(1, 1, 1)};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      npts[i] = gridsize;
    }

    amrex::Box domain(
      amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
      amrex::IntVect(AMREX_D_DECL(npts[0] - 1, npts[1] - 1, npts[2] - 1)));

    amrex::RealBox real_box(
      {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

    int coord = 0;

    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    amrex::Geometry geom(domain, real_box, coord, is_periodic);

    // Define BoxArray
    int max_size = 32;
    pp.query("max_size", max_size);
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);

    amrex::ParmParse ppa("amr");
    bool do_plot = true;
    ppa.query("do_plot", do_plot);
    std::string pltfile("plt");
    ppa.query("plot_file", pltfile);

    amrex::DistributionMapping dm{ba};
    int num_grow = 0;

    // Data MFs
    amrex::MultiFab mass_frac(ba, dm, NUM_SPECIES, num_grow);
    amrex::MultiFab temperature(ba, dm, 1, num_grow);
    amrex::MultiFab density(ba, dm, 1, num_grow);
    amrex::MultiFab energy(ba, dm, 1, num_grow);

    const auto geomdata = geom.data();
    {
      BL_PROFILE("Pele::init()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        auto const& Y_a = mass_frac.array(mfi);
        auto const& T_a = temperature.array(mfi);
        auto const& rho_a = density.array(mfi);
        auto const& e_a = energy.array(mfi);
        amrex::ParallelFor(
          bx, [Y_a, T_a, rho_a, e_a, geomdata,
               leosparm] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            initialize_data(i, j, k, Y_a, T_a, rho_a, e_a, geomdata, leosparm);
          });
      }
    }
    amrex::MultiFab cp(ba, dm, 1, num_grow);
    amrex::MultiFab wdot(ba, dm, NUM_SPECIES, num_grow);
    {
      BL_PROFILE("Pele::get_cp()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& box = mfi.tilebox();

        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& T_a = temperature.const_array(mfi);
        auto const& cp_a = cp.array(mfi);
        auto const& rho_a = density.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            get_cp(i, j, k, Y_a, T_a, rho_a, cp_a, leosparm);
          });
      }
    }
    {
      BL_PROFILE("Pele::get_wdot()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& box = mfi.tilebox();

        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& T_a = temperature.const_array(mfi);
        auto const& wdot_a = wdot.array(mfi);
        auto const& rho_a = density.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifndef USE_GAMMALAW_EOS
            get_wdot(i, j, k, Y_a, T_a, rho_a, wdot_a, leosparm);
#else
            // GammaLaw: Wdot not supported, just set to 0
            wdot_a(i,j,k,0) = 0.0;
#endif
          });
      }
    }

    {
      BL_PROFILE("Pele::get_T()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& box = mfi.tilebox();

        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& e_a = energy.const_array(mfi);
        auto const& T_a = temperature.array(mfi);
        auto const& rho_a = density.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            get_T_from_EY(i, j, k, Y_a, T_a, rho_a, e_a, leosparm);
          });
      }
    }
    if (do_plot) {
      amrex::MultiFab VarPlt(ba, dm, 4 + 2 * NUM_SPECIES, num_grow);
      amrex::MultiFab::Copy(VarPlt, density, 0, 0, 1, num_grow);
      amrex::MultiFab::Copy(VarPlt, mass_frac, 0, 1, NUM_SPECIES, num_grow);
      amrex::MultiFab::Copy(VarPlt, cp, 0, 1 + NUM_SPECIES, 1, num_grow);
      amrex::MultiFab::Copy(
        VarPlt, wdot, 0, 2 + NUM_SPECIES, NUM_SPECIES, num_grow);
      amrex::MultiFab::Copy(
        VarPlt, temperature, 0, 2 + 2 * NUM_SPECIES, 1, num_grow);
      amrex::MultiFab::Copy(
        VarPlt, energy, 0, 3 + 2 * NUM_SPECIES, 1, num_grow);

      std::string outfile = amrex::Concatenate(pltfile, 1);
      amrex::Vector<std::string> plt_VarsName, spec_VarsName;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
        spec_VarsName);

      plt_VarsName.push_back("rho");
      for (int n = 0; n < NUM_SPECIES; ++n) {
        plt_VarsName.push_back("Y(" + spec_VarsName[n] + ")");
      }
      plt_VarsName.push_back("cp");
      for (int n = 0; n < NUM_SPECIES; ++n) {
        plt_VarsName.push_back("wdot(" + spec_VarsName[n] + ")");
      }
      plt_VarsName.push_back("temperature");
      plt_VarsName.push_back("energy");

      amrex::WriteSingleLevelPlotfile(
        outfile, VarPlt, plt_VarsName, geom, 0.0, 0);
    }

    {
      BL_PROFILE("Pele::eos_checks()");
      amrex::Print() << "\n-- EOS consistency checks ("
                     << pele::physics::PhysicsType::eos_type::identifier()
                     << ") --\n";
      amrex::Real T_lo = 300.0;
      amrex::Real T_hi = 3000.0;
      amrex::ParmParse ppe("eos_eval");
      ppe.query("T_lo", T_lo);
      ppe.query("T_hi", T_hi);
      amrex::MultiFab err(ba, dm, eoscheck::NUM, num_grow);
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& err_a = err.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            eos_checks(i, j, k, Y_a, err_a, geomdata, T_lo, T_hi, leosparm);
          });
      }
      pele::testing::Checker check;
      for (const auto& spec : eos_check_specs) {
        check.within(
          eos_check_names[spec.comp], err.max(spec.comp), spec.tol, spec.tag);
      }
      status = check.summary("EosEval");
    }

    eos_parms.deallocate();
  }

  amrex::Finalize();

  return status;
}
