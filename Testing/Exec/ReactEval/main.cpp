#include <cmath>
#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <PltFileManager.H>

#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>
#include <ReactorBase.H>
#include <TestCheck.H>

namespace {
const std::string level_prefix{"Level_"};
}

void
GotoNextLine(std::istream& is)
{
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}

namespace {

// Integrates every level over dt with the given reactor, in place
void
advance(
  pele::physics::reactions::ReactorBase& reactor,
  const int finest_level,
  amrex::Vector<amrex::MultiFab>& mf,
  amrex::Vector<amrex::MultiFab>& mfE,
  amrex::Vector<amrex::MultiFab>& rY_source_ext,
  amrex::Vector<amrex::MultiFab>& rY_source_energy_ext,
  amrex::Vector<amrex::MultiFab>& fctCount,
  amrex::Vector<amrex::iMultiFab>& dummyMask,
  const int reactFunc,
  int ode_ncells,
  const amrex::Real dt,
  const int ndt,
  const bool subcycledt)
{
  BL_PROFILE_VAR_NS("React", ReactInLoop);
  BL_PROFILE_VAR_NS("Allocs", Allocs);
  BL_PROFILE_VAR_NS("Flatten", mainflatten);

  for (int lev = 0; lev <= finest_level; ++lev) {
    amrex::Real lvl_strt = amrex::ParallelDescriptor::second();
    BL_PROFILE_VAR("Advance_Level" + std::to_string(lev), Advance);
#ifdef AMREX_USE_OMP
    const auto tiling = amrex::MFItInfo().SetDynamic(true);
#pragma omp parallel
#else
    const bool tiling = amrex::TilingIfNotGPU();
#endif
    for (amrex::MFIter mfi(mf[lev], tiling); mfi.isValid(); ++mfi) {

      int omp_thread = 0;
#ifdef AMREX_USE_OMP
      omp_thread = omp_get_thread_num();
#endif
      const amrex::Box& box = mfi.tilebox();
      const int nc = box.numPts();
      int extra_cells = 0;

      auto const& rhoY = mf[lev].array(mfi);
      auto const& T = mf[lev].array(mfi, NUM_SPECIES);
      auto const& rhoE = mfE[lev].array(mfi);
      auto const& frcExt = rY_source_ext[lev].array(mfi);
      auto const& frcEExt = rY_source_energy_ext[lev].array(mfi);
      auto const& fc = fctCount[lev].array(mfi);
      auto const& mask = dummyMask[lev].array(mfi);

      // -------------------------------------------------------------
      // Integration with Array4 react function
      if (reactFunc == 1) {
        amrex::Real time = 0.0;
        amrex::Real dt_lev = (subcycledt) ? dt / std::pow(2, lev) : dt;
        amrex::Real dt_incr = dt_lev / ndt;
        if (omp_thread == 0) {
          amrex::Print() << "  [" << lev << "]"
                         << " integrating " << nc << " cells \n";
        }
        /* Solve */
        BL_PROFILE_VAR_START(ReactInLoop);
        for (int ii = 0; ii < ndt; ++ii) {
          reactor.react(
            box, rhoY, frcExt, T, rhoE, frcEExt, fc, mask, dt_incr, time
#ifdef AMREX_USE_GPU
            ,
            amrex::Gpu::gpuStream()
#endif
          );
          // react() returns the time actually integrated in dt_incr
          dt_incr = dt_lev / ndt;
          amrex::Gpu::Device::streamSynchronize();
        }
        BL_PROFILE_VAR_STOP(ReactInLoop);

        // -------------------------------------------------------------
        // Integration with 1dArray raw pointer react function
      } else if (reactFunc == 2) {

        // On GPU, integrate the entirely box at once
        // othewise use the user-input ode_ncells
#ifdef AMREX_USE_GPU
        ode_ncells = nc;
#endif
        // Padding up to a multiple of ode_ncells
        extra_cells = (ode_ncells - nc % ode_ncells) % ode_ncells;
        if (omp_thread == 0) {
          amrex::Print() << " Integrating " << nc << " cells with a "
                         << ode_ncells << " ode cell buffer ";
          amrex::Print() << "(" << extra_cells << " extra cells) \n";
        }

        BL_PROFILE_VAR_START(Allocs);
        int nCells = nc + extra_cells;

#ifdef AMREX_USE_GPU
        auto tmp_vect_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
          nCells * (NUM_SPECIES + 1) * sizeof(amrex::Real));
        auto tmp_src_vect_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
          nCells * NUM_SPECIES * sizeof(amrex::Real));
        auto tmp_vect_energy_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
          nCells * sizeof(amrex::Real));
        auto tmp_src_vect_energy_d =
          (amrex::Real*)amrex::The_Device_Arena()->alloc(
            nCells * sizeof(amrex::Real));
        auto tmp_fc_d = (long int*)amrex::The_Device_Arena()->alloc(
          nCells * sizeof(long int));
#endif

        auto tmp_vect = new amrex::Real[nCells * (NUM_SPECIES + 1)];
        auto tmp_src_vect = new amrex::Real[nCells * NUM_SPECIES];
        auto tmp_vect_energy = new amrex::Real[nCells];
        auto tmp_src_vect_energy = new amrex::Real[nCells];
        auto tmp_fc = new long int[nCells];

        BL_PROFILE_VAR_STOP(Allocs);

        BL_PROFILE_VAR_START(mainflatten);
#ifndef AMREX_USE_GPU
        reactor.flatten(
          box, nCells, rhoY, frcExt, T, rhoE, frcEExt, tmp_vect, tmp_src_vect,
          tmp_vect_energy, tmp_src_vect_energy);

        // The padding cells repeat the box's first cell
        const auto lo = amrex::lbound(box);
        for (int icell = nc; icell < nc + extra_cells; icell++) {
          for (int sp = 0; sp < NUM_SPECIES; sp++) {
            tmp_vect[icell * (NUM_SPECIES + 1) + sp] =
              rhoY(lo.x, lo.y, lo.z, sp);
            tmp_src_vect[icell * NUM_SPECIES + sp] =
              frcExt(lo.x, lo.y, lo.z, sp);
          }
          tmp_vect[icell * (NUM_SPECIES + 1) + NUM_SPECIES] =
            T(lo.x, lo.y, lo.z);
          tmp_vect_energy[icell] = rhoE(lo.x, lo.y, lo.z);
          tmp_src_vect_energy[icell] = frcEExt(lo.x, lo.y, lo.z);
        }
#else
        reactor.flatten(
          box, nCells, rhoY, frcExt, T, rhoE, frcEExt, tmp_vect_d,
          tmp_src_vect_d, tmp_vect_energy_d, tmp_src_vect_energy_d);

        amrex::Gpu::copy(
          amrex::Gpu::deviceToHost, tmp_vect_d,
          tmp_vect_d + nCells * (NUM_SPECIES + 1), tmp_vect);
        amrex::Gpu::copy(
          amrex::Gpu::deviceToHost, tmp_src_vect_d,
          tmp_src_vect_d + nCells * NUM_SPECIES, tmp_src_vect);
        amrex::Gpu::copy(
          amrex::Gpu::deviceToHost, tmp_vect_energy_d,
          tmp_vect_energy_d + nCells, tmp_vect_energy);
        amrex::Gpu::copy(
          amrex::Gpu::deviceToHost, tmp_src_vect_energy_d,
          tmp_src_vect_energy_d + nCells, tmp_src_vect_energy);
#endif
        BL_PROFILE_VAR_STOP(mainflatten);

        /* Solve */
        BL_PROFILE_VAR_START(ReactInLoop);
        for (int i = 0; i < nCells; i += ode_ncells) {
          tmp_fc[i] = 0;
          amrex::Real time = 0.0;
          amrex::Real dt_lev = (subcycledt) ? dt / std::pow(2, lev) : dt;
          amrex::Real dt_incr = dt_lev / ndt;
          for (int ii = 0; ii < ndt; ++ii) {
            tmp_fc[i] += reactor.react(
              &tmp_vect[i * (NUM_SPECIES + 1)], &tmp_src_vect[i * NUM_SPECIES],
              &tmp_vect_energy[i], &tmp_src_vect_energy[i], dt_incr, time,
              ode_ncells
#ifdef AMREX_USE_GPU
              ,
              amrex::Gpu::gpuStream()
#endif
            );

            // react() returns the time actually integrated in dt_incr
            dt_incr = dt_lev / ndt;
            for (int ic = i + 1; ic < i + ode_ncells; ++ic) {
              tmp_fc[ic] = tmp_fc[i];
            }
            amrex::Gpu::Device::streamSynchronize();
          }
        }
        BL_PROFILE_VAR_STOP(ReactInLoop);

        BL_PROFILE_VAR_START(mainflatten);
#ifndef AMREX_USE_GPU
        reactor.unflatten(
          box, nCells, rhoY, T, rhoE, frcEExt, fc, tmp_vect, tmp_vect_energy,
          tmp_fc, dt);
#else

        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, tmp_vect,
          tmp_vect + nCells * (NUM_SPECIES + 1), tmp_vect_d);
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, tmp_src_vect,
          tmp_src_vect + nCells * NUM_SPECIES, tmp_src_vect_d);
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, tmp_vect_energy, tmp_vect_energy + nCells,
          tmp_vect_energy_d);
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, tmp_src_vect_energy,
          tmp_src_vect_energy + nCells, tmp_src_vect_energy_d);
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, tmp_fc, tmp_fc + nCells, tmp_fc_d);

        reactor.unflatten(
          box, nCells, rhoY, T, rhoE, frcEExt, fc, tmp_vect_d,
          tmp_vect_energy_d, tmp_fc_d, dt);
#endif
        BL_PROFILE_VAR_STOP(mainflatten);

        delete[] tmp_vect;
        delete[] tmp_src_vect;
        delete[] tmp_vect_energy;
        delete[] tmp_src_vect_energy;
        delete[] tmp_fc;
#ifdef AMREX_USE_GPU
        amrex::The_Device_Arena()->free(tmp_vect_d);
        amrex::The_Device_Arena()->free(tmp_src_vect_d);
        amrex::The_Device_Arena()->free(tmp_vect_energy_d);
        amrex::The_Device_Arena()->free(tmp_src_vect_energy_d);
        amrex::The_Device_Arena()->free(tmp_fc_d);
#endif
      }
    }
    BL_PROFILE_VAR_STOP(Advance);
    amrex::Real lvl_run_time = amrex::ParallelDescriptor::second() - lvl_strt;
    amrex::ParallelDescriptor::ReduceRealMax(
      lvl_run_time, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "   >> Level " << lev << " advance: " << lvl_run_time
                   << "\n";
  }
}

// Largest differences over all cells between two final states, computed on
// the host from copies of the data (the solution may live on the device)
struct StateDiff
{
  amrex::Real T = 0.0;   // relative temperature difference
  amrex::Real Y = 0.0;   // absolute mass fraction difference
  amrex::Real rho = 0.0; // relative density difference
  amrex::Real Te = 0.0;  // relative difference between T and T(e or h, Y)
};

// Copy of src's valid data in pinned host memory
amrex::MultiFab
toHost(const amrex::MultiFab& src)
{
  amrex::MFInfo pinned;
  pinned.SetArena(amrex::The_Pinned_Arena());
  amrex::MultiFab h(
    src.boxArray(), src.DistributionMap(), src.nComp(), 0, pinned);
  amrex::dtoh_memcpy(h, src);
  amrex::Gpu::streamSynchronize();
  return h;
}

// Copy of each level of src (valid cells), in the default arena
amrex::Vector<amrex::MultiFab>
copyOf(const amrex::Vector<amrex::MultiFab>& src)
{
  amrex::Vector<amrex::MultiFab> dst(src.size());
  for (int lev = 0; lev < static_cast<int>(src.size()); ++lev) {
    dst[lev].define(
      src[lev].boxArray(), src[lev].DistributionMap(), src[lev].nComp(), 0);
    amrex::MultiFab::Copy(dst[lev], src[lev], 0, 0, src[lev].nComp(), 0);
  }
  return dst;
}

// mf, mf_ref: rhoY and T of the state and of the reference; mf0, mfE0: the
// initial rhoY and rho*e (reactor type 1) or rho*h (type 2), which a closed
// reactor with no forcing conserves.
StateDiff
compareStates(
  const int finest_level,
  const amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<amrex::MultiFab>& mf_ref,
  const amrex::Vector<amrex::MultiFab>& mf0,
  const amrex::Vector<amrex::MultiFab>& mfE0,
  const int reactor_type)
{
  const auto upd = [](amrex::Real& m, const amrex::Real v) {
    m = amrex::max(m, pele::testing::finiteOrInf(v));
  };
  auto eos = pele::physics::PhysicsType::eos();
  StateDiff d;
  for (int lev = 0; lev <= finest_level; ++lev) {
    const amrex::MultiFab s_h = toHost(mf[lev]);
    const amrex::MultiFab r_h = toHost(mf_ref[lev]);
    const amrex::MultiFab s0_h = toHost(mf0[lev]);
    const amrex::MultiFab e0_h = toHost(mfE0[lev]);
    for (amrex::MFIter mfi(s_h); mfi.isValid(); ++mfi) {
      const auto s = s_h.const_array(mfi);
      const auto r = r_h.const_array(mfi);
      const auto s0 = s0_h.const_array(mfi);
      const auto e0 = e0_h.const_array(mfi);
      amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
        amrex::Real rho = 0.0, rho_ref = 0.0, rho0 = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          rho += s(i, j, k, n);
          rho_ref += r(i, j, k, n);
          rho0 += s0(i, j, k, n);
        }
        amrex::Real Y[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
          Y[n] = s(i, j, k, n) / rho;
          upd(d.Y, std::abs(Y[n] - r(i, j, k, n) / rho_ref));
        }
        const amrex::Real T = s(i, j, k, NUM_SPECIES);
        const amrex::Real T_ref = r(i, j, k, NUM_SPECIES);
        upd(d.T, std::abs(T - T_ref) / T_ref);
        upd(d.rho, std::abs(rho - rho0) / rho0);
        amrex::Real Te = T;
        if (reactor_type == 1) {
          eos.REY2T(rho, e0(i, j, k) / rho, Y, Te);
        } else {
          eos.RHY2T(rho, e0(i, j, k) / rho, Y, Te);
        }
        upd(d.Te, std::abs(Te - T) / T);
      });
    }
  }
  amrex::Real dmax[4] = {d.T, d.Y, d.rho, d.Te};
  amrex::ParallelDescriptor::ReduceRealMax(dmax, 4);
  return {dmax[0], dmax[1], dmax[2], dmax[3]};
}

// Integrates the reference from the initial data (mf0, mfE0 and the
// forcing; overwritten) and checks the final state mf against it, and
// against the conservation of mass and energy. Returns the exit code.
int
checkFinalState(
  const std::string& chem_integrator,
  const int reactor_type,
  const int finest_level,
  const amrex::Vector<amrex::MultiFab>& mf,
  amrex::Vector<amrex::MultiFab>& mf0,
  amrex::Vector<amrex::MultiFab>& mfE0,
  amrex::Vector<amrex::MultiFab>& rY_source_ext0,
  amrex::Vector<amrex::MultiFab>& rY_source_energy_ext0,
  amrex::Vector<amrex::MultiFab>& fctCount,
  amrex::Vector<amrex::iMultiFab>& dummyMask,
  const amrex::Real dt,
  const int ndt,
  const bool subcycledt)
{
  amrex::ParmParse ppre("react_eval");

  // Reference: CVODE with a direct solver (GMRES on GPUs, which have no
  // dense solver) at tight tolerances, Array4 layout. The ParmParse entries
  // added here are read by the reference reactor's init only: the reactor
  // under test has already read its own.
  amrex::Real ref_rtol = 1.e-13;
  amrex::Real ref_atol = 1.e-15;
#ifdef AMREX_USE_GPU
  std::string ref_solve_type = "GMRES";
#else
  std::string ref_solve_type = "dense_direct";
#endif
  ppre.query("ref_rtol", ref_rtol);
  ppre.query("ref_atol", ref_atol);
  ppre.query("ref_solve_type", ref_solve_type);
  // One call over dt: CVODE's setup per call dominates the cost otherwise
  int ref_ndt = 1;
  ppre.query("ref_ndt", ref_ndt);
  {
    amrex::ParmParse ppode("ode");
    ppode.add("rtol", ref_rtol);
    ppode.add("atol", ref_atol);
    ppode.add("verbose", 0);
    amrex::ParmParse ppcv("cvode");
    ppcv.add("solve_type", ref_solve_type);
  }
  amrex::Print() << "\n-- Reference: ReactorCvode, " << ref_solve_type
                 << ", rtol " << ref_rtol << ", atol " << ref_atol << ", "
                 << ref_ndt << " step(s) --\n";
  auto ref = pele::physics::reactions::ReactorBase::create("ReactorCvode");
  ref->init(reactor_type, 1);
  amrex::Vector<amrex::MultiFab> mf_ref = copyOf(mf0);
  amrex::Vector<amrex::MultiFab> mfE_ref = copyOf(mfE0);
  advance(
    *ref, finest_level, mf_ref, mfE_ref, rY_source_ext0, rY_source_energy_ext0,
    fctCount, dummyMask, 1, 1, dt, ref_ndt, subcycledt);

  const StateDiff d =
    compareStates(finest_level, mf, mf_ref, mf0, mfE0, reactor_type);

  // Tolerances on the largest difference over all cells: about 10x the
  // largest measured over drm19 (CH4) and LiDryer (H2), both layouts and
  // both reactor types, with the inputs files at ode.dt = 1e-5, ode.ndt =
  // 100, ode.rtol = ode.atol = 1e-10 (as in the CI); override with
  // react_eval.T_tol etc. CVODE is the least accurate at these settings:
  // with LiDryer the radicals of the induction period are below atol (on
  // rhoY), which shifts the ignition (3e-4 in T, 9e-6 at atol = 1e-14).
  // BDF takes fixed substeps (ode.bdf_nsubsteps) and conserves mass to its
  // nonlinear tolerance only.
  amrex::Real T_tol = 0.0;
  amrex::Real Y_tol = 0.0;
  amrex::Real rho_tol = 1.e-12;
  if (chem_integrator == "ReactorCvode") { // any solve_type
    T_tol = 3.e-3;                         // 3.1e-4 (LiDryer), 3.5e-6 (drm19)
    Y_tol = 1.2e-3;                        // 1.2e-4 (LiDryer), 1.3e-6 (drm19)
  } else if (chem_integrator == "ReactorArkode") {
    T_tol = 1.e-5; // 7.3e-7 (LiDryer)
    Y_tol = 3.e-6; // 2.3e-7 (LiDryer)
  } else if (chem_integrator == "ReactorRK64") {
    T_tol = 3.e-6; // 2.2e-7 (LiDryer)
    Y_tol = 1.e-6; // 6.4e-8 (LiDryer)
  } else if (chem_integrator == "ReactorBDF") {
    T_tol = 4.e-4;   // 3.3e-5 (drm19)
    Y_tol = 1.5e-4;  // 1.3e-5 (drm19)
    rho_tol = 5.e-7; // 3.7e-8 (LiDryer)
  }
  // The inversion for T stops at a Newton step < 1e-6 K
  amrex::Real Te_tol = 1.e-8; // 0 (T is computed from e or h)
  ppre.query("T_tol", T_tol);
  ppre.query("Y_tol", Y_tol);
  ppre.query("rho_tol", rho_tol);
  ppre.query("Te_tol", Te_tol);

  pele::testing::Checker check;
  amrex::Print() << "\n-- " << chem_integrator << " vs reference --\n";
  check.within("T, max relative difference", d.T, T_tol);
  check.within("Y, max absolute difference", d.Y, Y_tol);
  check.within(
    "rho = sum(rhoY) conserved, max relative change", d.rho, rho_tol);
  check.within(
    reactor_type == 1 ? "T consistent with the conserved e"
                      : "T consistent with the conserved h",
    d.Te, Te_tol);
  return check.summary("ReactEval");
}

} // namespace

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  int status = 0;
  {

    amrex::Real strt_time = amrex::ParallelDescriptor::second();

    BL_PROFILE_VAR("main::main()", pmain);

    amrex::ParmParse pp;
    std::string fuel_name;
    pp.get("fuel_name", fuel_name);

    std::string chem_integrator = "";
    pp.get("chem_integrator", chem_integrator);

    std::string pltfile;
    bool do_plt = false;
    if (pp.countval("plotfile") > 0) {
      pp.get("plotfile", pltfile);
      do_plt = true;
    }

    /* initial data */
    // Can either generate single level data
    // or read in data from a plotfile like file
    int initFromPlt = 0;
    pp.query("initFromFile", initFromPlt);
    std::string pltfileIn = "";
    if (initFromPlt) {
      pp.query("initFile", pltfileIn);
    }

    // With react_eval.check = 1 the final state is checked against a
    // reference integrated from a copy of the initial data
    amrex::ParmParse ppre("react_eval");
    bool do_check = false;
    ppre.query("check", do_check);
    if (do_check && initFromPlt != 0) {
      // The checks assume a closed reactor; plotfile data carries forcing
      amrex::Abort(
        "react_eval.check = 1 requires generated initial data (no forcing): "
        "it cannot be combined with initFromFile = 1");
    }

    bool subcycledt = true;
    pp.query("subcycle_dt", subcycledt);

    /* react() function version */
    // 1 -> Array4 version of react()  (Default)
    // 2 -> 1d raw pointer version of react()
    std::string reactFormat = "Array4";
    int reactFunc;
    pp.query("reactFormat", reactFormat);
    if (reactFormat == "Array4") {
      reactFunc = 1;
    } else if (reactFormat == "1dArray") {
      reactFunc = 2;
    } else {
      amrex::Abort(" --> reactFormat can only be 'Array4' or '1dArray' !");
    }

    /* ODE inputs */
    amrex::ParmParse ppode("ode");
    int ode_ncells = 1;
    ppode.query("ode_ncells", ode_ncells); // number of cells to integrate per
                                           // call, used only if reactFunc = 2

    amrex::Real dt = 1.e-5;
    ppode.query("dt", dt);

    int ndt = 1;
    ppode.query("ndt", ndt); // number of solver calls per dt

    int ode_iE = -1;
    ppode.query(
      "reactor_type",
      ode_iE); // RHS type, 1: e (PeleC), !1: h (PeleLM)  <------ FIXME!

    amrex::Real rtol = 1e-10;
    ppode.query("rtol", rtol);

    amrex::Real atol = 1e-10;
    ppode.query("atol", atol);

    int use_typ_vals = 0;
    ppode.query("use_typ_vals", use_typ_vals);

    amrex::Print() << "ODE solver: " << chem_integrator << std::endl;
    amrex::Print() << "Type of reactor: "
                   << (ode_iE == 1 ? "e (PeleC)" : "h (PeleLM)")
                   << std::endl; // <---- FIXME
    amrex::Print() << "Fuel: " << fuel_name << ", Oxy: O2" << std::endl;

    /* Mixture info */
    int fuel_idx = -1;
    if (fuel_name == "H2") {
      fuel_idx = H2_ID;
#ifdef CH4_ID
    } else if (fuel_name == "CH4") {
      fuel_idx = CH4_ID;
#endif
#ifdef NC12H26_ID
    } else if (fuel_name == "NC12H26") {
      fuel_idx = NC12H26_ID;
#endif
#ifdef IC8H18_ID
    } else if (fuel_name == "IC8H18") {
      fuel_idx = IC8H18_ID;
#endif
    }

    /* Initialize reactor object inside OMP region, including tolerances */
    BL_PROFILE_VAR("main::reactor_info()", reactInfo);
    std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
      pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(ode_iE, ode_ncells);
    BL_PROFILE_VAR_STOP(reactInfo);

    // -----------------------------------------------------------------------------
    // Initialize geom/data
    // When initFromPlt = 0, default is single level data provided by
    // initialize_data with a -1:1 unit length realbox in each dir and cell
    // count provided by the user otherwise number of levels, geom and all are
    // read in from the pltfileIn
    // -----------------------------------------------------------------------------
    int finest_level = 0;
    amrex::Vector<amrex::Geometry> geoms;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmaps;

    BL_PROFILE_VAR("main::geometry_setup", GeomSetup);
    if (initFromPlt == 0) {
      // -----------------------------------------------------------------------------
      // Resize vectors
      // -----------------------------------------------------------------------------
      geoms.resize(finest_level + 1);
      grids.resize(finest_level + 1);
      dmaps.resize(finest_level + 1);

      // -----------------------------------------------------------------------------
      // Define the geometry
      // -----------------------------------------------------------------------------
      std::array<int, 3> ncells = {AMREX_D_DECL(1, 1, 1)};
      if (pp.countval("ncells") == 1) {
        pp.get("ncells", ncells[0]);
        ncells = {AMREX_D_DECL(ncells[0], 1, 1)};
      } else if (pp.countval("ncells") >= AMREX_SPACEDIM) {
        amrex::Vector<int> nc(AMREX_SPACEDIM);
        pp.getarr("ncells", nc, 0, AMREX_SPACEDIM);
        ncells = {AMREX_D_DECL(nc[0], nc[1], nc[2])};
      } else {
        amrex::Abort("ncells has to have length 1 or spacedim");
      }

      amrex::Box domain(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(
          AMREX_D_DECL(ncells[0] - 1, ncells[1] - 1, ncells[2] - 1)));

      amrex::RealBox real_box(
        {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

      int coord = 0;

      amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

      geoms[0] = amrex::Geometry(domain, real_box, coord, is_periodic);

      amrex::Print() << "Integrating " << domain.numPts()
                     << " cells for: " << dt << " seconds with " << ndt
                     << " substeps \n";

      // -----------------------------------------------------------------------------
      // Define amrex::BoxArray / Dmap
      // -----------------------------------------------------------------------------
      int max_grid_size = 16;
      pp.query("max_grid_size", max_grid_size);
      grids[0] = amrex::BoxArray(domain);
      grids[0].maxSize(max_grid_size);
      dmaps[0] = amrex::DistributionMapping(
        grids[0], amrex::ParallelDescriptor::NProcs());
    } else {
      // Use PltFileManager, only metadata loaded here
      pele::physics::pltfilemanager::PltFileManager pltData(pltfileIn);

      if (ode_iE == 1) {
        amrex::Abort(
          "The option to read in chemical data is currently only available "
          "with PeleLM data and requires ode_iE=2");
      }

      //--------------------------------------------
      // General info
      //--------------------------------------------
      finest_level = pltData.getNlev() - 1;
      dt = pltData.getTime(); // In LMeX we stored the dt in the pltfileIn time

      // -----------------------------------------------------------------------------
      // Resize vectors
      // -----------------------------------------------------------------------------
      geoms.resize(finest_level + 1);
      grids.resize(finest_level + 1);
      dmaps.resize(finest_level + 1);

      // -----------------------------------------------------------------------------
      // Define geoms, read amrex::BoxArray and define dmap
      // -----------------------------------------------------------------------------
      geoms[0] = pltData.getGeom(0);
      for (int lev = 1; lev <= finest_level; ++lev) {
        geoms[lev] = amrex::refine(geoms[lev - 1], 2); // Assumes ref_ratio = 2
      }

      for (int lev = 0; lev <= finest_level; ++lev) {
        // Set vector entries
        grids[lev] = pltData.getGrid(lev);
        dmaps[lev] = amrex::DistributionMapping(
          grids[lev], amrex::ParallelDescriptor::NProcs());
      }

      for (int lev = 0; lev <= finest_level; ++lev) {
        if (subcycledt) {
          amrex::Print() << "  Level " << lev << " integrating "
                         << grids[lev].numPts() << " cells on "
                         << grids[lev].size() << " boxes for "
                         << dt / std::pow(2, lev) << " seconds \n";
        } else {
          amrex::Print() << "  Level " << lev << " integrating "
                         << grids[lev].numPts() << " cells on "
                         << grids[lev].size() << " boxes for " << dt
                         << " seconds \n";
        }
      }
    }
    BL_PROFILE_VAR_STOP(GeomSetup);

    // -----------------------------------------------------------------------------
    // Create MFs and generate initial data
    // -----------------------------------------------------------------------------
    BL_PROFILE_VAR("main::initialize_data()", InitData);
    int num_grow = 0;
    amrex::Vector<amrex::MultiFab> mf(finest_level + 1);
    amrex::Vector<amrex::MultiFab> rY_source_ext(finest_level + 1);
    amrex::Vector<amrex::MultiFab> mfE(finest_level + 1);
    amrex::Vector<amrex::MultiFab> rY_source_energy_ext(finest_level + 1);
    amrex::Vector<amrex::MultiFab> fctCount(finest_level + 1);
    amrex::Vector<amrex::iMultiFab> dummyMask(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      mf[lev].define(grids[lev], dmaps[lev], NUM_SPECIES + 1, num_grow);
      rY_source_ext[lev].define(grids[lev], dmaps[lev], NUM_SPECIES, num_grow);
      mfE[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      rY_source_energy_ext[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      fctCount[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      dummyMask[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      dummyMask[lev].setVal(1);
    }

    amrex::FabArrayBase::mfiter_tile_size =
      amrex::IntVect(AMREX_D_DECL(1024, 1024, 1024));

    // -----------------------------------------------------------------------------
    // Initialize data
    // -----------------------------------------------------------------------------
    if (initFromPlt == 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        const auto geomdata = geoms[lev].data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mf[lev], amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

          const amrex::Box& box = mfi.tilebox();

          amrex::Array4<amrex::Real> const& rY_a = mf[lev].array(mfi);
          amrex::Array4<amrex::Real> const& rYs_a =
            rY_source_ext[lev].array(mfi);
          amrex::Array4<amrex::Real> const& E_a = mfE[lev].array(mfi);
          amrex::Array4<amrex::Real> const& rE_a =
            rY_source_energy_ext[lev].array(mfi);

          amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              initialize_data(
                i, j, k, fuel_idx, ode_iE, rY_a, rYs_a, E_a, rE_a, geomdata);
            });
        }
      }
    } else {
      pele::physics::pltfilemanager::PltFileManager pltData(pltfileIn);
      amrex::Vector<std::string> plt_vars = pltData.getVariableList();

      // Figure out the indices of reauired entries
      // Assuming PeleLMeX evaluate mode was used: rhoY + rhoH + T + F_rhoY +
      // F_rhoH
      amrex::Vector<std::string> spec_names;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
        spec_names);
      const std::string firstSpecVar = "rhoY(" + spec_names[0] + ")";
      const std::string firstSpecForceVar = "F_rhoY(" + spec_names[0] + ")";
      int idFirstSpec = -1, idFirstSpecForce = -1, idrhoH = -1, idT = -1,
          idrhoHForce = -1;
      for (int i = 0; i < plt_vars.size(); ++i) {
        if (plt_vars[i] == firstSpecVar) {
          idFirstSpec = i;
        } else if (plt_vars[i] == firstSpecForceVar) {
          idFirstSpecForce = i;
        } else if (plt_vars[i] == "rhoH") {
          idrhoH = i;
        } else if (plt_vars[i] == "Temp") {
          idT = i;
        } else if (plt_vars[i] == "F_rhoH") {
          idrhoHForce = i;
        }
      }
      if (idFirstSpec * idFirstSpecForce * idrhoH * idT * idrhoHForce < 0) {
        amrex::Abort(
          "All of 'rhoY(<spec>)', 'F_rhoY(<spec>)', 'rhoH', 'Temp' "
          "and 'F_rhoH' must be in the pltfile !");
      }

      // Load the data from pltfile
      for (int lev = 0; lev <= finest_level; ++lev) {
        // rhoYs
        pltData.fillPatchFromPlt(
          lev, geoms[lev], idFirstSpec, 0, NUM_SPECIES, mf[lev]);

        // Temperature
        pltData.fillPatchFromPlt(lev, geoms[lev], idT, NUM_SPECIES, 1, mf[lev]);

        // RhoH
        pltData.fillPatchFromPlt(lev, geoms[lev], idrhoH, 0, 1, mfE[lev]);

        // Force rhoYs
        pltData.fillPatchFromPlt(
          lev, geoms[lev], idFirstSpecForce, 0, NUM_SPECIES,
          rY_source_ext[lev]);

        // Force rhoH
        pltData.fillPatchFromPlt(
          lev, geoms[lev], idrhoHForce, 0, 1, rY_source_energy_ext[lev]);

        // Convert from MKS -> CGS since we have PeleLMeX data
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mf[lev], amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
          const amrex::Box& box = mfi.tilebox();

          auto const& rhoYs = mf[lev].array(mfi);
          auto const& temp = mf[lev].array(mfi, NUM_SPECIES);
          auto const& rhoH = mfE[lev].array(mfi);
          auto const& force_Y = rY_source_ext[lev].array(mfi);
          auto const& force_E = rY_source_energy_ext[lev].array(mfi);
          amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < NUM_SPECIES; n++) {
                rhoYs(i, j, k, n) *= 1.0e-3;
                force_Y(i, j, k, n) *= 1.0e-3;
              }
              rhoH(i, j, k) *= 10.0;
              force_E(i, j, k) *= 10.0;
            });
        }
      }
    }
    BL_PROFILE_VAR_STOP(InitData);

    // Print initial data if needed
    BL_PROFILE_VAR_NS("PlotFile", PlotFile);
    if (do_plt) {
      BL_PROFILE_VAR_START(PlotFile);
      std::string outfile = amrex::Concatenate(pltfile, 0);
      amrex::Vector<int> isteps(finest_level + 1, 0);
      amrex::Vector<amrex::IntVect> refRatios(
        finest_level, {AMREX_D_DECL(2, 2, 2)});
      amrex::Vector<std::string> plt_VarsName;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        plt_VarsName.push_back("SPEC" + std::to_string(k));
      }
      plt_VarsName.push_back("TEMP");

      amrex::WriteMultiLevelPlotfile(
        outfile, finest_level + 1, GetVecOfConstPtrs(mf), plt_VarsName, geoms,
        0.0, isteps, refRatios);
      BL_PROFILE_VAR_STOP(PlotFile);
    }

    // -----------------------------------------------------------------------------
    // Set typical values
    // -----------------------------------------------------------------------------
    if (
      chem_integrator == "ReactorCvode" || chem_integrator == "ReactorArkode") {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      {
        if (use_typ_vals) {
          amrex::Vector<amrex::Real> typ_vals(NUM_SPECIES + 1, 1.0e-10);
          if (ppode.contains("typ_vals")) {
            amrex::Print()
              << "Using user-defined typical values for the absolute "
                 "tolerances of the ode solver.\n";
            ppode.getarr("typ_vals", typ_vals, 0, NUM_SPECIES + 1);
            for (int i = 0; i < NUM_SPECIES; ++i) {
              typ_vals[i] = amrex::max<amrex::Real>(typ_vals[i], 1.e-10);
            }
          } else {
            amrex::Print()
              << "Using typical values from the initial data for the "
                 "absolute tolerances of the ode solver.\n";
            for (int lev = 0; lev <= finest_level; ++lev) {
              /*
                Pretty sure TypVal should be rhoYs in CGS, but keep that around
                just in case. amrex::MultiFab
                massFrac(grids[lev],dmaps[lev],NUM_SPECIES,0); #ifdef
                AMREX_USE_OMP #pragma omp parallel if
                (amrex::Gpu::notInLaunchRegion()) #endif for (amrex::MFIter
                mfi(mf[lev],amrex::TilingIfNotGPU()); mfi.isValid();
                ++mfi) { const amrex::Box& box = mfi.tilebox(); auto const&
                rhoYs = mf[lev].const_array(mfi); auto const& Ys    =
                massFrac.array(mfi); amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                amrex::Real rho = 0.0;
                for (int n = 0; n < NUM_SPECIES; n++) {
                rho += rhoYs(i,j,k,n);
                }
                amrex::Real rhoinv = 1.0/rho;
                for (int n = 0; n < NUM_SPECIES; n++) {
                Ys(i,j,k,n) = rhoYs(i,j,k,n) * rhoinv;
                }
                });
                }
              */
              for (int i = 0; i < NUM_SPECIES; ++i) {
                typ_vals[i] =
                  amrex::max<amrex::Real>(typ_vals[i], mf[lev].max(i));
              }
              typ_vals[NUM_SPECIES] = amrex::max<amrex::Real>(
                typ_vals[NUM_SPECIES], mf[lev].max(NUM_SPECIES));
            }
          }
          reactor->set_typ_vals_ode(typ_vals);
        }
      }
    }
    amrex::Print() << " \n STARTING THE ADVANCE \n";

    amrex::Vector<amrex::MultiFab> mf0, mfE0, rY_source_ext0,
      rY_source_energy_ext0;
    if (do_check) {
      mf0 = copyOf(mf);
      mfE0 = copyOf(mfE);
      rY_source_ext0 = copyOf(rY_source_ext);
      rY_source_energy_ext0 = copyOf(rY_source_energy_ext);
    }

    /* REACT */
    advance(
      *reactor, finest_level, mf, mfE, rY_source_ext, rY_source_energy_ext,
      fctCount, dummyMask, reactFunc, ode_ncells, dt, ndt, subcycledt);

    if (do_check) {
      status = checkFinalState(
        chem_integrator, ode_iE, finest_level, mf, mf0, mfE0, rY_source_ext0,
        rY_source_energy_ext0, fctCount, dummyMask, dt, ndt, subcycledt);
    }

    // TODO multilevel max.
    {
      amrex::Vector<double> typ_vals(NUM_SPECIES + 1);
      amrex::Print() << "ode.typ_vals= ";
      for (int i = 0; i < NUM_SPECIES + 1; ++i) {
        amrex::Print() << amrex::max<amrex::Real>(1.e-10, mf[0].max(i)) << " ";
      }
      amrex::Print() << std::endl;
    }

    if (do_plt) {
      BL_PROFILE_VAR_START(PlotFile);
      std::string outfile = amrex::Concatenate(pltfile, 1);
      // TODO: add fct count to this output
      amrex::Vector<int> isteps(finest_level + 1, 0);
      amrex::Vector<amrex::IntVect> refRatios(
        finest_level, {AMREX_D_DECL(2, 2, 2)});
      amrex::Vector<std::string> plt_VarsName;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        plt_VarsName.push_back("SPEC" + std::to_string(k));
      }
      plt_VarsName.push_back("TEMP");

      amrex::WriteMultiLevelPlotfile(
        outfile, finest_level + 1, GetVecOfConstPtrs(mf), plt_VarsName, geoms,
        0.0, isteps, refRatios);
      BL_PROFILE_VAR_STOP(PlotFile);
    }

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Real run_time = amrex::ParallelDescriptor::second() - strt_time;
    amrex::ParallelDescriptor::ReduceRealMax(
      run_time, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << " \n >> React::main() " << run_time << "\n\n";
  }

  amrex::Finalize();

  return status;
}
