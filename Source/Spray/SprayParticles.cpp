
#include "SprayParticles.H"
#include "Drag.H"
#include "SprayInterpolation.H"
#include "Transport.H"
#include "WallFunctions.H"
#include "TABBreakup.H"
#include "ReitzKHRT.H"
#include "WallFilm.H"
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

void
SprayParticleContainer::init_bcs(const amrex::BCRec& phys_bc)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (
      phys_bc.lo(dir) == amrex::PhysBCType::symmetry ||
      phys_bc.lo(dir) == amrex::PhysBCType::slipwall ||
      phys_bc.lo(dir) == amrex::PhysBCType::noslipwall) {
      reflect_lo[dir] = true;
    } else {
      reflect_lo[dir] = false;
    }
    if (
      phys_bc.hi(dir) == amrex::PhysBCType::symmetry ||
      phys_bc.hi(dir) == amrex::PhysBCType::slipwall ||
      phys_bc.hi(dir) == amrex::PhysBCType::noslipwall) {
      reflect_hi[dir] = true;
    } else {
      reflect_hi[dir] = false;
    }
  }
}

void
SprayParticleContainer::moveKick(
  MultiFab& state,
  MultiFab& source,
  const int level,
  const Real& dt,
  const Real time,
  const bool isVirtualPart,
  const bool isGhostPart,
  const int state_ghosts,
  const int source_ghosts,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev)
{
  bool do_move = false;
  moveKickDrift(
    state, source, level, dt, time, isVirtualPart, isGhostPart, state_ghosts,
    source_ghosts, do_move, ltransparm, spray_cfl_lev);
}

void
SprayParticleContainer::moveKickDrift(
  MultiFab& state,
  MultiFab& source,
  const int level,
  const Real& dt,
  const Real time,
  const bool isVirtualPart,
  const bool isGhostPart,
  const int state_ghosts,
  const int source_ghosts,
  const bool do_move,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev,
  const MultiFab* state_end,
  const bool rate_only,
  MultiFab* source_end)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(level >= 0);

  // If there are no particles at this level
  if (level >= this->GetParticles().size()) {
    return;
  }

  updateParticles(
    level, state, source, dt, time, state_ghosts, source_ghosts, isVirtualPart,
    isGhostPart, do_move, ltransparm, spray_cfl_lev, state_end, rate_only,
    source_end);
}

Real
SprayParticleContainer::estTimestep(int level) const
{
  BL_PROFILE("ParticleContainer::estTimestep()");
  Real dt = std::numeric_limits<Real>::max();
  if (level >= this->GetParticles().size() || m_sprayData->fixed_parts) {
    return -1.;
  }
  const Real cfl = SprayParticleContainer::spray_cfl;
  const auto dx = Geom(level).CellSizeArray();
  const auto dxi = Geom(level).InvCellSizeArray();
  {
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      for (MyParConstIter pti(*this, level); pti.isValid(); ++pti) {
        const AoS& pbox = pti.GetArrayOfStructs();
        const ParticleType* pstruct = pbox().data();
        const int n = pbox.numParticles();
        reduce_op.eval(
          n, reduce_data, [=] AMREX_GPU_DEVICE(const int i) -> ReduceTuple {
            const ParticleType& p = pstruct[i];
            if (p.id() > 0) {
              const Real max_mag_vdx = amrex::max<amrex::Real>(AMREX_D_DECL(
                std::abs(p.rdata(SprayComps::pstateVel)) * dxi[0],
                std::abs(p.rdata(SprayComps::pstateVel + 1)) * dxi[1],
                std::abs(p.rdata(SprayComps::pstateVel + 2)) * dxi[2]));
              Real dt_part = (max_mag_vdx > 0.) ? (cfl / max_mag_vdx) : 1.E50;
              return dt_part;
            }
            return -1.;
          });
      }
      ReduceTuple hv = reduce_data.value();
      Real ldt_cpu = amrex::get<0>(hv);
      dt = amrex::min<amrex::Real>(dt, ldt_cpu);
    }
  }
  // Check if the velocity of particles being injected is greater than existing
  // particle velocities
  for (const auto& sj : m_sprayJets) {
    if (sj->max_jet_vel() > 0.) {
      dt = amrex::min<amrex::Real>(dt, cfl * dx[0] / sj->max_jet_vel());
    }
  }
  ParallelDescriptor::ReduceRealMin(dt);
  return dt;
}

void
SprayParticleContainer::updateParticles(
  const int& level,
  MultiFab& state,
  MultiFab& source,
  const Real& flow_dt,
  const Real& /*time*/,
  const int state_ghosts,
  const int source_ghosts,
  const bool isVirt,
  const bool isGhost,
  const bool do_move,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev,
  const MultiFab* state_end,
  const bool rate_only,
  MultiFab* source_end)
{
  BL_PROFILE("SprayParticleContainer::updateParticles()");
  AMREX_ASSERT(OnSameGrids(level, state));
  AMREX_ALWAYS_ASSERT(!rate_only || (!do_move && state_end == nullptr));
  AMREX_ALWAYS_ASSERT(
    source_end == nullptr || (state_end != nullptr && advanceGivesEndRate()));
  const bool has_end = (state_end != nullptr);
  const bool has_rate_end = (source_end != nullptr);
  AMREX_ASSERT(OnSameGrids(level, source));
  bool isActive = !(isVirt || isGhost);
  bool do_splash = (m_sprayData->do_splash && isActive && do_move);
  bool do_breakup = (m_sprayData->do_breakup > 0) && !rate_only;
  Real B0 = m_khrtB0;
  Real B1 = m_khrtB1;
  Real C3 = m_khrtC3;
  Real max_ppp = m_maxNumPPP;
  if (m_verbose > 2 && ParallelDescriptor::IOProcessor()) {
    std::string move_string = "MK";
    if (do_move) {
      move_string = "MKD";
    }
    std::string part_type = "Active";
    if (isGhost) {
      part_type = "Ghost";
    } else if (isVirt) {
      part_type = "Virtual";
    }
    Print() << move_string << " on " << part_type << " particles on level "
            << level << std::endl;
  }
  const auto dxiarr = this->Geom(level).InvCellSizeArray();
  const auto dxarr = this->Geom(level).CellSizeArray();
  const auto ploarr = this->Geom(level).ProbLoArray();
  const auto phiarr = this->Geom(level).ProbHiArray();
  const RealVect dxi(AMREX_D_DECL(dxiarr[0], dxiarr[1], dxiarr[2]));
  const RealVect dx(AMREX_D_DECL(dxarr[0], dxarr[1], dxarr[2]));
  const RealVect plo(AMREX_D_DECL(ploarr[0], ploarr[1], ploarr[2]));
  const RealVect phi(AMREX_D_DECL(phiarr[0], phiarr[1], phiarr[2]));
  const auto domain = this->Geom(level).Domain();
#ifdef AMREX_USE_EB
  const auto& factory =
    dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
  const auto* cellcent = &(factory.getCentroid());
  const auto* bndrycent = &(factory.getBndryCent());
  const auto* bndrynorm = &(factory.getBndryNormal());
  const auto* volfrac = &(factory.getVolFrac());
#endif
  IntVect bndry_lo; // Designation for boundary types
  IntVect bndry_hi; // 0 - Periodic, 1 - Reflective, -1 - Non-reflective
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(level).isPeriodic(dir)) {
      if (reflect_lo[dir]) {
        bndry_lo[dir] = 1;
      } else {
        bndry_lo[dir] = -1;
      }
      if (reflect_hi[dir]) {
        bndry_hi[dir] = 1;
      } else {
        bndry_hi[dir] = -1;
      }
    } else {
      bndry_lo[dir] = 0;
      bndry_hi[dir] = 0;
    }
  }
  // Spray model in 2D assumes narrow domain in z-direction (Lz = dx)
  const Real vol = AMREX_D_PICK(
    dx[0] * dx[0] * dx[0], dx[0] * dx[1] * dx[0], dx[0] * dx[1] * dx[2]);
  const Real inv_vol = 1. / vol;
  // If particle subcycling is being done, determine the number of subcycles
  // Note: this is different than the AMR subcycling
  Real sub_cfl = 0.5; // CFL for each subcycle
  Real sub_dt = flow_dt;
  int num_iter = 1;
  if (do_move && spray_cfl_lev > sub_cfl) {
    num_iter = static_cast<int>(std::ceil(spray_cfl_lev / sub_cfl));
    sub_dt = flow_dt / static_cast<Real>(num_iter);
  }
  Real avg_inject_mass = 0.;
  if (isActive && m_sprayData->do_breakup == 2) {
    int numJets = static_cast<int>(m_sprayJets.size());
    for (int jindx = 0; jindx < numJets; ++jindx) {
      Real injDia = m_sprayJets[jindx]->get_avg_dia();
      Real injN = m_sprayJets[jindx]->num_ppp();
      avg_inject_mass +=
        injN * std::pow(injDia, 3) / static_cast<Real>(numJets);
    }
  }
  // Particle components indices
  SprayComps SPI = m_sprayIndx;
  // Start the ParIter, which loops over separate sets of particles in different
  // boxes
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
      const Box tile_box = pti.tilebox();
      const Box src_box = pti.growntilebox(source_ghosts);
      const Box state_box = pti.growntilebox(state_ghosts);
      bool at_bounds = tile_at_bndry(tile_box, bndry_lo, bndry_hi, domain);
      const int Np = pti.numParticles();
      if (Np == 0) {
        continue;
      }
      AoS& pbox = pti.GetArrayOfStructs();
      ParticleType* pstruct = pbox().data();
      const SprayData* fdat = d_sprayData;
      Array4<const Real> const& Tarr = state.array(pti, SPI.utempIndx);
      Array4<const Real> const& rhoYarr = state.array(pti, SPI.specIndx);
      Array4<const Real> const& rhoarr = state.array(pti, SPI.rhoIndx);
      Array4<const Real> const& momarr = state.array(pti, SPI.momIndx);
      Array4<const Real> const& engarr = state.array(pti, SPI.engIndx);
      // Gas state at the end of the step (only with state_end)
      Array4<const Real> Tarr_e, rhoYarr_e, rhoarr_e, momarr_e, engarr_e;
      if (has_end) {
        Tarr_e = state_end->const_array(pti, SPI.utempIndx);
        rhoYarr_e = state_end->const_array(pti, SPI.specIndx);
        rhoarr_e = state_end->const_array(pti, SPI.rhoIndx);
        momarr_e = state_end->const_array(pti, SPI.momIndx);
        engarr_e = state_end->const_array(pti, SPI.engIndx);
      }
      Array4<Real> const& rhoYSrcarr = source.array(pti, SPI.specSrcIndx);
      Array4<Real> const& rhoSrcarr = source.array(pti, SPI.rhoSrcIndx);
      Array4<Real> const& momSrcarr = source.array(pti, SPI.momSrcIndx);
      Array4<Real> const& engSrcarr = source.array(pti, SPI.engSrcIndx);
      // End-of-step exchange rate (only with source_end)
      Array4<Real> rhoYSrcarr_r, rhoSrcarr_r, momSrcarr_r, engSrcarr_r;
      if (has_rate_end) {
        rhoYSrcarr_r = source_end->array(pti, SPI.specSrcIndx);
        rhoSrcarr_r = source_end->array(pti, SPI.rhoSrcIndx);
        momSrcarr_r = source_end->array(pti, SPI.momSrcIndx);
        engSrcarr_r = source_end->array(pti, SPI.engSrcIndx);
      }
      bool eb_in_box = false;

#ifdef AMREX_USE_EB
      eb_in_box = true;
      const auto& interp_fab = static_cast<EBFArrayBox const&>(state[pti]);
      const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();
      Array4<const Real> ccent_fab;
      Array4<const Real> bcent_fab;
      Array4<const Real> bnorm_fab;
      Array4<const Real> volfrac_fab;
      const auto& flags_array = flags.array();
      if (flags.getType(state_box) == FabType::regular) {
        eb_in_box = false;
      } else {
        // Cell centroids
        ccent_fab = cellcent->array(pti);
        // Centroid of EB
        bcent_fab = bndrycent->array(pti);
        // Normal of EB
        bnorm_fab = bndrynorm->array(pti);
        volfrac_fab = volfrac->array(pti);
      }
#endif
      bool do_splash_box = (do_splash && (eb_in_box || at_bounds));
      FArrayBox wf_fab;
      Array4<Real> wf_arr;
      if (do_splash_box) {
        wf_fab.resize(src_box, 1, The_Async_Arena());
        wf_fab.setVal<RunOn::Device>(0.);
        wf_arr = wf_fab.array();
        // TODO: Adjust this for EB faces
        Real face_area = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
        amrex::ParallelFor(Np, [=] AMREX_GPU_DEVICE(int pid) noexcept {
          ParticleType& p = pstruct[pid];
          if (p.id() > 0 && p.rdata(SprayComps::pstateFilmHght) > 0.) {
            fillFilmFab(wf_arr, p, face_area, plo, dx);
          }
        });
      }
      // Data structures for creating new particles during splashing/breakup
      Gpu::HostVector<splash_breakup> N_SB_h;
      Gpu::DeviceVector<splash_breakup> N_SB_d;
      SBVects refv;
      SBPtrs rf_d;
      bool make_new_drops =
        ((do_breakup || do_splash_box) && isActive && do_move);
      if (make_new_drops) {
        N_SB_h.assign(Np, splash_breakup::no_change);
        N_SB_d.resize(Np);
        Gpu::copyAsync(
          Gpu::hostToDevice, N_SB_h.begin(), N_SB_h.end(), N_SB_d.begin());
        refv.build(Np);
        refv.fillPtrs_d(rf_d);
      }
      auto* N_SB = N_SB_d.dataPtr();
      amrex::ParallelFor(Np, [=] AMREX_GPU_DEVICE(int pid) noexcept {
        ParticleType& p = pstruct[pid];
        if (p.id() > 0) {
          auto eos = pele::physics::PhysicsType::eos(fdat->eosparm);
          GasPhaseVals gpv;
          GpuArray<Real, SPRAY_FUEL_NUM>
            cBoilT; // Boiling temperature at current pressure
          eos.molecular_weight(gpv.mw.data());
          for (int n = 0; n < NUM_CHEM_SPECIES; ++n) {
            gpv.mw[n] *= SprayUnits::mass_conv;
          }
          GpuArray<IntVect, AMREX_D_PICK(2, 4, 8)>
            indx_array; // array of adjacent cells
          GpuArray<Real, AMREX_D_PICK(2, 4, 8)>
            weights; // array of corresponding weights
          RealVect lx = (p.pos() - plo) * dxi + 0.5;
          IntVect ijk = lx.floor(); // Upper cell center
          RealVect lxc = (p.pos() - plo) * dxi;
          IntVect ijkc = lxc.floor(); // Cell with particle
          IntVect bflags(IntVect::TheZeroVector());
          if (at_bounds) {
            // Check if particle has left the domain or is boundary adjacent
            bool left_dom =
              check_bounds(p.pos(), plo, phi, dx, bndry_lo, bndry_hi, bflags);
            if (left_dom) {
              Abort("Particle has incorrectly left the domain");
            }
          }
          // Used for ETAB breakup model
          Real Utan_total = 0.;
          Real Reyn_d = 0.;
          // Subcycle loop
          for (int cur_iter = 0; cur_iter < num_iter && p.id() > 0;
               ++cur_iter) {
            bool is_film = false;
            // Gather wall film values
            if (p.rdata(SprayComps::pstateFilmHght) > 0.) {
              is_film = true;
            }
            // Flag for whether we are near EB boundaries
            bool do_fe_interp = false;
#ifdef AMREX_USE_EB
            if (eb_in_box) {
              do_fe_interp = eb_interp(
                p, ijkc, ijk, dx, dxi, lx, plo, bflags, flags_array, ccent_fab,
                bcent_fab, bnorm_fab, volfrac_fab, fdat->min_eb_vfrac,
                indx_array.data(), weights.data());
            } else
#endif
            {
              trilinear_interp(
                ijk, lx, indx_array.data(), weights.data(), bflags);
            }
            if (is_film && rate_only) {
              // No instantaneous rate for wall film particles
              continue;
            }
            // spray_SDC (has_end or rate_only): second order along the
            // particle path, away from walls, EB and the domain boundaries.
            // path_weights: interpolation weights at pos, if pos is away from
            // them too.
            const bool sdc_path = (has_end || rate_only) && !do_fe_interp &&
                                  bflags == IntVect::TheZeroVector();
            auto path_weights =
              [&](const RealVect& pos, IntVect* idx, Real* wts) -> bool {
              if (at_bounds) {
                IntVect bf(IntVect::TheZeroVector());
                if (
                  check_bounds(pos, plo, phi, dx, bndry_lo, bndry_hi, bf) ||
                  bf != IntVect::TheZeroVector()) {
                  return false;
                }
              }
              const RealVect lpos = (pos - plo) * dxi + 0.5;
              trilinear_interp(
                lpos.floor(), lpos, idx, wts, IntVect::TheZeroVector());
#ifdef AMREX_USE_EB
              if (eb_in_box) {
                for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
                  if (!flags_array(idx[aindx]).isRegular()) {
                    return false;
                  }
                }
              }
#endif
              return true;
            };
            const RealVect vel0(AMREX_D_DECL(
              p.rdata(SprayComps::pstateVel),
              p.rdata(SprayComps::pstateVel + 1),
              p.rdata(SprayComps::pstateVel + 2)));
            const bool moves =
              has_end && do_move && !fdat->fixed_parts && !is_film;
            // Interpolate fluid state
            gpv.reset();
            InterpolateGasPhase(
              gpv, state_box, rhoarr, rhoYarr, Tarr, momarr, engarr,
              indx_array.data(), weights.data(), fdat->eosparm, fdat);
            // With has_end, the gas varies from theta = 0 (state) to 1
            // (state_end) over the step; the subcycle starts at theta0 and
            // ends at theta1, where gpv_1 is the gas at the particle position
            // predicted by Euler
            const Real theta0 =
              static_cast<Real>(cur_iter) / static_cast<Real>(num_iter);
            const Real theta1 =
              static_cast<Real>(cur_iter + 1) / static_cast<Real>(num_iter);
            GasPhaseVals gpv_1;
            // End-of-step exchange rate, from the last subcycle
            GasPhaseVals gpv_r;
            gpv_r.reset();
            if (has_end) {
              GasPhaseVals gpv_e;
              gpv_e.reset();
              InterpolateGasPhase(
                gpv_e, state_box, rhoarr_e, rhoYarr_e, Tarr_e, momarr_e,
                engarr_e, indx_array.data(), weights.data(), fdat->eosparm,
                fdat);
              const GasPhaseVals gpv_s = gpv;
              // Wall film: gas held at the middle of the subcycle
              gpv.setFluidState(
                gpv_s, gpv_e, is_film ? 0.5 * (theta0 + theta1) : theta0);
              GpuArray<IntVect, AMREX_D_PICK(2, 4, 8)> indx_1;
              GpuArray<Real, AMREX_D_PICK(2, 4, 8)> weights_1;
              if (
                moves && sdc_path &&
                path_weights(
                  p.pos() + sub_dt * vel0, indx_1.data(), weights_1.data())) {
                GasPhaseVals gpv_1s;
                gpv_1s.reset();
                InterpolateGasPhase(
                  gpv_1s, state_box, rhoarr, rhoYarr, Tarr, momarr, engarr,
                  indx_1.data(), weights_1.data(), fdat->eosparm, fdat);
                gpv_e.reset();
                InterpolateGasPhase(
                  gpv_e, state_box, rhoarr_e, rhoYarr_e, Tarr_e, momarr_e,
                  engarr_e, indx_1.data(), weights_1.data(), fdat->eosparm,
                  fdat);
                gpv_1.setFluidState(gpv_1s, gpv_e, theta1);
              } else {
                gpv_1.setFluidState(gpv_s, gpv_e, theta1);
              }
            }
            // Solve for avg mw and pressure at droplet location
            fdat->calcBoilT(gpv, cBoilT.data());
            if (is_film) {
              calculateFilmSource(
                sub_dt, gpv, *fdat, p, cBoilT.data(), ltransparm);
            } else {
              Reyn_d = calculateSpraySource(
                sub_dt, gpv, *fdat, p, cBoilT.data(), ltransparm,
                has_end ? &gpv_1 : nullptr, rate_only,
                (has_rate_end && cur_iter == num_iter - 1) ? &gpv_r : nullptr);
            }
            // Trapezoid in time for the position with has_end
            RealVect vel_move(AMREX_D_DECL(
              p.rdata(SprayComps::pstateVel),
              p.rdata(SprayComps::pstateVel + 1),
              p.rdata(SprayComps::pstateVel + 2)));
            if (has_end) {
              vel_move = 0.5 * (vel0 + vel_move);
            }
            if (p.id() > 0 && do_breakup) {
              // Update breakup variables and determine if breakup occurs
              if (fdat->do_breakup == 1) {
                Utan_total += updateBreakupTAB(
                  Reyn_d, sub_dt, cBoilT.data(), gpv, *fdat, p);
              }
              if (cur_iter == num_iter - 1) {
                if (fdat->do_breakup == 1 && make_new_drops) {
                  // Determine if parcel must be split into multiple parcels
                  splitDropletTAB(pid, p, max_ppp, N_SB, rf_d, Utan_total);
                } else {
                  // Update breakup for KH-RT model
                  updateBreakupKHRT(
                    pid, p, Reyn_d, fdat->dtmod * flow_dt, cBoilT.data(),
                    avg_inject_mass, B0, B1, C3, gpv, *fdat, N_SB, rf_d,
                    make_new_drops);
                }
              }
            }
            // Deposit the exchange. On the spray_SDC path it is spread over
            // the cells around the middle of the subcycle path with the
            // interpolation weights, which is second order in time for a
            // moving particle; otherwise all of it goes to the particle cell.
            auto deposit_to = [&](
                                Array4<Real> const& momA,
                                Array4<Real> const& rhoA,
                                Array4<Real> const& rhoYA,
                                Array4<Real> const& engA, const GasPhaseVals& g,
                                const IntVect& cur_indx, const Real wt) {
              Real cvol = inv_vol;
#ifdef AMREX_USE_EB
              if (flags_array(cur_indx).isSingleValued()) {
                cvol *= 1. / (volfrac_fab(cur_indx));
              }
#endif
              if (!src_box.contains(cur_indx)) {
                if (!isGhost) {
                  Abort(
                    "SprayParticleContainer::updateParticles() -- source box "
                    "too small");
                }
                return;
              }
              const Real cur_coef = -wt * cvol;
              if (fdat->mom_trans) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                  Gpu::Atomic::Add(
                    &momA(cur_indx, dir), cur_coef * g.fluid_mom_src[dir]);
                }
              }
              if (fdat->mass_trans) {
                Gpu::Atomic::Add(&rhoA(cur_indx), cur_coef * g.fluid_mass_src);
                for (int pc = 0; pc < fdat->N_pc; ++pc) {
                  const int pcspec = fdat->pc_indx[pc];
                  Gpu::Atomic::Add(
                    &rhoYA(cur_indx, pc), cur_coef * g.fluid_Y_dot[pcspec]);
                }
              }
              Gpu::Atomic::Add(&engA(cur_indx), cur_coef * g.fluid_eng_src);
            };
            // The step source: this subcycle's share of the step
            auto deposit = [&](const IntVect& cur_indx, const Real wt) {
              deposit_to(
                momSrcarr, rhoSrcarr, rhoYSrcarr, engSrcarr, gpv, cur_indx,
                wt * sub_dt / flow_dt);
            };
            GpuArray<IntVect, AMREX_D_PICK(2, 4, 8)> indx_d;
            GpuArray<Real, AMREX_D_PICK(2, 4, 8)> weights_d;
            const RealVect pos_mid =
              moves ? RealVect(p.pos() + 0.5 * sub_dt * vel_move) : p.pos();
            if (
              sdc_path &&
              path_weights(pos_mid, indx_d.data(), weights_d.data())) {
              for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
                if (weights_d[aindx] > 0.) {
                  deposit(indx_d[aindx], weights_d[aindx]);
                }
              }
            } else {
              deposit(ijkc, 1.);
            }
            // Real new_time = static_cast<Real>(cur_iter + 1) * sub_dt;
            // Modify particle position by whole time step
            if (do_move && !fdat->fixed_parts && p.id() > 0 && !is_film) {
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                p.pos(dir) += sub_dt * vel_move[dir];
              }
              if (at_bounds || do_fe_interp) {
                // First check if particle has exited the domain through a
                // Cartesian boundary
                bool left_dom = check_bounds(
                  p.pos(), plo, phi, dx, bndry_lo, bndry_hi, bflags);
                if (left_dom) {
                  p.id() = -1;
                } else {
                  Real film_h = 0.;
                  if (do_splash_box) {
                    film_h = wf_arr(ijkc, 0);
                  }
                  // Next reflect particles off BC or EB walls if necessary
                  impose_wall(
                    do_splash_box, pid, p, *fdat, dx, plo, phi, bflags,
                    cBoilT.data(), gpv.p_fluid, eb_in_box,
#ifdef AMREX_USE_EB
                    flags_array, bcent_fab, bnorm_fab,
#endif
                    ijkc, N_SB, rf_d, film_h);
                }
              } // if (at_bounds || fe_interp)
              // Update indices
              lx = (p.pos() - plo) * dxi + 0.5;
              ijk = lx.floor();
              lxc = (p.pos() - plo) * dxi;
              ijkc = lxc.floor(); // New cell center
            }
            if (isGhost && !src_box.contains(ijkc)) {
              p.id() = -1;
            }
            if (has_rate_end && cur_iter == num_iter - 1 && p.id() > 0) {
              // End-of-step rate, spread around the end position
              if (path_weights(p.pos(), indx_d.data(), weights_d.data())) {
                for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
                  if (weights_d[aindx] > 0.) {
                    deposit_to(
                      momSrcarr_r, rhoSrcarr_r, rhoYSrcarr_r, engSrcarr_r,
                      gpv_r, indx_d[aindx], weights_d[aindx]);
                  }
                }
              } else {
                deposit_to(
                  momSrcarr_r, rhoSrcarr_r, rhoYSrcarr_r, engSrcarr_r, gpv_r,
                  ijkc, 1.);
              }
            }
          } // End of subcycle loop
        } // End of p.id() > 0 check
      }); // End of loop over particles
      if (make_new_drops) {
        Gpu::copy(
          Gpu::deviceToHost, N_SB_d.begin(), N_SB_d.end(), N_SB_h.begin());
        bool get_new_parts = false;
        for (int n = 0; n < Np; n++) {
          if (N_SB_h[n] != splash_breakup::no_change) {
            get_new_parts = true;
          }
        }
        if (get_new_parts) {
          refv.retrieve_data();
          SBPtrs rfh;
          refv.fillPtrs_h(rfh);
          CreateSBDroplets(Np, sub_dt, N_SB_h.data(), rfh, level);
        }
      }
      Gpu::streamSynchronize();
    } // for (int MyParIter pti..
  }
}
