#ifndef GPU_MISC_H
#define GPU_MISC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>
#include <TestCheck.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
initialize_data(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& mf,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& rho,
  amrex::Array4<amrex::Real> const& energy,
  amrex::GeometryData const& geomdata,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{
  const amrex::Real* plo = geomdata.ProbLo();
  const amrex::Real* phi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::Real dTemp = 5.0;
  amrex::Real dRho = 0.005;
  amrex::Real denergy = 1000000.0;
  amrex::Real y = plo[1] + (j + 0.5) * dx[1];
  amrex::Real x = plo[0] + (i + 0.5) * dx[0];
  amrex::Real pi = 3.1415926535897932;
  amrex::GpuArray<amrex::Real, 3> L;
  amrex::GpuArray<amrex::Real, 3> P;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_lo;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_hi;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> ei;

  for (int n = 0; n < 3; n++) {
    L[n] = phi[n] - plo[n];
    P[n] = L[n] / 4.0;
  }
#ifdef USE_MANIFOLD_EOS
  for (int n = 0; n < NUM_SPECIES - 1; n++) {
    Y_lo[n] = 0.0;
    Y_hi[n] = 1.0 / (NUM_SPECIES - 1);
  }
  Y_lo[NUM_SPECIES - 1] = 1.0;
  Y_hi[NUM_SPECIES - 1] = 1.0;
#else
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y_lo[n] = 0.0;
    Y_hi[n] = 1.0 / NUM_SPECIES;
  }
  Y_lo[0] = 1.0;
#endif

  // T, Yk, rho
#if (AMREX_SPACEDIM == 1)
  temp(i, j, k) = 500.0;
  rho(i, j, k) = 0.01;
#else
  temp(i, j, k) = 500.0 + dTemp * std::sin(2.0 * pi * (y - plo[1]) / P[1]);
  rho(i, j, k) = 0.01 + dRho * std::sin(2.0 * pi * (y - plo[1]) / P[1]);
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf(i, j, k, n) = Y_lo[n] + (Y_hi[n] - Y_lo[n]) * (x - plo[0]) / L[0];
  }

  // compute real E + perturb it for future tests
  auto eos = pele::physics::PhysicsType::eos(eosparm);
  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

#ifndef USE_MANIFOLD_EOS
  eos.RTY2E(rho(i, j, k), temp(i, j, k), mf_pt.arr, energy(i, j, k));
#endif
#if (AMREX_SPACEDIM != 1)
  energy(i, j, k) += denergy * std::sin(2.0 * pi * (y - plo[1]) / P[1]);
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cp(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<const amrex::Real> const& rho,
  amrex::Array4<amrex::Real> const& cp,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{

  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos(eosparm);
  eos.RTY2Cp(rho(i, j, k), temp(i, j, k), mf_pt.arr, cp(i, j, k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_wdot(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<const amrex::Real> const& rho,
  amrex::Array4<amrex::Real> const& wdot,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{

  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_pt, rxn_pt;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos(eosparm);
  eos.RTY2WDOT(rho(i, j, k), temp(i, j, k), mf_pt.arr, rxn_pt.arr);

  for (int n = 0; n < NUM_SPECIES; n++) {
    wdot(i, j, k, n) = rxn_pt[n];
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_T_from_EY(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<const amrex::Real> const& rho,
  amrex::Array4<const amrex::Real> const& energy,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{

  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos(eosparm);
  eos.REY2T(rho(i, j, k), energy(i, j, k), mf_pt.arr, temp(i, j, k));
}

// Thermodynamic consistency checks, one relative error per check (main.cpp
// has the list and the tolerances). T spans [T_lo, T_hi] along y and the
// composition is the one initialize_data set (it varies along x); each cell
// is evaluated at the pressures below and keeps its largest error. Checks
// that do not apply to the EOS stay zero.
namespace eoscheck {
enum { TET = 0, THT, PRP, RPT, REPT, REPP, HE, CVFD, CPFD, CPCV, NUM };
} // namespace eoscheck

// Relative error, with a non-finite one (e.g. a diverged inversion) mapped
// to +inf so that the max over the checks and cells keeps it
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
eos_rel_err(const amrex::Real a, const amrex::Real b) noexcept
{
  return pele::testing::finiteOrInf(pele::testing::relErr(a, b));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
eos_checks(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<amrex::Real> const& err,
  amrex::GeometryData const& geomdata,
  const amrex::Real T_lo,
  const amrex::Real T_hi,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{
  auto eos = pele::physics::PhysicsType::eos(eosparm);
  amrex::Real Y[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y[n] = mf(i, j, k, n);
  }
  for (int c = 0; c < eoscheck::NUM; c++) {
    err(i, j, k, c) = 0.0;
  }

#ifdef USE_MANIFOLD_EOS
  // The density comes from the table: p(rho(Y)) is the nominal pressure
  amrex::ignore_unused(geomdata, T_lo, T_hi);
  const amrex::Real T = 0.0; // not used by the manifold EOS
  const amrex::Real p = eosparm->Pnom_cgs;
  amrex::Real rho = 0.0, p2 = 0.0;
  eos.PYT2R(p, Y, T, rho);
  eos.RTY2P(rho, T, Y, p2);
  err(i, j, k, eoscheck::PRP) = eos_rel_err(p2, p);
#else
  const amrex::Real ny = geomdata.Domain().length(1);
  const amrex::Real T = T_lo + (T_hi - T_lo) * (j + 0.5) / ny;
  // Central differences: truncation ~(dT/T)^2, round-off ~eps*|e|/(cv*dT)
  const amrex::Real dT = 1.e-4 * T;
  const amrex::Real patm = pele::physics::Constants::PATM;
  const amrex::Real pressures[3] = {patm, 10.0 * patm, 50.0 * patm};
  for (const amrex::Real p : pressures) {
    amrex::Real cell_err[eoscheck::NUM] = {0.0};

    amrex::Real rho = 0.0, p2 = 0.0, T2 = 0.0;
    eos.PYT2R(p, Y, T, rho);
    eos.RTY2P(rho, T, Y, p2);
    cell_err[eoscheck::PRP] = eos_rel_err(p2, p);
    eos.RYP2T(rho, Y, p, T2);
    cell_err[eoscheck::RPT] = eos_rel_err(T2, T);

    // The inversions iterate from the T passed in: start 25% off
    amrex::Real e = 0.0;
    eos.RTY2E(rho, T, Y, e);
    T2 = 1.25 * T;
    eos.REY2T(rho, e, Y, T2);
    cell_err[eoscheck::TET] = eos_rel_err(T2, T);
    T2 = 1.25 * T;
    eos.RYET2P(rho, Y, e, T2, p2);
    cell_err[eoscheck::REPT] = eos_rel_err(T2, T);
    cell_err[eoscheck::REPP] = eos_rel_err(p2, p);

    amrex::Real h = 0.0;
    eos.RTY2H(rho, T, Y, h);
    cell_err[eoscheck::HE] = eos_rel_err(h - e, p / rho);
#ifndef USE_SRK_EOS
    // SRK has no h inversion (RHY2T aborts)
    T2 = 1.25 * T;
    eos.RHY2T(rho, h, Y, T2);
    cell_err[eoscheck::THT] = eos_rel_err(T2, T);
#endif

    amrex::Real cv = 0.0, cp = 0.0;
    eos.RTY2Cv(rho, T, Y, cv);
    eos.RTY2Cp(rho, T, Y, cp);

    amrex::Real ep = 0.0, em = 0.0;
    eos.RTY2E(rho, T + dT, Y, ep);
    eos.RTY2E(rho, T - dT, Y, em);
    cell_err[eoscheck::CVFD] = eos_rel_err((ep - em) / (2.0 * dT), cv);

    amrex::Real rhop = 0.0, rhom = 0.0, hp = 0.0, hm = 0.0;
    eos.PYT2R(p, Y, T + dT, rhop);
    eos.PYT2R(p, Y, T - dT, rhom);
    eos.RTY2H(rhop, T + dT, Y, hp);
    eos.RTY2H(rhom, T - dT, Y, hm);
    cell_err[eoscheck::CPFD] = eos_rel_err((hp - hm) / (2.0 * dT), cp);

    // Ideal gas only: SRK has no spec for it
    amrex::Real wbar = 0.0;
    eos.Y2WBAR(Y, wbar);
    cell_err[eoscheck::CPCV] =
      eos_rel_err(cp - cv, pele::physics::Constants::RU / wbar);

    for (int c = 0; c < eoscheck::NUM; c++) {
      err(i, j, k, c) = amrex::max(err(i, j, k, c), cell_err[c]);
    }
  }
#endif
}

#endif
