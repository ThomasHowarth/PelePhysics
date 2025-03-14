#ifndef GPU_MISC_H
#define GPU_MISC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
initialize_data(
  int i,
  int j,
  int k,
  int fuel_id,
  int energy_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& frcExt,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Real t0,
  amrex::Real equiv_ratio,
  amrex::Real p,
  amrex::GeometryData const& geomdata,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{
  const amrex::Real* plo = geomdata.ProbLo();
  const amrex::Real* phi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::Real Temp_lo = 2000.0;
  amrex::Real Temp_hi = 2500.0;
  amrex::Real dTemp = 100.0;
  amrex::Real pressure = p;
  amrex::Real density = 0.0;
  amrex::Real energy = 0.0;
  amrex::Real temp = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> X;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y;

  amrex::Real y = plo[1] + (j + 0.5) * dx[1];
  // amrex::Real x = plo[0] + (i+0.5)*dx[0];
  amrex::Real pi = 3.1415926535897932;
  amrex::GpuArray<amrex::Real, 3> L;
  amrex::GpuArray<amrex::Real, 3> P;

  for (int n = 0; n < 3; n++) {
    L[n] = phi[n] - plo[n];
    P[n] = L[n] / 4.0;
  }
  // Y
  auto eos = pele::physics::PhysicsType::eos(eosparm);

#ifndef USE_MANIFOLD_EOS
  for (int n = 0; n < NUM_SPECIES; n++) {
    X[n] = 0.0;
  }
  // Assumes n-Dodecane
  X[O2_ID] = 1.0 / (equiv_ratio / 18.5 + 1 + 0.79 / 0.21);
  X[N2_ID] = 0.79 * X[O2_ID] / 0.21;
  X[fuel_id] = 1.0 - X[O2_ID] - X[N2_ID];
  eos.X2Y(&X[0], &Y[0]);
#else
  Y[0] = 0.0; // Progress variable
  Y[1] = 1.0; // Density
#endif

  // T
  temp = t0;
  // temp =  Temp_lo + (Temp_hi-Temp_lo)*y/L[1] + dTemp *
  // std::sin(2.0*pi*y/P[1]);
  if (energy_type == 1) {
    // get rho and E
    eos.PYT2R(pressure, &Y[0], temp, density);
    eos.RTY2E(density, temp, &Y[0], energy);
  } else {
    // get rho and H
    eos.PYT2R(pressure, &Y[0], temp, density);
    eos.TY2H(temp, &Y[0], energy);
  }
  // Fill vect
  for (int n = 0; n < NUM_SPECIES; n++) {
    rhoY(i, j, k, n) = Y[n] * density;
    frcExt(i, j, k, n) = 0.0;
  }
  rhoY(i, j, k, NUM_SPECIES) = temp;

  rhoE(i, j, k) = energy * density;

  frcEExt(i, j, k) = 0.0;
}

#endif
