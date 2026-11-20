#ifndef SURFACECHEMISTRY_H
#define SURFACECHEMISTRY_H

#include <AMReX_Algorithm.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <limits>
#include "mechanism.H"

// Surface chemistry at a catalytic wall, built on the routines CEPTR
// generates for heterogeneous mechanisms (comp_qfqr_surface,
// surfaceProductionRate, get_site_info).
//
// Inputs and outputs are in CGS, like the EOS: gas density in g/cm^3, gas
// species mass fluxes in g/cm^2/s. Coverages theta_k are site fractions of
// the NUM_SURFACE_SPECIES surface species and sum to one.
//
// Coverages can be handled in two ways:
// - transient: the caller stores the coverages and advances them with
//   coverageRates / surfaceRates;
// - quasi-steady: steadyCoverages finds the coverages at which the surface
//   species are at steady state for the given gas state.
//
// Everything is only defined for mechanisms with surface species.

namespace pele::physics::surface {

#if NUM_SURFACE_SPECIES > 0

namespace detail {

// Gas concentrations [mol/m^3] in sc[0 .. NUM_GAS_SPECIES-1].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gasConcentrations(
  const amrex::Real rho, const amrex::Real* Y, amrex::Real* sc) noexcept
{
  amrex::Real imw[NUM_GAS_SPECIES];
  get_imw(imw);
  for (int k = 0; k < NUM_GAS_SPECIES; ++k) {
    // mol/cm^3 -> mol/m^3
    sc[k] = 1.0e6 * rho * Y[k] * imw[k];
  }
}

// Surface concentrations [mol/m^2] in sc[NUM_GAS_SPECIES ..].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
surfaceConcentrations(
  const amrex::Real* theta,
  const amrex::Real site_density,
  const amrex::Real* site_size,
  amrex::Real* sc) noexcept
{
  for (int k = 0; k < NUM_SURFACE_SPECIES; ++k) {
    sc[NUM_GAS_SPECIES + k] = theta[k] * site_density / site_size[k];
  }
}

// dtheta/dt [1/s] for the gas concentrations already in sc.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
coverageRatesFromConc(
  const amrex::Real T,
  const amrex::Real* theta,
  const amrex::Real site_density,
  const amrex::Real* site_size,
  amrex::Real* sc,
  amrex::Real* dthetadt) noexcept
{
  surfaceConcentrations(theta, site_density, site_size, sc);
  amrex::Real sdot[NUM_ALL_SPECIES];
  surfaceProductionRate(sdot, sc, T);
  for (int k = 0; k < NUM_SURFACE_SPECIES; ++k) {
    dthetadt[k] = sdot[NUM_GAS_SPECIES + k] * site_size[k] / site_density;
  }
}

// Solve A x = b in place (b becomes x) by Gaussian elimination with partial
// pivoting; A is n x n, row major, and is overwritten. Returns false if A is
// singular.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
denseSolve(amrex::Real* A, amrex::Real* b) noexcept
{
  constexpr int n = NUM_SURFACE_SPECIES;
  for (int c = 0; c < n; ++c) {
    int p = c;
    for (int r = c + 1; r < n; ++r) {
      if (amrex::Math::abs(A[r * n + c]) > amrex::Math::abs(A[p * n + c])) {
        p = r;
      }
    }
    if (A[p * n + c] == 0.0) {
      return false;
    }
    if (p != c) {
      for (int j = 0; j < n; ++j) {
        amrex::Real tmp = A[c * n + j];
        A[c * n + j] = A[p * n + j];
        A[p * n + j] = tmp;
      }
      amrex::Real tmp = b[c];
      b[c] = b[p];
      b[p] = tmp;
    }
    for (int r = c + 1; r < n; ++r) {
      const amrex::Real f = A[r * n + c] / A[c * n + c];
      for (int j = c; j < n; ++j) {
        A[r * n + j] -= f * A[c * n + j];
      }
      b[r] -= f * b[c];
    }
  }
  for (int r = n - 1; r >= 0; --r) {
    for (int j = r + 1; j < n; ++j) {
      b[r] -= A[r * n + j] * b[j];
    }
    b[r] /= A[r * n + r];
  }
  return true;
}

// One backward Euler step theta = theta0 + dt * g(theta) by Newton's method
// with a finite-difference Jacobian, converged when no coverage changes by
// more than tol. Site conservation makes dg/dtheta singular along
// sum(theta), where round-off in g would otherwise accumulate, so the
// equation of the most covered species is replaced by sum(theta) = 1.
// theta holds theta0 on input. Returns false (theta unchanged) if Newton
// fails to converge or leaves a coverage negative.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
backwardEulerStep(
  const amrex::Real T,
  const amrex::Real dt,
  const amrex::Real tol,
  const amrex::Real site_density,
  const amrex::Real* site_size,
  amrex::Real* sc,
  amrex::Real* theta) noexcept
{
  constexpr int n = NUM_SURFACE_SPECIES;
  constexpr int max_newton = 20;
  amrex::Real th[n], g[n], gp[n], res[n], J[n * n];
  int kmax = 0;
  for (int k = 0; k < n; ++k) {
    th[k] = theta[k];
    if (theta[k] > theta[kmax]) {
      kmax = k;
    }
  }
  for (int it = 0; it < max_newton; ++it) {
    coverageRatesFromConc(T, th, site_density, site_size, sc, g);
    amrex::Real sum = 0.0;
    for (int k = 0; k < n; ++k) {
      res[k] = -(th[k] - theta[k] - dt * g[k]);
      sum += th[k];
    }
    res[kmax] = 1.0 - sum;
    // J = I - dt dg/dtheta
    for (int j = 0; j < n; ++j) {
      const amrex::Real h = 1.0e-8 * amrex::max<amrex::Real>(th[j], 1.0e-6);
      const amrex::Real saved = th[j];
      th[j] += h;
      coverageRatesFromConc(T, th, site_density, site_size, sc, gp);
      th[j] = saved;
      for (int k = 0; k < n; ++k) {
        J[k * n + j] = ((k == j) ? 1.0 : 0.0) - dt * (gp[k] - g[k]) / h;
      }
      J[kmax * n + j] = 1.0;
    }
    if (!denseSolve(J, res)) {
      return false;
    }
    amrex::Real dmax = 0.0;
    for (int k = 0; k < n; ++k) {
      const amrex::Real d = amrex::Math::abs(res[k]);
      // an update of a whole coverage, or NaN (which max() would drop),
      // means Newton has diverged
      if (!(d < 1.0)) {
        return false;
      }
      th[k] += res[k];
      dmax = amrex::max(dmax, d);
    }
    if (dmax < tol) {
      for (int k = 0; k < n; ++k) {
        if (th[k] < -1.0e-10) {
          return false;
        }
      }
      // remove round-off: clip and renormalize the site fractions
      amrex::Real sum = 0.0;
      for (int k = 0; k < n; ++k) {
        th[k] = amrex::max<amrex::Real>(th[k], 0.0);
        sum += th[k];
      }
      for (int k = 0; k < n; ++k) {
        theta[k] = th[k] / sum;
      }
      return true;
    }
  }
  return false;
}

} // namespace detail

/// Coverage rates of change dtheta/dt [1/s] of the surface species.
///
/// @param rho   gas density [g/cm^3]
/// @param T     temperature [K]
/// @param Y     gas mass fractions (NUM_GAS_SPECIES)
/// @param theta coverages (NUM_SURFACE_SPECIES)
/// @param dthetadt output (NUM_SURFACE_SPECIES)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
coverageRates(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real* Y,
  const amrex::Real* theta,
  amrex::Real* dthetadt) noexcept
{
  amrex::Real site_density, site_size[NUM_SURFACE_SPECIES];
  get_site_info(site_density, site_size);
  amrex::Real sc[NUM_ALL_SPECIES];
  detail::gasConcentrations(rho, Y, sc);
  detail::coverageRatesFromConc(
    T, theta, site_density, site_size, sc, dthetadt);
}

/// Gas species mass fluxes from the wall and coverage rates of change.
///
/// @param rho   gas density [g/cm^3]
/// @param T     temperature [K]
/// @param Y     gas mass fractions (NUM_GAS_SPECIES)
/// @param theta coverages (NUM_SURFACE_SPECIES)
/// @param flux  output: gas species mass production by the surface
///              [g/cm^2/s] (NUM_GAS_SPECIES), positive into the gas
/// @param dthetadt output: [1/s] (NUM_SURFACE_SPECIES)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
surfaceRates(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real* Y,
  const amrex::Real* theta,
  amrex::Real* flux,
  amrex::Real* dthetadt) noexcept
{
  amrex::Real site_density, site_size[NUM_SURFACE_SPECIES];
  get_site_info(site_density, site_size);
  amrex::Real sc[NUM_ALL_SPECIES];
  detail::gasConcentrations(rho, Y, sc);
  detail::surfaceConcentrations(theta, site_density, site_size, sc);
  amrex::Real sdot[NUM_ALL_SPECIES];
  surfaceProductionRate(sdot, sc, T);
  amrex::Real mw[NUM_GAS_SPECIES];
  get_mw(mw);
  for (int k = 0; k < NUM_GAS_SPECIES; ++k) {
    // mol/m^2/s * g/mol -> g/cm^2/s
    flux[k] = 1.0e-4 * sdot[k] * mw[k];
  }
  for (int k = 0; k < NUM_SURFACE_SPECIES; ++k) {
    dthetadt[k] = sdot[NUM_GAS_SPECIES + k] * site_size[k] / site_density;
  }
}

/// Quasi-steady coverages for a fixed gas state.
///
/// Integrates dtheta/dt with backward Euler steps of growing size until a
/// step of max_dt [s] changes no coverage by more than atol. theta is the
/// initial guess on input (a previous solution makes a good warm start) and
/// the steady coverages on output; it is left unchanged on failure.
///
/// @return true on convergence
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
steadyCoverages(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real* Y,
  amrex::Real* theta,
  const amrex::Real atol = 1.0e-12,
  const amrex::Real max_dt = 1.0e10,
  const int max_steps = 500) noexcept
{
  constexpr int n = NUM_SURFACE_SPECIES;
  amrex::Real site_density, site_size[n];
  get_site_info(site_density, site_size);
  amrex::Real sc[NUM_ALL_SPECIES];
  detail::gasConcentrations(rho, Y, sc);

  amrex::Real th[n], th_old[n], g[n];
  for (int k = 0; k < n; ++k) {
    th[k] = theta[k];
  }
  // first step changes the coverages by about 0.1
  detail::coverageRatesFromConc(T, th, site_density, site_size, sc, g);
  amrex::Real gmax = 0.0;
  for (int k = 0; k < n; ++k) {
    gmax = amrex::max(gmax, amrex::Math::abs(g[k]));
  }
  if (!(gmax < std::numeric_limits<amrex::Real>::max())) {
    return false;
  }
  amrex::Real dt =
    (gmax > 0.0) ? amrex::min<amrex::Real>(0.1 / gmax, max_dt) : max_dt;

  for (int step = 0; step < max_steps; ++step) {
    for (int k = 0; k < n; ++k) {
      th_old[k] = th[k];
    }
    if (!detail::backwardEulerStep(
          T, dt, atol, site_density, site_size, sc, th)) {
      dt *= 0.25;
      if (dt < 1.0e-30) {
        return false;
      }
      continue;
    }
    amrex::Real dmax = 0.0;
    for (int k = 0; k < n; ++k) {
      dmax = amrex::max(dmax, amrex::Math::abs(th[k] - th_old[k]));
    }
    if (dmax < atol) {
      if (dt >= max_dt) {
        for (int k = 0; k < n; ++k) {
          theta[k] = th[k];
        }
        return true;
      }
      // (nearly) steady: confirm with a step of max_dt
      dt = max_dt;
    } else {
      dt = amrex::min<amrex::Real>(3.0 * dt, max_dt);
    }
  }
  return false;
}

#endif

} // namespace pele::physics::surface

#endif
