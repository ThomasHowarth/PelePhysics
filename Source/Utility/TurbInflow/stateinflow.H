#ifndef STATEINFLOW_H
#define STATEINFLOW_H

// State inflow: named scalar fields (e.g. temp, Y(N2)) carried by a
// TurbInflow file next to the three velocity components, declared by a
// SCALARS_V1 trailer in its HDR (see turbinflow.cpp).  The file, its plane
// chunk and its through-plane position are managed by TurbInflow;
// TurbInflow::scalar_sampler() returns the device-callable sampler below
// for one face and time.  The solver decides what the scalars mean.

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX.H>

#include <cmath>

#include "../MeshMap/MeshMapEvaluator.H"

// Feature macro: TurbInflow files may carry scalar fields and
// TurbInflow::scalar_sampler() / scalar_names() exist.
#define PELEPHYSICS_TURBINFLOW_HAS_SCALARS 1

namespace pele::physics::turbinflow {

//! Coordinate in a turbulence file of a position x (already scaled to file
//! units) along transverse direction \p idim: the identity for a file that
//! is uniform in physical position, the inverse of the file's map for a
//! mesh-mapped one (with the tile_periodic wrap applied first).  Shared by
//! the host-side velocity path and the device-side scalar sampler.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
turb_file_coordinate(
  const bool has_map,
  const pele::physics::MeshMapEvaluator& map,
  const amrex::Real xi_lo,
  const amrex::Real xi_hi,
  const bool tile_periodic,
  const int periodic,
  const amrex::Real pboxlo,
  const amrex::Real pboxsize,
  const amrex::Real dx,
  const amrex::Real dxinv,
  const int idim,
  const amrex::Real x_file_phys) noexcept
{
  if (!has_map) {
    return x_file_phys;
  }
  const amrex::Real x_lo = map.x_phys_from_xi(idim, xi_lo, xi_lo, xi_hi);
  const amrex::Real x_hi = map.x_phys_from_xi(idim, xi_hi, xi_lo, xi_hi);

  amrex::Real x = x_file_phys;
  if (tile_periodic && periodic != 0) {
    // Tiled: bring x into the file's physical period before inverting.
    // The inverse clamps to [x_lo, x_hi], which would otherwise collapse
    // every tiled copy onto the edge of the file.
    const amrex::Real period = x_hi - x_lo;
    x -= std::floor((x - x_lo) / period) * period;
  } else if (x < x_lo || x >= x_hi) {
    // Outside the file's physical extent.  A physically-uniform file leaves
    // such positions unfilled (the range test of the interpolation); keep
    // that behaviour rather than let the clamped inverse pin them to the
    // edge row, by returning a Xi coordinate that is outside the box
    // whatever turb_center shifted pboxlo to (|shift| <= pboxsize/2).
    return (x < x_lo) ? xi_lo - pboxsize : xi_hi + pboxsize;
  }

  amrex::Real xi = map.xi_from_x_phys(idim, x, xi_lo, xi_hi);

  // When the target's map and Xi grid coincide with the file's, xi lands on
  // a file cell centre up to the round-off of the inverse (log/atanh, or
  // the Newton tolerance for InteriorStretch).  Snap it so that the
  // interpolation weights collapse exactly and the same-map case
  // reproduces the file bit for bit.
  const amrex::Real idx = (xi - pboxlo) * dxinv + 0.5;
  const amrex::Real idx_round = std::round(idx);
  constexpr amrex::Real snap_tol = 1.0e-8;
  if (std::abs(idx - idx_round) < snap_tol) {
    xi = pboxlo + (idx_round - 0.5) * dx;
  }
  return xi;
}

/**
 * \brief Device-callable view of the scalar fields of the turbulence file
 *        feeding one face, at one time.
 *
 * Obtained from TurbInflow::scalar_sampler() right before use: it points at
 * the plane chunk currently loaded in memory and holds the through-plane
 * interpolation weights for the requested time, so it is invalidated by the
 * next call to add_turb() or scalar_sampler() for that face.  Trivially
 * copyable, so it can be captured by value in a boundary-fill functor.
 */
struct StateInflowSampler
{
  int nscal = 0; //!< 0: no scalar data on this face
  amrex::Array4<const amrex::Real> sd;
  int scomp = 0; //!< Component of the first scalar in sd

  amrex::GpuArray<int, 2> npboxcells = {{0, 0}};
  amrex::GpuArray<amrex::Real, 2> pboxlo = {{0.0, 0.0}};
  amrex::GpuArray<amrex::Real, 2> pboxsize = {{0.0, 0.0}};
  amrex::GpuArray<amrex::Real, 2> dx = {{0.0, 0.0}};
  amrex::GpuArray<amrex::Real, 2> dxinv = {{0.0, 0.0}};
  amrex::GpuArray<int, 2> periodicity = {{1, 1}};
  bool tile_periodic = false;
  bool extrap_nonperiodic = false;
  bool lininterp = false;
  amrex::Real scale_loc = 1.0;

  // Through-plane stencil: planes k0, k0+1, k0+2 with weights cz
  int k0 = 0;
  amrex::GpuArray<amrex::Real, 3> cz = {{0.0, 0.0, 0.0}};

  bool has_map = false;
  pele::physics::MeshMapEvaluator map{};
  amrex::GpuArray<amrex::Real, 2> map_xi_lo = {{0.0, 0.0}};
  amrex::GpuArray<amrex::Real, 2> map_xi_hi = {{0.0, 0.0}};

  //! Transverse interpolation stencil of one position, shared by all the
  //! scalars sampled there
  struct Stencil
  {
    amrex::GpuArray<int, 2> i0{};
    amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 2> c{};
  };

  /**
   * \brief Locate the physical transverse position (\p x1, \p x2), in case
   *        units along the face's transverse directions
   *        TurbInflow::transverseDirs(dir, ...), in the file.
   *
   * Returns false when it lies outside the file (and the file is not tiled
   * there), so that the caller can fall back to its own boundary value.
   * Locate once, then interp() each scalar: the coordinate inversion and the
   * weights are the expensive part on a mapped file.
   */
  AMREX_GPU_HOST_DEVICE bool
  locate(const amrex::Real x1, const amrex::Real x2, Stencil& st) const noexcept
  {
    const amrex::GpuArray<amrex::Real, 2> xp = {x1, x2};
    for (int d = 0; d < 2; ++d) {
      const amrex::Real xf = turb_file_coordinate(
        has_map, map, map_xi_lo[d], map_xi_hi[d], tile_periodic, periodicity[d],
        pboxlo[d], pboxsize[d], dx[d], dxinv[d], d, xp[d] * scale_loc);
      const amrex::Real from_lo = xf - pboxlo[d];
      const bool inside = (from_lo >= 0.0 && from_lo < pboxsize[d]) ||
                          (tile_periodic && periodicity[d] != 0);
      if (!inside) {
        return false;
      }
      if (
        (from_lo < 0.5 * dx[d] || from_lo > pboxsize[d] - dx[d]) &&
        (periodicity[d] == 0) && !extrap_nonperiodic) {
        amrex::Error(
          "TurbInflow scalar interp stencil touches ghost cell for "
          "nonperiodic direction and extrap_nonperiodic option not used");
      }
      // From the first cell center (ghost cell, 1/2 dx outside of pbox)
      amrex::Real ff = from_lo * dxinv[d] + 0.5;
      int i0 = static_cast<int>(std::floor(ff));
      ff -= static_cast<amrex::Real>(i0);
      if (tile_periodic) {
        i0 = (i0 % npboxcells[d] + npboxcells[d]) % npboxcells[d];
      }
      st.i0[d] = i0;
      st.c[d][0] = lininterp ? 1.0 - ff : 0.5 * (ff - 1.0) * (ff - 2.0);
      st.c[d][1] = lininterp ? ff : ff * (2.0 - ff);
      st.c[d][2] = lininterp ? 0.0 : 0.5 * ff * (ff - 1.0);
    }
    return true;
  }

  //! Value of scalar \p s (0 <= s < nscal) at a located position
  [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Real
  interp(const Stencil& st, const int s) const noexcept
  {
    const int n = scomp + s;
    amrex::Real v = 0.0;
    for (int ii = 0; ii <= 2; ++ii) {
      for (int jj = 0; jj <= 2; ++jj) {
        const int i = st.i0[0] + ii;
        const int j = st.i0[1] + jj;
        const amrex::Real zdata = cz[0] * sd(i, j, k0, n) +
                                  cz[1] * sd(i, j, k0 + 1, n) +
                                  cz[2] * sd(i, j, k0 + 2, n);
        v += st.c[0][ii] * st.c[1][jj] * zdata;
      }
    }
    return v;
  }

  /**
   * \brief Value of scalar \p s at the physical transverse position
   *        (\p x1, \p x2): locate() then interp().
   *
   * Returns false, leaving \p val untouched, when the position lies outside
   * the file or \p s is out of range.
   */
  AMREX_GPU_HOST_DEVICE bool sample(
    const amrex::Real x1,
    const amrex::Real x2,
    const int s,
    amrex::Real& val) const noexcept
  {
    Stencil st;
    if (s < 0 || s >= nscal || !locate(x1, x2, st)) {
      return false;
    }
    val = interp(st, s);
    return true;
  }
};
} // namespace pele::physics::turbinflow

#endif
