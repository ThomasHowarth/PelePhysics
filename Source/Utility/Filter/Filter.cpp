#include <cmath>
#include <numeric>

#include <AMReX_GpuAsyncArray.H>
#include "Filter.H"

namespace {

// Ratios of the optimized 3 and 5 point approximations for fgr = 1..10,
// Table I in Sagaut & Grohens (1999) Int. J. Num. Meth. Fluids
constexpr int n_optimized = 10;
constexpr amrex::Real box_3pt_ratio[n_optimized] = {
  0.079, 0.274, 1.377, -2.375, -1.000, -0.779, -0.680, -0.627, -0.596, -0.575};
constexpr amrex::Real box_5pt_ratios[n_optimized][2] = {
  {0.0886, -0.0169}, {0.3178, -0.0130}, {1.0237, 0.0368},  {2.4414, 0.5559},
  {0.2949, 0.7096},  {-0.5276, 0.4437}, {-0.6708, 0.3302}, {-0.7003, 0.2767},
  {-0.7077, 0.2532}, {-0.6996, 0.2222}};
constexpr amrex::Real gaussian_3pt_ratio[n_optimized] = {
  0.0763, 0.2527, 1.1160, -3.144, -1.102,
  -0.809, -0.696, -0.638, -0.604, -0.581};
constexpr amrex::Real gaussian_5pt_ratios[n_optimized][2] = {
  {0.0871, -0.0175}, {0.2596, -0.0021}, {0.4740, 0.0785},  {0.1036, 0.2611},
  {-0.4252, 0.3007}, {-0.6134, 0.2696}, {-0.6679, 0.2419}, {-0.6836, 0.2231},
  {-0.6873, 0.2103}, {-0.6870, 0.2014}};

bool
has_optimized_ratios(const int fgr)
{
  return fgr >= 1 && fgr <= n_optimized;
}

// One 1D pass of the filter along the unit vector e
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
run_filter_1d(
  const amrex::IntVect& iv,
  const int n,
  const amrex::IntVect& e,
  const int ng,
  const amrex::Real* w,
  amrex::Array4<const amrex::Real> const& q)
{
  amrex::Real sum = 0.0;
  for (int m = -ng; m <= ng; m++) {
    sum += w[m + ng] * q(iv + e * m, n);
  }
  return sum;
}

// Applies the filter with 1D weights w (2 * ng + 1 of them, in device memory)
// to components [nstart, nstart + ncnt) of q on box, writing qh, as one 1D
// pass per direction. Each pass covers box grown by ng in the directions not
// yet filtered, so the last pass only needs box itself and the same input
// cells are read as by run_filter. A free function rather than a private
// member so that it can contain GPU lambdas under CUDA.
void
run_filter_separable(
  const amrex::Box& box,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& qh,
  const int nstart,
  const int ncnt,
  const int ng,
  const amrex::Real* w)
{
  AMREX_ALWAYS_ASSERT(
    nstart >= 0 && ncnt >= 0 && nstart + ncnt <= q.nComp() &&
    nstart + ncnt <= qh.nComp());

  // Intermediate results; the async arena frees them once the kernels are done
  amrex::Array<amrex::FArrayBox, 2> tmp;
  amrex::Array4<const amrex::Real> src = q;
  int src_comp = nstart;

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::Box pbox = box;
    for (int d = dir + 1; d < AMREX_SPACEDIM; ++d) {
      pbox.grow(d, ng);
    }

    amrex::Array4<amrex::Real> dst = qh;
    int dst_comp = nstart;
    if (dir < AMREX_SPACEDIM - 1) {
      tmp[dir % 2].resize(pbox, ncnt, amrex::The_Async_Arena());
      dst = tmp[dir % 2].array();
      dst_comp = 0;
    }

    const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
    amrex::ParallelFor(
      pbox, ncnt, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        amrex::ignore_unused(j, k); // unused below 3D
        dst(iv, n + dst_comp) = run_filter_1d(iv, n + src_comp, e, ng, w, src);
      });

    src = dst;
    src_comp = dst_comp;
  }
}

} // namespace

// Set the filter weights for the standard box filter
void
Filter::set_box_weights()
{
  AMREX_ASSERT(_fgr % 2 == 0 || _fgr == 1);
  resize_weights(_fgr / 2);

  for (auto& w : _weights) {
    w = 1.0 / _fgr;
  }

  // Only half the cell is used at the ends
  if (_fgr > 1) {
    _weights.front() = 0.5 * _weights.front();
    _weights.back() = _weights.front();
  }
}

// Set the filter weights for the standard gaussian filter
void
Filter::set_gaussian_weights()
{
  AMREX_ASSERT(_fgr % 2 == 0);
  resize_weights(_fgr / 2);

  amrex::Real gamma = 6.0;
  amrex::Real sigma = std::sqrt(1.0 / (2.0 * gamma)) * _fgr;
  for (int i = 0; i < static_cast<int>(_weights.size()); i++) {
    // equivalent to: std::sqrt(gamma / (PI * _fgr * _fgr)) * std::exp((-gamma
    // * (i-_ngrow) * (i-_ngrow)) / (_fgr * _fgr));
    _weights[i] =
      1.0 / (std::sqrt(2.0 * amrex::Math::pi<amrex::Real>()) * sigma) *
      std::exp((-(i - _ngrow) * (i - _ngrow)) / (2 * sigma * sigma));
  }
  // normalize to ensure it all sums to one
  amrex::Real sum = std::accumulate(_weights.begin(), _weights.end(), 0.0);
  for (auto& w : _weights) {
    w /= sum;
  }
}

// Set the filter weights for the 3pt polynomial truncation
// approximation of the box filter. See Eq. 26 in Sagaut & Grohens
// (1999) Int. J. Num. Meth. Fluids.
void
Filter::set_box_3pt_approx_weights()
{
  resize_weights(1);
  _weights[0] = _fgr * _fgr / 24.0;
  _weights[1] = (12.0 - _fgr * _fgr) / 12.0;
  _weights[2] = _weights[0];
}

// Set the filter weights for the 5pt polynomial truncation
// approximation of the box filter. See Eq. 27 in Sagaut & Grohens
// (1999) Int. J. Num. Meth. Fluids (though there are typos).
void
Filter::set_box_5pt_approx_weights()
{
  resize_weights(2);
  const int fgr2 = _fgr * _fgr;
  const int fgr4 = fgr2 * fgr2;
  _weights[0] = (3.0 * fgr4 - 20.0 * fgr2) / 5760.0;
  _weights[1] = (80.0 * fgr2 - 3.0 * fgr4) / 1440.0;
  _weights[2] = (3.0 * fgr4 - 100.0 * fgr2 + 960.0) / 960.0;
  _weights[3] = _weights[1];
  _weights[4] = _weights[0];
}

// Set the filter weights for the 5pt polynomial truncation
// approximation of the gaussian filter. See Eq. 29 in Sagaut &
// Grohens (1999) Int. J. Num. Meth. Fluids.
void
Filter::set_gaussian_5pt_approx_weights()
{
  resize_weights(2);
  const int fgr2 = _fgr * _fgr;
  const int fgr4 = fgr2 * fgr2;
  _weights[0] = (fgr4 - 4.0 * fgr2) / 1152.0;
  _weights[1] = (16.0 * fgr2 - fgr4) / 288.0;
  _weights[2] = (fgr4 - 20.0 * fgr2 + 192.0) / 192.0;
  _weights[3] = _weights[1];
  _weights[4] = _weights[0];
}

void
Filter::set_3pt_weights_from_ratio(const amrex::Real ratio)
{
  resize_weights(1);
  _weights[0] = ratio / (1 + 2.0 * ratio);
  _weights[1] = 1.0 - 2.0 * _weights[0];
  _weights[2] = _weights[0];
}

void
Filter::set_5pt_weights_from_ratios(
  const amrex::Real ratio1, const amrex::Real ratio2)
{
  resize_weights(2);
  _weights[0] = ratio2 / (1 + 2.0 * ratio1 + 2.0 * ratio2);
  _weights[1] = ratio1 / ratio2 * _weights[0];
  _weights[2] = 1.0 - 2.0 * _weights[0] - 2.0 * _weights[1];
  _weights[3] = _weights[1];
  _weights[4] = _weights[0];
}

// Optimized approximations (Table I in Sagaut & Grohens (1999)). Beyond the
// table, the box ones fall back to the standard box filter and the gaussian
// ones to the polynomial approximations.
void
Filter::set_box_3pt_optimized_approx_weights()
{
  if (has_optimized_ratios(_fgr)) {
    set_3pt_weights_from_ratio(box_3pt_ratio[_fgr - 1]);
  } else {
    set_box_weights();
  }
}

void
Filter::set_box_5pt_optimized_approx_weights()
{
  if (has_optimized_ratios(_fgr)) {
    set_5pt_weights_from_ratios(
      box_5pt_ratios[_fgr - 1][0], box_5pt_ratios[_fgr - 1][1]);
  } else {
    set_box_weights();
  }
}

void
Filter::set_gaussian_3pt_optimized_approx_weights()
{
  if (has_optimized_ratios(_fgr)) {
    set_3pt_weights_from_ratio(gaussian_3pt_ratio[_fgr - 1]);
  } else {
    set_box_3pt_approx_weights(); // same as the 3pt gaussian approximation
  }
}

void
Filter::set_gaussian_5pt_optimized_approx_weights()
{
  if (has_optimized_ratios(_fgr)) {
    set_5pt_weights_from_ratios(
      gaussian_5pt_ratios[_fgr - 1][0], gaussian_5pt_ratios[_fgr - 1][1]);
  } else {
    set_gaussian_5pt_approx_weights();
  }
}

// Run the filtering operation on a MultiFab
void
Filter::apply_filter(const amrex::MultiFab& in, amrex::MultiFab& out) const
{
  apply_filter(in, out, 0, out.nComp());
}

void
Filter::apply_filter(
  const amrex::MultiFab& in,
  amrex::MultiFab& out,
  const int nstart,
  const int ncnt) const
{
  BL_PROFILE("Filter::apply_filter()");
  AMREX_ASSERT(in.nGrow() >= out.nGrow() + _ngrow);

  // The device copy of the weights is freed once the kernels are done
  amrex::Gpu::AsyncArray<amrex::Real> weights(_weights.data(), _weights.size());
  const amrex::Real* w = weights.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    run_filter_separable(
      mfi.growntilebox(), in.const_array(mfi), out.array(mfi), nstart, ncnt,
      _ngrow, w);
  }
}

// Run the filtering operation on a FAB
void
Filter::apply_filter(
  const amrex::Box& cbox,
  const amrex::FArrayBox& in,
  amrex::FArrayBox& out) const
{
  apply_filter(cbox, in, out, 0, out.nComp());
}

void
Filter::apply_filter(
  const amrex::Box& cbox,
  const amrex::FArrayBox& in,
  amrex::FArrayBox& out,
  const int nstart,
  const int ncnt) const
{
  BL_PROFILE("Filter::apply_filter()");
  amrex::Gpu::AsyncArray<amrex::Real> weights(_weights.data(), _weights.size());
  run_filter_separable(
    cbox, in.const_array(), out.array(), nstart, ncnt, _ngrow, weights.data());
}

void
Filter::apply_filter(
  const amrex::Box& box,
  const amrex::FArrayBox& in,
  amrex::FArrayBox& out,
  const int nstart,
  const int ncnt,
  const int /*ncomp*/) const
{
  apply_filter(box, in, out, nstart, ncnt);
}

void
Filter::apply_filter_1d(
  const int dir,
  const amrex::MultiFab& src,
  const int scomp,
  amrex::MultiFab& dst,
  const int dcomp,
  const int ncnt) const
{
  BL_PROFILE("Filter::apply_filter_1d()");
  AMREX_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM);
  AMREX_ASSERT(src.nGrowVect()[dir] >= _ngrow);

  // The device copy of the weights is freed once the kernels are done
  amrex::Gpu::AsyncArray<amrex::Real> weights(_weights.data(), _weights.size());
  const amrex::Real* w = weights.data();
  const int ng = _ngrow;
  const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

  const auto& s = src.const_arrays();
  const auto& d = dst.arrays();
  amrex::ParallelFor(
    dst, amrex::IntVect(0), ncnt,
    [=] AMREX_GPU_DEVICE(int b, int i, int j, int k, int n) noexcept {
      const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
      amrex::ignore_unused(j, k); // unused below 3D
      d[b](iv, n + dcomp) = run_filter_1d(iv, n + scomp, e, ng, w, s[b]);
    });
  amrex::Gpu::streamSynchronize();
}

void
apply_filter_low_memory(
  const amrex::Vector<Filter>& filters,
  const amrex::Vector<const amrex::MultiFab*>& in,
  const amrex::Vector<amrex::MultiFab*>& out,
  const int nstart,
  const int ncnt,
  const FilterFillFn& fill)
{
  BL_PROFILE("apply_filter_low_memory()");
  const int nlev = static_cast<int>(in.size());
  AMREX_ALWAYS_ASSERT(
    static_cast<int>(out.size()) == nlev &&
    static_cast<int>(filters.size()) == nlev);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    // Data filtered along the previous directions, with ghost cells along dir
    amrex::Vector<amrex::MultiFab> tmp(nlev);
    for (int lev = 0; lev < nlev; ++lev) {
      const amrex::MultiFab& src = (dir == 0) ? *in[lev] : *out[lev];
      amrex::IntVect ng(0);
      ng[dir] = filters[lev].get_filter_ngrow();
      tmp[lev].define(
        out[lev]->boxArray(), out[lev]->DistributionMap(), ncnt, ng);
      amrex::MultiFab::Copy(tmp[lev], src, nstart, 0, ncnt, 0);
    }
    fill(tmp, nstart);
    for (int lev = 0; lev < nlev; ++lev) {
      filters[lev].apply_filter_1d(dir, tmp[lev], 0, *out[lev], nstart, ncnt);
    }
  }
}
