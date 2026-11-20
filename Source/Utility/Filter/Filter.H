#ifndef FILTER_H
#define FILTER_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

// Filter types
enum filter_types {
  no_filter = 0,
  box,                           // 1
  gaussian,                      // 2
  box_3pt_approx,                // 3
  box_5pt_approx,                // 4
  box_3pt_optimized_approx,      // 5
  box_5pt_optimized_approx,      // 6
  gaussian_3pt_approx,           // 7
  gaussian_5pt_approx,           // 8
  gaussian_3pt_optimized_approx, // 9
  gaussian_5pt_optimized_approx, // 10
  num_filter_types
};

// Direct multi-dimensional convolution with the 1D weights w at one cell,
// accumulated into qh. Filter::apply_filter does not use it (it applies the
// same filter as one 1D pass per direction); it is kept as a reference.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
run_filter(
  const int i,
  const int j,
  const int k,
  const int nc,
  const int ng,
  const int nstart,
  const amrex::Real* w,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& qh)
{
  // Only filter along the directions that exist
  const int ngy = (AMREX_SPACEDIM > 1) ? ng : 0;
  const int ngz = (AMREX_SPACEDIM > 2) ? ng : 0;
  for (int n = -ngz; n <= ngz; n++) {
    const amrex::Real wz = (AMREX_SPACEDIM > 2) ? w[n + ng] : 1.0;
    for (int m = -ngy; m <= ngy; m++) {
      const amrex::Real wy = (AMREX_SPACEDIM > 1) ? w[m + ng] : 1.0;
      for (int l = -ng; l <= ng; l++) {
        qh(i, j, k, nc + nstart) +=
          w[l + ng] * wy * wz * q(i + l, j + m, k + n, nc + nstart);
      }
    }
  }
}

// One 1D pass of the filter along the direction (di, dj, dk)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
run_filter_1d(
  const int i,
  const int j,
  const int k,
  const int n,
  const int di,
  const int dj,
  const int dk,
  const int ng,
  const amrex::Real* w,
  amrex::Array4<const amrex::Real> const& q)
{
  amrex::Real sum = 0.0;
  for (int m = -ng; m <= ng; m++) {
    sum += w[m + ng] * q(i + m * di, j + m * dj, k + m * dk, n);
  }
  return sum;
}

class Filter
{

public:
  // Every filter type is a tensor product of the same 1D weights in each
  // direction, so the filter is applied as one 1D pass per direction: the
  // cost per cell is AMREX_SPACEDIM * (2 * ngrow + 1) rather than
  // (2 * ngrow + 1)^AMREX_SPACEDIM.
  explicit Filter(const int type = box, const int fgr = 2)
    : _type(type), _fgr(fgr)
  {

    switch (_type) {

    case box:
      set_box_weights();
      break;

    case gaussian:
      set_gaussian_weights();
      break;

    case box_3pt_approx:
    case gaussian_3pt_approx: // same as box_3pt_approx
      set_box_3pt_approx_weights();
      break;

    case box_5pt_approx:
      set_box_5pt_approx_weights();
      break;

    case box_3pt_optimized_approx:
      set_box_3pt_optimized_approx_weights();
      break;

    case box_5pt_optimized_approx:
      set_box_5pt_optimized_approx_weights();
      break;

    case gaussian_5pt_approx:
      set_gaussian_5pt_approx_weights();
      break;

    case gaussian_3pt_optimized_approx:
      set_gaussian_3pt_optimized_approx_weights();
      break;

    case gaussian_5pt_optimized_approx:
      set_gaussian_5pt_optimized_approx_weights();
      break;

    case no_filter:
    default:
      _fgr = 1;
      _ngrow = 0;
      _nweights = 2 * _ngrow + 1;
      _weights.resize(_nweights);
      _weights[0] = 1.;
      break;

    } // end switch
  }

  // Default destructor
  ~Filter() = default;

  int get_filter_ngrow() const { return _ngrow; }

  const amrex::Vector<amrex::Real>& get_weights() const { return _weights; }

  void apply_filter(const amrex::MultiFab& in, amrex::MultiFab& out);

  void apply_filter(
    const amrex::MultiFab& in,
    amrex::MultiFab& out,
    const int nstart,
    const int ncnt);

  void apply_filter(
    const amrex::Box& cbox, const amrex::FArrayBox& in, amrex::FArrayBox& out);

  void apply_filter(
    const amrex::Box& cbox,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out,
    const int nstart,
    const int ncnt);

  void apply_filter(
    const amrex::Box& box,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out,
    const int nstart,
    const int ncnt,
    const int ncomp);

private:
  int _type;
  int _fgr;
  int _ngrow;
  int _nweights;
  amrex::Vector<amrex::Real> _weights;

  void apply_filter_box(
    const amrex::Box& box,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& qh,
    const int nstart,
    const int ncnt) const;

  void set_box_weights();

  void set_gaussian_weights();

  void set_box_3pt_approx_weights();

  void set_box_5pt_approx_weights();

  void set_box_3pt_optimized_approx_weights();

  void set_box_5pt_optimized_approx_weights();

  void set_gaussian_5pt_approx_weights();

  void set_gaussian_3pt_optimized_approx_weights();

  void set_gaussian_5pt_optimized_approx_weights();
};

#endif /*_FILTER_H*/
