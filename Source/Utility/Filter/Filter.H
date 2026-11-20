#ifndef FILTER_H
#define FILTER_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>

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
  const int ngy = AMREX_D_PICK(0, ng, ng);
  const int ngz = AMREX_D_PICK(0, 0, ng);
  for (int n = -ngz; n <= ngz; n++) {
    const amrex::Real wz = AMREX_D_PICK(1.0, 1.0, w[n + ng]);
    for (int m = -ngy; m <= ngy; m++) {
      const amrex::Real wy = AMREX_D_PICK(1.0, w[m + ng], w[m + ng]);
      for (int l = -ng; l <= ng; l++) {
        qh(i, j, k, nc + nstart) +=
          w[l + ng] * wy * wz * q(i + l, j + m, k + n, nc + nstart);
      }
    }
  }
}

// Every filter type is a tensor product of the same 1D weights in each
// direction, so the filter is applied as one 1D pass per direction: the cost
// per cell is AMREX_SPACEDIM * (2 * ngrow + 1) rather than
// (2 * ngrow + 1)^AMREX_SPACEDIM.
class Filter
{

public:
  explicit Filter(const int type = box, const int fgr = 2) : _fgr(fgr)
  {
    switch (type) {
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
      resize_weights(0);
      _weights[0] = 1.;
      break;
    }
  }

  int get_filter_ngrow() const { return _ngrow; }

  const amrex::Vector<amrex::Real>& get_weights() const { return _weights; }

  void apply_filter(const amrex::MultiFab& in, amrex::MultiFab& out) const;

  // Filters components [nstart, nstart + ncnt) of in into out, which may
  // have ghost cells; in needs get_filter_ngrow() more ghost cells than out.
  void apply_filter(
    const amrex::MultiFab& in,
    amrex::MultiFab& out,
    const int nstart,
    const int ncnt) const;

  void apply_filter(
    const amrex::Box& cbox,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out) const;

  // Filters components [nstart, nstart + ncnt) of in into out on cbox; in
  // needs get_filter_ngrow() ghost cells around cbox in every direction.
  void apply_filter(
    const amrex::Box& cbox,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out,
    const int nstart,
    const int ncnt) const;

  // ncomp is unused, kept for compatibility
  void apply_filter(
    const amrex::Box& box,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out,
    const int nstart,
    const int ncnt,
    const int ncomp) const;

private:
  int _fgr;
  int _ngrow{0};
  amrex::Vector<amrex::Real> _weights;

  void resize_weights(const int ngrow)
  {
    _ngrow = ngrow;
    _weights.resize(2 * ngrow + 1);
  }

  // Symmetric 3 and 5 point weights from the ratios of Sagaut & Grohens
  // (1999), Table I
  void set_3pt_weights_from_ratio(const amrex::Real ratio);

  void set_5pt_weights_from_ratios(
    const amrex::Real ratio1, const amrex::Real ratio2);

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
