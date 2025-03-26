#ifndef _stateinflow_H_
#define _stateinflow_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <PeleLM_Index.H>
#include "mechanism.h"

namespace pele {
namespace physics {
namespace stateinflow {

struct StateParm
{
  // state data file
  std::string m_state_file = "";

  // Geometry of the state data
  amrex::GpuArray<int, AMREX_SPACEDIM> npboxcells = {{0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pboxlo = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pboxsize = {{0.0}};
  int dir;                       // Direction where to use this TurbParm
  amrex::Orientation::Side side; // Side on which to use this TurbParm
  amrex::Real time_shift =
    0.0; // Enable to offset the location in the turb file

  int nplane = 32; // Number of turb planes stored at once in memory
  amrex::FArrayBox* sdata = nullptr; // Currently loaded chunk of turb data
  amrex::Real szlo = -1.0e12;        // Position of the current chunk low plane
  amrex::Real szhi = -1.0e11;        // Position of the current chunk high plane

  bool isswirltype = false; // Unused: enable rotating velocity
  amrex::Real state_scale_loc =
    1.; // Spatial scale relating turbfile geom -> actual case scale
  //amrex::Real turb_scale_vel = 1.; // Velocity scaling
  amrex::Real state_conv_vel =
    1.; // Mean inlet velocity, used to move through data planes

  amrex::Vector<long> m_offset;
  int kmax;     // Number of plane in Turbfile
  long* offset; // Binary offset of each data plane in Turbfile
  long offset_size;
};

struct StateInflow
{
public:
  StateInflow() = default;

  ~StateInflow() = default;

  void init(amrex::Geometry const& geom);

  void add_state(
    amrex::Box const& bx,
    amrex::FArrayBox& data,
    const int dcomp,
    amrex::Geometry const& geom,
    const amrex::Real time,
    const int dir,
    const amrex::Orientation::Side& side, int ncomp);

  bool is_initialized() const { return stateinflow_initialized; }

  static void read_state_planes(StateParm& a_sp, amrex::Real z);

  static void read_one_state_plane(StateParm& a_sp, int iplane, int k);

  static void fill_state_plane(
    StateParm& a_sp,
    const amrex::Vector<amrex::Real>& x,
    const amrex::Vector<amrex::Real>& y,
    amrex::Real z,
    amrex::FArrayBox& s);

  static void set_state(
    int normDir,
    int transDir1,
    int transDir2,
    amrex::FArrayBox& s,
    amrex::FArrayBox& data,
    const int dcomp, int ncomp);

private:
  amrex::Vector<StateParm> sp;
  bool stateinflow_initialized = false;
};
} // namespace stateinflow

} // namespace physics
} // namespace pele

#endif
