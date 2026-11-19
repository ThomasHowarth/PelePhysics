#ifndef SURFACEDECIMATION_H
#define SURFACEDECIMATION_H

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

// Decimation of iso-surfaces (3D triangles) or iso-lines (2D segments)
// extracted on the dual grid of an AMR hierarchy, by vertex clustering with
// quadrics (P. Lindstrom, "Out-of-core simplification of large polygonal
// models", SIGGRAPH 2000).
//
// Every surface node lies on a dual-grid edge, identified by a key: its level,
// the low cell of the edge and the edge axis. A node belongs to the cluster
// holding that low cell, on blocks of clusterSize^D cells of the level's index
// grid, so the clustering only depends on the mesh and not on the domain
// decomposition. Each cluster accumulates, over the elements touching it, the
// quadric of their planes, an incidence count and the sums of their node
// values. Elements collapsing onto the same clusters are combined into one,
// weighted by their net orientation. All of these are additive: the clustered
// surfaces of several ranks are merged exactly by summing (mergeClusters), and
// the merged clusters are then turned into the decimated surface
// (finalizeClusters). Each cluster becomes one vertex, placed at the minimum of
// its quadric (restricted to the cluster), carrying the cluster average of the
// node values beyond the coordinates.
//
// Open boundaries of the surface (domain boundaries, filter cuts, coarse-fine
// gaps) would otherwise retreat into the surface, since nothing constrains
// the vertices along it. Periodic boundaries count as open: the surface
// extracted by DiagIsosurface is not stitched across them, its facets there
// carry keys a period apart and never pair up, so both sides keep meeting.
// As in Garland & Heckbert (SIGGRAPH 1997), each boundary facet (an edge in
// 3D, an end point in 2D) adds the quadric of the plane through it
// perpendicular to its element, weighted by boundaryWeight times the element
// measure. A facet used by a single element of a rank may still be interior
// to the global surface, its neighbour lying on another rank: such candidates
// are carried along with the clusters, and only those left unpaired once all
// ranks are merged contribute their quadric.

namespace pele::physics::diagnostics {

// Dual-edge key: level, low cell, axis. Cluster keys use the same layout with
// the cluster index in place of the low cell and -1 in place of the axis.
inline constexpr int surfaceKeySize = AMREX_SPACEDIM + 2;
inline constexpr int surfaceNodesPerElt = AMREX_SPACEDIM;
inline constexpr int surfaceFacetKeySize =
  (surfaceNodesPerElt - 1) * surfaceKeySize;
using SurfaceKey = std::array<int, surfaceKeySize>;

struct IntArrayHash
{
  template <std::size_t N>
  std::size_t operator()(const std::array<int, N>& a_key) const noexcept
  {
    std::uint64_t h = 0;
    for (const int v : a_key) {
      amrex::hash_combine(h, v);
    }
    return static_cast<std::size_t>(h);
  }
};

using SurfaceKeyMap = std::unordered_map<SurfaceKey, int, IntArrayHash>;

// Surface whose records (nodes, or clusters) are identified by global keys,
// so that the pieces built on several ranks can be concatenated and merged.
struct KeyedSurface
{
  std::vector<int> keys;         // surfaceKeySize per record
  std::vector<amrex::Real> vals; // fixed number of values per record
  std::vector<int> elts;         // surfaceNodesPerElt record ids per element
  // Clustered surfaces only: net orientation of each element (its record ids
  // being sorted), and boundary facet candidates given by the node keys of
  // the facet (sorted) and the quadric of its boundary plane
  std::vector<int> eltWeights;
  std::vector<int> bkeys;
  std::vector<amrex::Real> bquad;
};

// Number of values per cluster for nodes carrying a_nNodeComp values
// (coordinates first).
int clusterAccSize(int a_nNodeComp);

// Merge the records sharing a key, summing their a_recSize values (a_sum) or
// keeping those of the first one, and renumber the elements accordingly.
// Returns the map from key to merged record.
SurfaceKeyMap mergeByKey(KeyedSurface& a_surf, int a_recSize, bool a_sum);

// Cluster a surface whose records are nodes carrying a_nNodeComp values
// (coordinates first), keyed by their dual edge. a_boundaryWeight scales the
// boundary quadrics (0 disables boundary preservation).
KeyedSurface clusterSurface(
  const KeyedSurface& a_surf,
  int a_nNodeComp,
  int a_clusterSize,
  amrex::Real a_boundaryWeight);

// Merge clustered surfaces concatenated from several ranks: sum the clusters
// sharing a key, combine the elements, and add the quadrics of the boundary
// candidates left unpaired to the clusters of their nodes.
void mergeClusters(KeyedSurface& a_surf, int a_nNodeComp, int a_clusterSize);

// Build the decimated surface from merged clusters: a_nodes receives
// a_nNodeComp values per vertex, a_elts the elements.
void finalizeClusters(
  const KeyedSurface& a_surf,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  int a_clusterSize,
  int a_nNodeComp,
  std::vector<amrex::Real>& a_nodes,
  std::vector<int>& a_elts);

// Measure (area in 3D, length in 2D) of the element a_ids, or of the whole
// surface, whose nodes carry a_nNodeComp values, coordinates first.
amrex::Real elementMeasure(
  const std::vector<amrex::Real>& a_nodes, const int* a_ids, int a_nNodeComp);
amrex::Real surfaceMeasure(
  const std::vector<amrex::Real>& a_nodes,
  const std::vector<int>& a_elts,
  int a_nNodeComp);

} // namespace pele::physics::diagnostics

#endif
