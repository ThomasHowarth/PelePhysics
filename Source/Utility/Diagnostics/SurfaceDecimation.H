#ifndef SURFACEDECIMATION_H
#define SURFACEDECIMATION_H

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <cstddef>
#include <functional>
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
// decomposition. Each cluster accumulates the quadric of the element planes
// touching it, its node count and the sums of its node values. All of these
// are additive: clusters split between ranks are merged exactly by summing
// their accumulators (mergeClusters), and the merged clusters are then turned
// into the decimated surface (finalizeClusters). Each cluster becomes one
// vertex, placed at the minimum of its quadric (restricted to the cluster),
// carrying the cluster average of the node values beyond the coordinates.
//
// Open boundaries of the surface (domain boundaries, filter cuts, coarse-fine
// gaps) would otherwise retreat into the surface, since nothing constrains
// the vertices along it. Periodic boundaries count as open: the surface
// extracted by DiagIsosurface is not stitched across them, its facets there
// carry keys a period apart and never pair up, so both sides keep meeting. As
// in Garland & Heckbert (SIGGRAPH 1997), each boundary facet (an edge in 3D, an
// end point in 2D) adds the quadric of the plane through it perpendicular to
// its element, weighted by boundaryWeight times the element measure. A facet
// used by a single element of a rank may still be interior to the global
// surface, its neighbour lying on another rank: such candidates are carried
// along with the clusters, and only those left unpaired once all ranks are
// merged contribute their quadric.

namespace pele::physics::diagnostics {

// Dual-edge key: level, low cell, axis. Cluster keys use the same layout with
// the cluster index in place of the low cell and -1 in place of the axis.
inline constexpr int surfaceKeySize = AMREX_SPACEDIM + 2;
inline constexpr int surfaceNodesPerElt = AMREX_SPACEDIM;
using SurfaceKey = std::array<int, surfaceKeySize>;

struct IntArrayHash
{
  template <std::size_t N>
  std::size_t operator()(const std::array<int, N>& a_key) const noexcept
  {
    std::size_t h = 0;
    for (const int v : a_key) {
      h ^= std::hash<int>{}(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

// Clustered surface: one key and one accumulator per cluster, elements given
// as cluster indices.
struct ClusteredSurface
{
  std::vector<int> keys;        // surfaceKeySize per cluster
  std::vector<amrex::Real> acc; // clusterAccSize(nNodeComp) per cluster
  std::vector<int> elts;        // surfaceNodesPerElt per element
  // Boundary facet candidates: node keys of the facet (sorted), and the
  // quadric of its boundary plane (clusterQuadricSize() values)
  std::vector<int> bkeys;
  std::vector<amrex::Real> bquad;
};

// Number of ints identifying a boundary facet candidate in bkeys
inline constexpr int surfaceFacetKeySize =
  (surfaceNodesPerElt - 1) * surfaceKeySize;

// Number of values of a quadric
int clusterQuadricSize();

// Number of accumulated values per cluster for nodes carrying a_nNodeComp
// values (coordinates first).
int clusterAccSize(int a_nNodeComp);

// Cluster a surface. a_nodes holds a_nNodeComp values per node (coordinates
// first), a_keys the dual-edge key of each node, a_elts the elements as node
// indices. Elements collapsing onto fewer than surfaceNodesPerElt clusters
// are dropped. a_boundaryWeight scales the boundary quadrics (0 disables
// boundary preservation).
ClusteredSurface clusterSurface(
  const std::vector<amrex::Real>& a_nodes,
  const std::vector<int>& a_keys,
  const std::vector<int>& a_elts,
  int a_nNodeComp,
  int a_clusterSize,
  amrex::Real a_boundaryWeight);

// Merge clusters sharing the same key, e.g. after gathering the clustered
// surfaces of several ranks, and renumber the elements accordingly. Boundary
// candidates found twice are interior and dropped; the quadric of the others
// is added to the clusters of their nodes. The result has no candidates left.
ClusteredSurface mergeClusters(
  const ClusteredSurface& a_surf, int a_nNodeComp, int a_clusterSize);

// Build the decimated surface from merged clusters: a_nodes receives
// a_nNodeComp values per vertex, a_elts the unique, non-degenerate elements.
// Returns the measure (area in 3D, length in 2D) of the decimated surface.
amrex::Real finalizeClusters(
  const ClusteredSurface& a_surf,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  int a_clusterSize,
  int a_nNodeComp,
  std::vector<amrex::Real>& a_nodes,
  std::vector<int>& a_elts);

// Measure (area in 3D, length in 2D) of the element a_ids of a surface whose
// nodes carry a_nNodeComp values, coordinates first.
amrex::Real elementMeasure(
  const std::vector<amrex::Real>& a_nodes, const int* a_ids, int a_nNodeComp);

} // namespace pele::physics::diagnostics

#endif
