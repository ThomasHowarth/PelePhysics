#include "SurfaceDecimation.H"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace pele::physics::diagnostics {

namespace {

constexpr int nodesPerElt = surfaceNodesPerElt;
constexpr int keySize = surfaceKeySize;

// Accumulator layout: quadric (packed symmetric A, then b, then c, with
// error(x) = x.A.x + 2 b.x + c), node count, sums of the node values.
constexpr int nSym = AMREX_SPACEDIM * (AMREX_SPACEDIM + 1) / 2;
constexpr int accB = nSym;
constexpr int accC = nSym + AMREX_SPACEDIM;
constexpr int accCount = accC + 1;
constexpr int accSums = accCount + 1;

// Eigenvalues below this fraction of the largest are treated as zero when
// minimizing a cluster quadric: along those directions, which the planes of
// the cluster do not constrain, the vertex stays at the cluster mean.
constexpr amrex::Real quadricRelTol = 1.0e-3;

using KeyMap = std::unordered_map<SurfaceKey, int, IntArrayHash>;
constexpr int facetKeySize = surfaceFacetKeySize;
constexpr int nodesPerFacet = nodesPerElt - 1;
using FacetKey = std::array<int, facetKeySize>;

constexpr int
symIdx(int a_i, int a_j)
{
  const int r = std::min(a_i, a_j);
  const int c = std::max(a_i, a_j);
  return r * AMREX_SPACEDIM - r * (r - 1) / 2 + (c - r);
}

int
floorDiv(int a_num, int a_den)
{
  return (a_num >= 0) ? a_num / a_den : -((-a_num + a_den - 1) / a_den);
}

// Key of the cluster holding the node with dual-edge key a_nodeKey
SurfaceKey
clusterKeyOf(const int* a_nodeKey, int a_clusterSize)
{
  SurfaceKey ckey;
  ckey[0] = a_nodeKey[0];
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    ckey[1 + d] = floorDiv(a_nodeKey[1 + d], a_clusterSize);
  }
  ckey[keySize - 1] = -1;
  return ckey;
}

// Add a_weight times the quadric of the plane a_normal.x + a_dist = 0
void
addPlaneQuadric(
  amrex::Real* a_quad,
  const amrex::Real a_normal[AMREX_SPACEDIM],
  amrex::Real a_dist,
  amrex::Real a_weight)
{
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    for (int j = i; j < AMREX_SPACEDIM; ++j) {
      a_quad[symIdx(i, j)] += a_weight * a_normal[i] * a_normal[j];
    }
    a_quad[accB + i] += a_weight * a_normal[i] * a_dist;
  }
  a_quad[accC] += a_weight * a_dist * a_dist;
}

// Eigen-decomposition of a symmetric DxD matrix by cyclic Jacobi rotations:
// on exit a_eval holds the eigenvalues and the columns of a_evec the
// corresponding eigenvectors. a_mat is overwritten.
void
symEigen(
  amrex::Real a_mat[AMREX_SPACEDIM][AMREX_SPACEDIM],
  amrex::Real a_eval[AMREX_SPACEDIM],
  amrex::Real a_evec[AMREX_SPACEDIM][AMREX_SPACEDIM])
{
  constexpr int D = AMREX_SPACEDIM;
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < D; ++j) {
      a_evec[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
  for (int sweep = 0; sweep < 50; ++sweep) {
    amrex::Real off = 0.0;
    amrex::Real diag = 0.0;
    for (int p = 0; p < D; ++p) {
      diag += a_mat[p][p] * a_mat[p][p];
      for (int q = p + 1; q < D; ++q) {
        off += a_mat[p][q] * a_mat[p][q];
      }
    }
    if (off <= 1.0e-30 * diag) {
      break;
    }
    for (int p = 0; p < D; ++p) {
      for (int q = p + 1; q < D; ++q) {
        if (a_mat[p][q] == 0.0) {
          continue;
        }
        const amrex::Real theta =
          (a_mat[q][q] - a_mat[p][p]) / (2.0 * a_mat[p][q]);
        const amrex::Real t =
          std::copysign(1.0, theta) /
          (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        const amrex::Real c = 1.0 / std::sqrt(t * t + 1.0);
        const amrex::Real s = t * c;
        for (int k = 0; k < D; ++k) {
          const amrex::Real akp = a_mat[k][p];
          const amrex::Real akq = a_mat[k][q];
          a_mat[k][p] = c * akp - s * akq;
          a_mat[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < D; ++k) {
          const amrex::Real apk = a_mat[p][k];
          const amrex::Real aqk = a_mat[q][k];
          a_mat[p][k] = c * apk - s * aqk;
          a_mat[q][k] = s * apk + c * aqk;
        }
        for (int k = 0; k < D; ++k) {
          const amrex::Real vkp = a_evec[k][p];
          const amrex::Real vkq = a_evec[k][q];
          a_evec[k][p] = c * vkp - s * vkq;
          a_evec[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  for (int i = 0; i < D; ++i) {
    a_eval[i] = a_mat[i][i];
  }
}

// Unit normal of an element and its measure. Returns false for a zero-measure
// element, whose plane is undefined.
bool
elementPlane(
  const std::vector<amrex::Real>& a_nodes,
  const int* a_ids,
  int a_nNodeComp,
  amrex::Real a_normal[AMREX_SPACEDIM],
  amrex::Real& a_measure)
{
  const amrex::Real* p0 =
    &a_nodes[static_cast<std::size_t>(a_ids[0]) * a_nNodeComp];
  const amrex::Real* p1 =
    &a_nodes[static_cast<std::size_t>(a_ids[1]) * a_nNodeComp];
#if AMREX_SPACEDIM == 2
  a_normal[0] = -(p1[1] - p0[1]);
  a_normal[1] = p1[0] - p0[0];
  const amrex::Real len =
    std::sqrt(a_normal[0] * a_normal[0] + a_normal[1] * a_normal[1]);
  a_measure = len;
#else
  const amrex::Real* p2 =
    &a_nodes[static_cast<std::size_t>(a_ids[2]) * a_nNodeComp];
  const amrex::Real ax = p1[0] - p0[0];
  const amrex::Real ay = p1[1] - p0[1];
  const amrex::Real az = p1[2] - p0[2];
  const amrex::Real bx = p2[0] - p0[0];
  const amrex::Real by = p2[1] - p0[1];
  const amrex::Real bz = p2[2] - p0[2];
  a_normal[0] = ay * bz - az * by;
  a_normal[1] = az * bx - ax * bz;
  a_normal[2] = ax * by - ay * bx;
  const amrex::Real len = std::sqrt(
    a_normal[0] * a_normal[0] + a_normal[1] * a_normal[1] +
    a_normal[2] * a_normal[2]);
  a_measure = 0.5 * len;
#endif
  if (len <= 0.0) {
    return false;
  }
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    a_normal[d] /= len;
  }
  return true;
}

} // namespace

int
clusterQuadricSize()
{
  return accCount;
}

int
clusterAccSize(int a_nNodeComp)
{
  return accSums + a_nNodeComp;
}

amrex::Real
elementMeasure(
  const std::vector<amrex::Real>& a_nodes, const int* a_ids, int a_nNodeComp)
{
  amrex::Real normal[AMREX_SPACEDIM];
  amrex::Real measure = 0.0;
  elementPlane(a_nodes, a_ids, a_nNodeComp, normal, measure);
  return measure;
}

ClusteredSurface
clusterSurface(
  const std::vector<amrex::Real>& a_nodes,
  const std::vector<int>& a_keys,
  const std::vector<int>& a_elts,
  int a_nNodeComp,
  int a_clusterSize,
  amrex::Real a_boundaryWeight)
{
  AMREX_ALWAYS_ASSERT(a_clusterSize >= 1);
  const bool doBoundary = a_boundaryWeight > 0.0;
  const int nAcc = clusterAccSize(a_nNodeComp);
  const int nNodes = static_cast<int>(a_keys.size()) / keySize;

  ClusteredSurface out;
  KeyMap clusterMap;
  std::vector<int> nodeCluster(nNodes);
  for (int g = 0; g < nNodes; ++g) {
    const SurfaceKey ckey = clusterKeyOf(
      &a_keys[static_cast<std::size_t>(g) * keySize], a_clusterSize);
    const auto res =
      clusterMap.emplace(ckey, static_cast<int>(clusterMap.size()));
    if (res.second) {
      out.keys.insert(out.keys.end(), ckey.begin(), ckey.end());
      out.acc.resize(out.acc.size() + nAcc, 0.0);
    }
    const int cl = res.first->second;
    nodeCluster[g] = cl;
    amrex::Real* acc = &out.acc[static_cast<std::size_t>(cl) * nAcc];
    const amrex::Real* node =
      &a_nodes[static_cast<std::size_t>(g) * a_nNodeComp];
    acc[accCount] += 1.0;
    for (int n = 0; n < a_nNodeComp; ++n) {
      acc[accSums + n] += node[n];
    }
  }

  // Boundary facet candidates: facet -> (index in out.bquad, use count)
  std::unordered_map<FacetKey, std::pair<int, int>, IntArrayHash> facets;
  constexpr int nQ = accCount;

  const std::size_t nElts = a_elts.size() / nodesPerElt;
  for (std::size_t e = 0; e < nElts; ++e) {
    const int* ids = &a_elts[e * nodesPerElt];
    std::array<int, nodesPerElt> cls{};
    for (int k = 0; k < nodesPerElt; ++k) {
      cls[k] = nodeCluster[ids[k]];
    }

    // Area-weighted (length-weighted in 2D) plane quadric, added to the
    // cluster of each node of the element
    amrex::Real normal[AMREX_SPACEDIM];
    amrex::Real w = 0.0;
    if (elementPlane(a_nodes, ids, a_nNodeComp, normal, w)) {
      const amrex::Real* p0 =
        &a_nodes[static_cast<std::size_t>(ids[0]) * a_nNodeComp];
      amrex::Real dpl = 0.0;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dpl -= normal[d] * p0[d];
      }
      for (const int cl : cls) {
        addPlaneQuadric(
          &out.acc[static_cast<std::size_t>(cl) * nAcc], normal, dpl, w);
      }
    }

    // Facets of the element: edges (3D) or end points (2D). The facet
    // starting at local node f holds nodes f, ..., f + nodesPerFacet - 1.
    for (int f = 0; doBoundary && f < nodesPerElt; ++f) {
      std::array<int, nodesPerFacet> fn{};
      for (int k = 0; k < nodesPerFacet; ++k) {
        fn[k] = ids[(f + k) % nodesPerElt];
      }
      std::array<const int*, nodesPerFacet> fk{};
      for (int k = 0; k < nodesPerFacet; ++k) {
        fk[k] = &a_keys[static_cast<std::size_t>(fn[k]) * keySize];
      }
      std::sort(fk.begin(), fk.end(), [](const int* a, const int* b) {
        return std::lexicographical_compare(a, a + keySize, b, b + keySize);
      });
      FacetKey fkey;
      for (int k = 0; k < nodesPerFacet; ++k) {
        std::copy_n(fk[k], keySize, fkey.begin() + k * keySize);
      }
      const auto res = facets.emplace(
        fkey, std::make_pair(static_cast<int>(facets.size()), 0));
      ++res.first->second.second;
      if (!res.second) {
        continue;
      }

      // Plane through the facet, perpendicular to the element
      amrex::Real bn[AMREX_SPACEDIM] = {0.0};
      const amrex::Real* pa =
        &a_nodes[static_cast<std::size_t>(fn[0]) * a_nNodeComp];
#if AMREX_SPACEDIM == 2
      // In 2D the facet is an end point and the plane is normal to the
      // segment
      if (w > 0.0) {
        bn[0] = normal[1];
        bn[1] = -normal[0];
      }
#else
      const amrex::Real* pb =
        &a_nodes[static_cast<std::size_t>(fn[1]) * a_nNodeComp];
      const amrex::Real tx = pb[0] - pa[0];
      const amrex::Real ty = pb[1] - pa[1];
      const amrex::Real tz = pb[2] - pa[2];
      if (w > 0.0) {
        bn[0] = ty * normal[2] - tz * normal[1];
        bn[1] = tz * normal[0] - tx * normal[2];
        bn[2] = tx * normal[1] - ty * normal[0];
        const amrex::Real bl =
          std::sqrt(bn[0] * bn[0] + bn[1] * bn[1] + bn[2] * bn[2]);
        for (auto& v : bn) {
          v = (bl > 0.0) ? v / bl : 0.0;
        }
      }
#endif
      amrex::Real bd = 0.0;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bd -= bn[d] * pa[d];
      }
      out.bkeys.insert(out.bkeys.end(), fkey.begin(), fkey.end());
      out.bquad.resize(out.bquad.size() + nQ, 0.0);
      addPlaneQuadric(
        &out.bquad[out.bquad.size() - nQ], bn, bd, a_boundaryWeight * w);
    }

    // Cluster keys are global, so an element collapsing here collapses in
    // the merged surface too
    bool degenerate = false;
    for (int a = 0; a < nodesPerElt; ++a) {
      for (int b = a + 1; b < nodesPerElt; ++b) {
        degenerate = degenerate || (cls[a] == cls[b]);
      }
    }
    if (!degenerate) {
      out.elts.insert(out.elts.end(), cls.begin(), cls.end());
    }
  }

  // Keep the facets used once on this rank: the global boundary, plus the
  // facets shared with another rank, paired up in mergeClusters
  std::vector<int> bkeys;
  std::vector<amrex::Real> bquad;
  for (const auto& [fkey, info] : facets) {
    if (info.second == 1) {
      const auto k0 = static_cast<std::ptrdiff_t>(info.first);
      bkeys.insert(
        bkeys.end(), out.bkeys.begin() + k0 * facetKeySize,
        out.bkeys.begin() + (k0 + 1) * facetKeySize);
      bquad.insert(
        bquad.end(), out.bquad.begin() + k0 * nQ,
        out.bquad.begin() + (k0 + 1) * nQ);
    }
  }
  out.bkeys = std::move(bkeys);
  out.bquad = std::move(bquad);
  return out;
}

ClusteredSurface
mergeClusters(
  const ClusteredSurface& a_surf, int a_nNodeComp, int a_clusterSize)
{
  const int nAcc = clusterAccSize(a_nNodeComp);
  const int nClusters = static_cast<int>(a_surf.keys.size()) / keySize;

  ClusteredSurface out;
  KeyMap keyMap;
  keyMap.reserve(nClusters);
  std::vector<int> remap(nClusters);
  for (int cl = 0; cl < nClusters; ++cl) {
    SurfaceKey key;
    std::copy_n(
      &a_surf.keys[static_cast<std::size_t>(cl) * keySize], keySize,
      key.begin());
    const auto res = keyMap.emplace(key, static_cast<int>(keyMap.size()));
    const amrex::Real* src = &a_surf.acc[static_cast<std::size_t>(cl) * nAcc];
    if (res.second) {
      out.keys.insert(out.keys.end(), key.begin(), key.end());
      out.acc.insert(out.acc.end(), src, src + nAcc);
    } else {
      amrex::Real* dst =
        &out.acc[static_cast<std::size_t>(res.first->second) * nAcc];
      for (int n = 0; n < nAcc; ++n) {
        dst[n] += src[n];
      }
    }
    remap[cl] = res.first->second;
  }

  out.elts.resize(a_surf.elts.size());
  for (std::size_t n = 0; n < a_surf.elts.size(); ++n) {
    out.elts[n] = remap[a_surf.elts[n]];
  }

  // Pair up the boundary candidates: a facet found twice is shared by two
  // ranks and interior to the surface. The others are on the boundary and
  // constrain the clusters of their nodes.
  constexpr int nQ = accCount;
  const int nCand = static_cast<int>(a_surf.bkeys.size()) / facetKeySize;
  std::unordered_map<FacetKey, std::pair<int, int>, IntArrayHash> facets;
  facets.reserve(nCand);
  for (int c = 0; c < nCand; ++c) {
    FacetKey fkey;
    std::copy_n(
      &a_surf.bkeys[static_cast<std::size_t>(c) * facetKeySize], facetKeySize,
      fkey.begin());
    ++facets.emplace(fkey, std::make_pair(c, 0)).first->second.second;
  }
  for (const auto& [fkey, info] : facets) {
    if (info.second != 1) {
      continue;
    }
    const amrex::Real* q =
      &a_surf.bquad[static_cast<std::size_t>(info.first) * nQ];
    for (int k = 0; k < nodesPerFacet; ++k) {
      const auto found =
        keyMap.find(clusterKeyOf(&fkey[k * keySize], a_clusterSize));
      AMREX_ALWAYS_ASSERT(found != keyMap.end());
      amrex::Real* dst =
        &out.acc[static_cast<std::size_t>(found->second) * nAcc];
      for (int n = 0; n < nQ; ++n) {
        dst[n] += q[n];
      }
    }
  }
  return out;
}

amrex::Real
finalizeClusters(
  const ClusteredSurface& a_surf,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  int a_clusterSize,
  int a_nNodeComp,
  std::vector<amrex::Real>& a_nodes,
  std::vector<int>& a_elts)
{
  const int nAcc = clusterAccSize(a_nNodeComp);
  const int nClusters = static_cast<int>(a_surf.keys.size()) / keySize;

  // Elements mapping onto the same clusters are merged, summing their
  // orientations: a pair with opposite orientations is a fold pinched by the
  // clustering and cancels out, and the net orientation of the others is
  // kept. Elements are emitted in order of first occurrence and clusters
  // numbered in order of first use, so the output only depends on the input
  // order.
  using EltKey = std::array<int, nodesPerElt>;
  std::unordered_map<EltKey, int, IntArrayHash> eltIndex;
  std::vector<EltKey> uniqueElts;
  std::vector<int> netOrient;
  const std::size_t nElts = a_surf.elts.size() / nodesPerElt;
  for (std::size_t e = 0; e < nElts; ++e) {
    EltKey ids;
    std::copy_n(&a_surf.elts[e * nodesPerElt], nodesPerElt, ids.begin());
    // Sort, tracking the parity of the permutation
    int parity = 1;
    for (int a = 0; a < nodesPerElt; ++a) {
      for (int b = 0; b + 1 < nodesPerElt - a; ++b) {
        if (ids[b] > ids[b + 1]) {
          std::swap(ids[b], ids[b + 1]);
          parity = -parity;
        }
      }
    }
    bool degenerate = false;
    for (int k = 1; k < nodesPerElt; ++k) {
      degenerate = degenerate || (ids[k] == ids[k - 1]);
    }
    if (degenerate) {
      continue;
    }
    const auto res = eltIndex.emplace(ids, static_cast<int>(uniqueElts.size()));
    if (res.second) {
      uniqueElts.push_back(ids);
      netOrient.push_back(0);
    }
    netOrient[res.first->second] += parity;
  }

  std::vector<int> newID(nClusters, -1);
  int nUsed = 0;
  a_elts.clear();
  for (std::size_t u = 0; u < uniqueElts.size(); ++u) {
    if (netOrient[u] == 0) {
      continue;
    }
    EltKey ids = uniqueElts[u];
    if (netOrient[u] < 0) {
      std::swap(ids[0], ids[1]);
    }
    for (const int id : ids) {
      if (newID[id] < 0) {
        newID[id] = nUsed++;
      }
      a_elts.push_back(newID[id]);
    }
  }

  a_nodes.assign(static_cast<std::size_t>(nUsed) * a_nNodeComp, 0.0);
  for (int cl = 0; cl < nClusters; ++cl) {
    if (newID[cl] < 0) {
      continue;
    }
    const amrex::Real* acc = &a_surf.acc[static_cast<std::size_t>(cl) * nAcc];
    amrex::Real* node =
      &a_nodes[static_cast<std::size_t>(newID[cl]) * a_nNodeComp];
    const amrex::Real count = acc[accCount];
    for (int n = 0; n < a_nNodeComp; ++n) {
      node[n] = acc[accSums + n] / count;
    }

    // Minimize the quadric about the cluster mean x0:
    // x = x0 + A^+ (-b - A x0), A^+ dropping the near-null eigen-directions
    amrex::Real mat[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real rhs[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      rhs[i] = -acc[accB + i];
      for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        mat[i][j] = acc[symIdx(i, j)];
        rhs[i] -= acc[symIdx(i, j)] * node[j];
      }
    }
    amrex::Real eval[AMREX_SPACEDIM];
    amrex::Real evec[AMREX_SPACEDIM][AMREX_SPACEDIM];
    symEigen(mat, eval, evec);
    amrex::Real lmax = 0.0;
    for (const amrex::Real l : eval) {
      lmax = std::max(lmax, l);
    }
    for (int k = 0; k < AMREX_SPACEDIM; ++k) {
      if (lmax <= 0.0 || eval[k] <= quadricRelTol * lmax) {
        continue;
      }
      amrex::Real proj = 0.0;
      for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        proj += evec[i][k] * rhs[i];
      }
      for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        node[i] += evec[i][k] * proj / eval[k];
      }
    }

    // Keep the vertex within the region spanned by the cluster's dual edges:
    // from the center of its first cell to the center of the first cell of
    // the next cluster
    const int* key = &a_surf.keys[static_cast<std::size_t>(cl) * keySize];
    const amrex::Geometry& geom = a_geoms[key[0]];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      const amrex::Real dx = geom.CellSize(d);
      const amrex::Real lo =
        geom.ProbLo(d) + (key[1 + d] * a_clusterSize + 0.5) * dx;
      node[d] = std::clamp(node[d], lo, lo + a_clusterSize * dx);
    }
  }

  amrex::Real measure = 0.0;
  for (std::size_t e = 0; e < a_elts.size(); e += nodesPerElt) {
    measure += elementMeasure(a_nodes, &a_elts[e], a_nNodeComp);
  }
  return measure;
}

} // namespace pele::physics::diagnostics
