#include "SurfaceDecimation.H"

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>

#include <algorithm>
#include <cmath>

namespace pele::physics::diagnostics {

namespace {

constexpr int nodesPerElt = surfaceNodesPerElt;
constexpr int nodesPerFacet = nodesPerElt - 1;
constexpr int keySize = surfaceKeySize;
constexpr int facetKeySize = surfaceFacetKeySize;
using EltKey = std::array<int, nodesPerElt>;
using FacetKey = std::array<int, facetKeySize>;

// Cluster record layout: quadric (packed symmetric A, then b, then c, with
// error(x) = x.A.x + 2 b.x + c), incidence count, sums of the node values.
constexpr int nSym = AMREX_SPACEDIM * (AMREX_SPACEDIM + 1) / 2;
constexpr int accB = nSym;
constexpr int accC = nSym + AMREX_SPACEDIM;
constexpr int nQuad = accC + 1;
constexpr int accCount = nQuad;
constexpr int accSums = accCount + 1;

// Eigenvalues below this fraction of the largest are treated as zero when
// minimizing a cluster quadric: along those directions, which the planes of
// the cluster do not constrain, the vertex stays at the cluster mean.
constexpr amrex::Real quadricRelTol = 1.0e-3;

constexpr int
symIdx(int a_i, int a_j)
{
  const int r = std::min(a_i, a_j);
  const int c = std::max(a_i, a_j);
  return r * AMREX_SPACEDIM - r * (r - 1) / 2 + (c - r);
}

// Key of the cluster holding the node with dual-edge key a_nodeKey
SurfaceKey
clusterKeyOf(const int* a_nodeKey, int a_clusterSize)
{
  SurfaceKey ckey;
  ckey[0] = a_nodeKey[0];
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    ckey[1 + d] = amrex::coarsen(a_nodeKey[1 + d], a_clusterSize);
  }
  ckey[keySize - 1] = -1;
  return ckey;
}

// Add a_weight times the quadric of the plane a_normal.x + a_dist = 0
void
addPlaneQuadric(
  amrex::Real* a_quad,
  const amrex::RealVect& a_normal,
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

amrex::RealVect
nodePos(const std::vector<amrex::Real>& a_nodes, int a_id, int a_nNodeComp)
{
  return amrex::RealVect(
    &a_nodes[static_cast<std::size_t>(a_id) * a_nNodeComp]);
}

// Unit normal of an element and its measure. The normal is zero for a
// zero-measure element, whose plane is undefined.
amrex::Real
elementPlane(
  const std::vector<amrex::Real>& a_nodes,
  const int* a_ids,
  int a_nNodeComp,
  amrex::RealVect& a_normal)
{
  const amrex::RealVect p0 = nodePos(a_nodes, a_ids[0], a_nNodeComp);
  const amrex::RealVect t = nodePos(a_nodes, a_ids[1], a_nNodeComp) - p0;
#if AMREX_SPACEDIM == 2
  a_normal = amrex::RealVect(-t[1], t[0]);
  const amrex::Real len = a_normal.vectorLength();
  const amrex::Real measure = len;
#else
  a_normal = t.crossProduct(nodePos(a_nodes, a_ids[2], a_nNodeComp) - p0);
  const amrex::Real len = a_normal.vectorLength();
  const amrex::Real measure = 0.5 * len;
#endif
  if (len > 0.0) {
    a_normal /= len;
  }
  return measure;
}

// Combine the elements spanning the same records: each element's record ids
// are sorted, its weight (1 if a_weights is empty) multiplied by the parity of
// that permutation, and the weights of identical elements summed. Opposite
// orientations thus cancel out, e.g. the two sides of a fold pinched by the
// clustering. Elements with repeated records or a zero net weight are
// dropped; the others are kept in order of first occurrence.
void
combineElements(std::vector<int>& a_elts, std::vector<int>& a_weights)
{
  struct Item
  {
    EltKey ids;
    int weight;
    int first;
  };
  const int nElts = static_cast<int>(a_elts.size()) / nodesPerElt;
  std::vector<Item> items;
  items.reserve(nElts);
  for (int e = 0; e < nElts; ++e) {
    Item it{{}, a_weights.empty() ? 1 : a_weights[e], e};
    std::copy_n(
      &a_elts[static_cast<std::size_t>(e) * nodesPerElt], nodesPerElt,
      it.ids.begin());
    for (int a = 0; a < nodesPerElt; ++a) {
      for (int b = 0; b + 1 < nodesPerElt - a; ++b) {
        if (it.ids[b] > it.ids[b + 1]) {
          std::swap(it.ids[b], it.ids[b + 1]);
          it.weight = -it.weight;
        }
      }
    }
    if (std::adjacent_find(it.ids.begin(), it.ids.end()) == it.ids.end()) {
      items.push_back(it);
    }
  }

  std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
    return (a.ids != b.ids) ? a.ids < b.ids : a.first < b.first;
  });
  std::size_t nRuns = 0;
  for (std::size_t i = 0; i < items.size(); ++i) {
    if (nRuns > 0 && items[nRuns - 1].ids == items[i].ids) {
      items[nRuns - 1].weight += items[i].weight;
    } else {
      items[nRuns++] = items[i];
    }
  }
  items.resize(nRuns);
  std::erase_if(items, [](const Item& a) { return a.weight == 0; });
  std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
    return a.first < b.first;
  });

  a_elts.resize(items.size() * nodesPerElt);
  a_weights.resize(items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    std::copy_n(items[i].ids.begin(), nodesPerElt, &a_elts[i * nodesPerElt]);
    a_weights[i] = items[i].weight;
  }
}

// Boundary facet candidates of a node-keyed surface: the facets (edges in 3D,
// end points in 2D) used by a single element, with the quadric of the plane
// through them perpendicular to their element. Facets are first counted on
// the rank-local node ids, which are unique within the rank, and only the
// candidates get their global key and quadric.
void
addBoundaryCandidates(
  const KeyedSurface& a_surf,
  int a_nNodeComp,
  amrex::Real a_boundaryWeight,
  KeyedSurface& a_out)
{
  const int nElts = static_cast<int>(a_surf.elts.size()) / nodesPerElt;
  // Facet starting at local node f of element e: nodes f, ...,
  // f + nodesPerFacet - 1, packed into one integer
  auto facetId = [&](int a_e, int a_f) {
    const int* ids = &a_surf.elts[static_cast<std::size_t>(a_e) * nodesPerElt];
#if AMREX_SPACEDIM == 2
    return static_cast<std::uint64_t>(ids[a_f]);
#else
    const auto a = static_cast<std::uint64_t>(ids[a_f]);
    const auto b = static_cast<std::uint64_t>(ids[(a_f + 1) % nodesPerElt]);
    return (std::min(a, b) << 32) | std::max(a, b);
#endif
  };

  std::unordered_map<std::uint64_t, int> uses;
  uses.reserve(static_cast<std::size_t>(nElts) * nodesPerElt);
  for (int e = 0; e < nElts; ++e) {
    for (int f = 0; f < nodesPerElt; ++f) {
      ++uses[facetId(e, f)];
    }
  }

  for (int e = 0; e < nElts; ++e) {
    const int* ids = &a_surf.elts[static_cast<std::size_t>(e) * nodesPerElt];
    for (int f = 0; f < nodesPerElt; ++f) {
      if (uses[facetId(e, f)] != 1) {
        continue;
      }
      std::array<SurfaceKey, nodesPerFacet> nodeKeys;
      for (int k = 0; k < nodesPerFacet; ++k) {
        std::copy_n(
          &a_surf.keys
             [static_cast<std::size_t>(ids[(f + k) % nodesPerElt]) * keySize],
          keySize, nodeKeys[k].begin());
      }
      std::sort(nodeKeys.begin(), nodeKeys.end());
      for (const auto& key : nodeKeys) {
        a_out.bkeys.insert(a_out.bkeys.end(), key.begin(), key.end());
      }

      // Plane through the facet, perpendicular to the element: normal to
      // the segment through its end point in 2D, containing the edge and the
      // element normal in 3D
      amrex::RealVect normal;
      const amrex::Real w = elementPlane(a_surf.vals, ids, a_nNodeComp, normal);
      const amrex::RealVect pa = nodePos(a_surf.vals, ids[f], a_nNodeComp);
#if AMREX_SPACEDIM == 2
      amrex::RealVect bn(normal[1], -normal[0]);
#else
      amrex::RealVect bn =
        (nodePos(a_surf.vals, ids[(f + 1) % nodesPerElt], a_nNodeComp) - pa)
          .crossProduct(normal);
      if (const amrex::Real bl = bn.vectorLength(); bl > 0.0) {
        bn /= bl;
      }
#endif
      a_out.bquad.resize(a_out.bquad.size() + nQuad, 0.0);
      addPlaneQuadric(
        &a_out.bquad[a_out.bquad.size() - nQuad], bn, -bn.dotProduct(pa),
        a_boundaryWeight * w);
    }
  }
}

} // namespace

int
clusterAccSize(int a_nNodeComp)
{
  return accSums + a_nNodeComp;
}

amrex::Real
elementMeasure(
  const std::vector<amrex::Real>& a_nodes, const int* a_ids, int a_nNodeComp)
{
  amrex::RealVect normal;
  return elementPlane(a_nodes, a_ids, a_nNodeComp, normal);
}

amrex::Real
surfaceMeasure(
  const std::vector<amrex::Real>& a_nodes,
  const std::vector<int>& a_elts,
  int a_nNodeComp)
{
  amrex::Real measure = 0.0;
  for (std::size_t e = 0; e < a_elts.size(); e += nodesPerElt) {
    measure += elementMeasure(a_nodes, &a_elts[e], a_nNodeComp);
  }
  return measure;
}

SurfaceKeyMap
mergeByKey(KeyedSurface& a_surf, int a_recSize, bool a_sum)
{
  const int nRec = static_cast<int>(a_surf.keys.size()) / keySize;
  SurfaceKeyMap keyMap;
  keyMap.reserve(nRec);
  std::vector<int> remap(nRec);
  // Compact in place: the record merged into slot n is never before slot n
  int nUnique = 0;
  for (int r = 0; r < nRec; ++r) {
    SurfaceKey key;
    std::copy_n(
      &a_surf.keys[static_cast<std::size_t>(r) * keySize], keySize,
      key.begin());
    const auto [it, inserted] = keyMap.try_emplace(key, nUnique);
    amrex::Real* src = &a_surf.vals[static_cast<std::size_t>(r) * a_recSize];
    amrex::Real* dst =
      &a_surf.vals[static_cast<std::size_t>(it->second) * a_recSize];
    if (inserted) {
      if (nUnique != r) {
        std::copy_n(
          key.begin(), keySize,
          &a_surf.keys[static_cast<std::size_t>(nUnique) * keySize]);
        std::copy_n(src, a_recSize, dst);
      }
      ++nUnique;
    } else if (a_sum) {
      for (int n = 0; n < a_recSize; ++n) {
        dst[n] += src[n];
      }
    }
    remap[r] = it->second;
  }
  a_surf.keys.resize(static_cast<std::size_t>(nUnique) * keySize);
  a_surf.vals.resize(static_cast<std::size_t>(nUnique) * a_recSize);
  for (auto& id : a_surf.elts) {
    id = remap[id];
  }
  return keyMap;
}

KeyedSurface
clusterSurface(
  const KeyedSurface& a_surf,
  int a_nNodeComp,
  int a_clusterSize,
  amrex::Real a_boundaryWeight)
{
  AMREX_ALWAYS_ASSERT(a_clusterSize >= 1);
  const int nAcc = clusterAccSize(a_nNodeComp);
  const int nNodes = static_cast<int>(a_surf.keys.size()) / keySize;
  const int nElts = static_cast<int>(a_surf.elts.size()) / nodesPerElt;

  KeyedSurface out;
  SurfaceKeyMap clusterMap;
  clusterMap.reserve(nNodes / a_clusterSize + 1);
  std::vector<int> nodeCluster(nNodes, -1);
  auto clusterOf = [&](int a_node) {
    if (nodeCluster[a_node] < 0) {
      const SurfaceKey ckey = clusterKeyOf(
        &a_surf.keys[static_cast<std::size_t>(a_node) * keySize],
        a_clusterSize);
      const auto [it, inserted] =
        clusterMap.try_emplace(ckey, static_cast<int>(clusterMap.size()));
      if (inserted) {
        out.keys.insert(out.keys.end(), ckey.begin(), ckey.end());
        out.vals.resize(out.vals.size() + nAcc, 0.0);
      }
      nodeCluster[a_node] = it->second;
    }
    return nodeCluster[a_node];
  };

  // Accumulate over the element incidences rather than over the nodes: nodes
  // on box boundaries exist on every rank touching them, while every element
  // exists once, so the result does not depend on the decomposition.
  out.elts.reserve(a_surf.elts.size());
  for (int e = 0; e < nElts; ++e) {
    const int* ids = &a_surf.elts[static_cast<std::size_t>(e) * nodesPerElt];
    amrex::RealVect normal;
    const amrex::Real w = elementPlane(a_surf.vals, ids, a_nNodeComp, normal);
    const amrex::Real dist =
      -normal.dotProduct(nodePos(a_surf.vals, ids[0], a_nNodeComp));
    for (int k = 0; k < nodesPerElt; ++k) {
      const int cl = clusterOf(ids[k]);
      amrex::Real* acc = &out.vals[static_cast<std::size_t>(cl) * nAcc];
      addPlaneQuadric(acc, normal, dist, w);
      acc[accCount] += 1.0;
      const amrex::Real* node =
        &a_surf.vals[static_cast<std::size_t>(ids[k]) * a_nNodeComp];
      for (int n = 0; n < a_nNodeComp; ++n) {
        acc[accSums + n] += node[n];
      }
      out.elts.push_back(cl);
    }
  }
  combineElements(out.elts, out.eltWeights);

  if (a_boundaryWeight > 0.0) {
    addBoundaryCandidates(a_surf, a_nNodeComp, a_boundaryWeight, out);
  }
  return out;
}

void
mergeClusters(KeyedSurface& a_surf, int a_nNodeComp, int a_clusterSize)
{
  const int nAcc = clusterAccSize(a_nNodeComp);
  const SurfaceKeyMap keyMap = mergeByKey(a_surf, nAcc, true);
  combineElements(a_surf.elts, a_surf.eltWeights);

  // A candidate found twice is a facet shared by two ranks, interior to the
  // surface. The others are on the boundary and constrain the clusters of
  // their nodes.
  const int nCand = static_cast<int>(a_surf.bkeys.size()) / facetKeySize;
  auto facetKey = [&](int a_c) {
    FacetKey fkey;
    std::copy_n(
      &a_surf.bkeys[static_cast<std::size_t>(a_c) * facetKeySize], facetKeySize,
      fkey.begin());
    return fkey;
  };
  std::unordered_map<FacetKey, int, IntArrayHash> uses;
  uses.reserve(nCand);
  for (int c = 0; c < nCand; ++c) {
    ++uses[facetKey(c)];
  }
  for (int c = 0; c < nCand; ++c) {
    const FacetKey fkey = facetKey(c);
    if (uses[fkey] != 1) {
      continue;
    }
    const amrex::Real* q = &a_surf.bquad[static_cast<std::size_t>(c) * nQuad];
    for (int k = 0; k < nodesPerFacet; ++k) {
      const auto found =
        keyMap.find(clusterKeyOf(&fkey[k * keySize], a_clusterSize));
      AMREX_ALWAYS_ASSERT(found != keyMap.end());
      amrex::Real* acc =
        &a_surf.vals[static_cast<std::size_t>(found->second) * nAcc];
      for (int n = 0; n < nQuad; ++n) {
        acc[n] += q[n];
      }
    }
  }
  a_surf.bkeys = {};
  a_surf.bquad = {};
}

void
finalizeClusters(
  const KeyedSurface& a_surf,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  int a_clusterSize,
  int a_nNodeComp,
  std::vector<amrex::Real>& a_nodes,
  std::vector<int>& a_elts)
{
  const int nAcc = clusterAccSize(a_nNodeComp);
  const int nClusters = static_cast<int>(a_surf.keys.size()) / keySize;
  const std::size_t nElts = a_surf.eltWeights.size();
  AMREX_ALWAYS_ASSERT(a_surf.elts.size() == nElts * nodesPerElt);

  // Elements with their net orientation, clusters numbered in order of
  // first use
  std::vector<int> newID(nClusters, -1);
  int nUsed = 0;
  a_elts.resize(a_surf.elts.size());
  for (std::size_t e = 0; e < nElts; ++e) {
    EltKey ids;
    std::copy_n(&a_surf.elts[e * nodesPerElt], nodesPerElt, ids.begin());
    if (a_surf.eltWeights[e] < 0) {
      std::swap(ids[0], ids[1]);
    }
    for (int k = 0; k < nodesPerElt; ++k) {
      if (newID[ids[k]] < 0) {
        newID[ids[k]] = nUsed++;
      }
      a_elts[e * nodesPerElt + k] = newID[ids[k]];
    }
  }

  a_nodes.assign(static_cast<std::size_t>(nUsed) * a_nNodeComp, 0.0);
  for (int cl = 0; cl < nClusters; ++cl) {
    if (newID[cl] < 0) {
      continue;
    }
    const amrex::Real* acc = &a_surf.vals[static_cast<std::size_t>(cl) * nAcc];
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
    const amrex::Real lmax = *std::max_element(eval, eval + AMREX_SPACEDIM);
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
}

} // namespace pele::physics::diagnostics
