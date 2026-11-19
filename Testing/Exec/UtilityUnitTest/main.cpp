#include <AMReX.H>
#include <AMReX_Print.H>
#include "Utilities.H"
#include "Table.H"
#include "BlackBoxFunction.H"
#include "SurfaceDecimation.H"

#include <cmath>
#include <vector>

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

#if AMREX_SPACEDIM == 2
namespace {

namespace sd = pele::physics::diagnostics;

// Polylines carry, per node, coordinates and a constant field, keyed by the
// dual cell holding them (as DiagIsosurface keys its iso-line nodes)
constexpr int nNC = 3;
constexpr amrex::Real fieldVal = 7.0;

void
addNode(
  sd::KeyedSurface& a_p, amrex::Real a_x, amrex::Real a_y, amrex::Real a_dx)
{
  a_p.vals.insert(a_p.vals.end(), {a_x, a_y, fieldVal});
  a_p.keys.insert(
    a_p.keys.end(), {0, static_cast<int>(std::floor(a_x / a_dx - 0.5)),
                     static_cast<int>(std::floor(a_y / a_dx - 0.5)), 0});
}

// Nodes are spaced by more than sqrt(2) dx, so that no two share a dual cell
sd::KeyedSurface
makeCircle(amrex::Real a_r, amrex::Real a_dx)
{
  sd::KeyedSurface p;
  const int n = static_cast<int>(2.0 * M_PI * a_r / (1.6 * a_dx));
  for (int i = 0; i < n; ++i) {
    const amrex::Real t = 2.0 * M_PI * i / n;
    addNode(p, 0.5 + a_r * std::cos(t), 0.5 + a_r * std::sin(t), a_dx);
    p.elts.insert(p.elts.end(), {i, (i + 1) % n});
  }
  return p;
}

sd::KeyedSurface
makeLine(
  amrex::Real a_x0,
  amrex::Real a_y0,
  amrex::Real a_x1,
  amrex::Real a_y1,
  amrex::Real a_dx)
{
  sd::KeyedSurface p;
  const amrex::Real len = std::hypot(a_x1 - a_x0, a_y1 - a_y0);
  const int n = static_cast<int>(len / (1.6 * a_dx));
  for (int i = 0; i <= n; ++i) {
    const amrex::Real t = static_cast<amrex::Real>(i) / n;
    addNode(p, a_x0 + t * (a_x1 - a_x0), a_y0 + t * (a_y1 - a_y0), a_dx);
    if (i < n) {
      p.elts.insert(p.elts.end(), {i, i + 1});
    }
  }
  return p;
}

// Elements [a_eBegin, a_eEnd) of a_p with the nodes they use, renumbered
sd::KeyedSurface
subPolyline(const sd::KeyedSurface& a_p, int a_eBegin, int a_eEnd)
{
  sd::KeyedSurface sub;
  std::vector<int> local(a_p.keys.size() / sd::surfaceKeySize, -1);
  for (int e = a_eBegin; e < a_eEnd; ++e) {
    for (int k = 0; k < 2; ++k) {
      const int g = a_p.elts[2 * e + k];
      if (local[g] < 0) {
        local[g] = static_cast<int>(sub.keys.size()) / sd::surfaceKeySize;
        sub.vals.insert(
          sub.vals.end(), a_p.vals.begin() + g * nNC,
          a_p.vals.begin() + (g + 1) * nNC);
        sub.keys.insert(
          sub.keys.end(), a_p.keys.begin() + g * sd::surfaceKeySize,
          a_p.keys.begin() + (g + 1) * sd::surfaceKeySize);
      }
      sub.elts.push_back(local[g]);
    }
  }
  return sub;
}

// Append a_b to a_a, as the IO rank does with the clusters of several ranks
void
appendClusters(sd::KeyedSurface& a_a, const sd::KeyedSurface& a_b)
{
  const int offset = static_cast<int>(a_a.keys.size()) / sd::surfaceKeySize;
  a_a.keys.insert(a_a.keys.end(), a_b.keys.begin(), a_b.keys.end());
  a_a.vals.insert(a_a.vals.end(), a_b.vals.begin(), a_b.vals.end());
  for (const int id : a_b.elts) {
    a_a.elts.push_back(id + offset);
  }
  a_a.eltWeights.insert(
    a_a.eltWeights.end(), a_b.eltWeights.begin(), a_b.eltWeights.end());
  a_a.bkeys.insert(a_a.bkeys.end(), a_b.bkeys.begin(), a_b.bkeys.end());
  a_a.bquad.insert(a_a.bquad.end(), a_b.bquad.begin(), a_b.bquad.end());
}

struct Decimated
{
  std::vector<amrex::Real> nodes;
  std::vector<int> elts;
  amrex::Real length{0.0};
};

void
testSurfaceDecimation()
{
  constexpr int nCell = 64;
  constexpr amrex::Real dx = 1.0 / nCell;
  constexpr int clusterSize = 2;
  constexpr amrex::Real boundaryWeight = 10.0;
  const amrex::Vector<amrex::Geometry> geoms{amrex::Geometry(
    amrex::Box(amrex::IntVect(0), amrex::IntVect(nCell - 1)),
    amrex::RealBox({0.0, 0.0}, {1.0, 1.0}), 0, {0, 0})};

  // Merge and finalize clusters, as on the IO rank
  auto finalize = [&](sd::KeyedSurface a_clusters, int a_clusterSize) {
    Decimated out;
    sd::mergeClusters(a_clusters, nNC, a_clusterSize);
    sd::finalizeClusters(
      a_clusters, geoms, a_clusterSize, nNC, out.nodes, out.elts);
    out.length = sd::surfaceMeasure(out.nodes, out.elts, nNC);
    return out;
  };
  auto decimate = [&](
                    const sd::KeyedSurface& a_p, int a_clusterSize,
                    amrex::Real a_boundaryWeight) {
    return finalize(
      sd::clusterSurface(a_p, nNC, a_clusterSize, a_boundaryWeight),
      a_clusterSize);
  };

  // Closed circle: fewer segments, length and radius kept, still a closed
  // manifold curve, fields averaged
  const amrex::Real r = 0.3;
  const sd::KeyedSurface circle = makeCircle(r, dx);
  const Decimated dec = decimate(circle, clusterSize, boundaryWeight);
  const int nDecNodes = static_cast<int>(dec.nodes.size()) / nNC;
  AMREX_ALWAYS_ASSERT(dec.elts.size() < circle.elts.size());
  AMREX_ALWAYS_ASSERT(
    std::abs(dec.length - 2.0 * M_PI * r) < 1.0e-2 * 2.0 * M_PI * r);
  std::vector<int> valence(nDecNodes, 0);
  for (const int id : dec.elts) {
    ++valence[id];
  }
  for (int i = 0; i < nDecNodes; ++i) {
    AMREX_ALWAYS_ASSERT(valence[i] == 2);
    const amrex::Real* node = &dec.nodes[i * nNC];
    AMREX_ALWAYS_ASSERT(
      std::abs(std::hypot(node[0] - 0.5, node[1] - 0.5) - r) < 0.5 * dx);
    AMREX_ALWAYS_ASSERT(std::abs(node[2] - fieldVal) < 1.0e-12);
  }

  // Same circle split in two pieces clustered separately, as on two ranks,
  // the node at each cut existing in both: the merged result must match the
  // single piece
  const int nSeg = static_cast<int>(circle.elts.size()) / 2;
  const sd::KeyedSurface pieceA = subPolyline(circle, 0, nSeg / 2);
  const sd::KeyedSurface pieceB = subPolyline(circle, nSeg / 2, nSeg);
  sd::KeyedSurface split =
    sd::clusterSurface(pieceA, nNC, clusterSize, boundaryWeight);
  appendClusters(
    split, sd::clusterSurface(pieceB, nNC, clusterSize, boundaryWeight));
  const Decimated decSplit = finalize(split, clusterSize);
  AMREX_ALWAYS_ASSERT(decSplit.elts == dec.elts);
  AMREX_ALWAYS_ASSERT(decSplit.nodes.size() == dec.nodes.size());
  for (std::size_t n = 0; n < dec.nodes.size(); ++n) {
    AMREX_ALWAYS_ASSERT(std::abs(decSplit.nodes[n] - dec.nodes[n]) < 1.0e-12);
  }

  // Open straight line: vertices stay on the line, and the boundary quadrics
  // keep the end points, hence the length
  const amrex::Real x0 = 0.1, y0 = 0.2, x1 = 0.9, y1 = 0.6;
  const sd::KeyedSurface line = makeLine(x0, y0, x1, y1, dx);
  const Decimated decLine = decimate(line, clusterSize, boundaryWeight);
  AMREX_ALWAYS_ASSERT(decLine.elts.size() < line.elts.size());
  AMREX_ALWAYS_ASSERT(
    std::abs(decLine.length - sd::surfaceMeasure(line.vals, line.elts, nNC)) <
    1.0e-12);
  for (std::size_t i = 0; i < decLine.nodes.size(); i += nNC) {
    const amrex::Real dist = (y1 - y0) * (decLine.nodes[i] - x0) -
                             (x1 - x0) * (decLine.nodes[i + 1] - y0);
    AMREX_ALWAYS_ASSERT(std::abs(dist) < 1.0e-12);
  }

  // Line starting 0.1 dx past the first cell center of a cluster of 4^2
  // cells: its first two nodes share that cluster. Without boundary
  // quadrics, the cluster vertex is their mean and the line end retreats; with
  // them, the end point is kept.
  const sd::KeyedSurface line4 = makeLine(8.6 * dx, 8.6 * dx, x1, y1, dx);
  const amrex::Real lenKept = decimate(line4, 4, boundaryWeight).length;
  const amrex::Real lenFree = decimate(line4, 4, 0.0).length;
  AMREX_ALWAYS_ASSERT(
    std::abs(lenKept - sd::surfaceMeasure(line4.vals, line4.elts, nNC)) <
    1.0e-12);
  AMREX_ALWAYS_ASSERT(lenFree < lenKept - 0.5 * dx);
}

} // namespace
#endif

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  // Test table interpolation
  pele::physics::TabulatedFunctionParams tabfunc_par;
  tabfunc_par.host_only_parm().parm_parse_prefix = "test";
  tabfunc_par.initialize();
  pele::physics::TabulatedFunction tabfunc(
    static_cast<const pele::physics::TabulatedFunctionData*>(
      &tabfunc_par.host_parm()));
  amrex::Real outval;
  const amrex::GpuArray<amrex::Real, 2> invals{1.0, 1.0};
  const std::string varname = "T";
  const int idxT =
    pele::physics::get_var_index(varname.c_str(), &tabfunc_par.host_parm());
  tabfunc.get_value(idxT, invals.data(), outval);
  AMREX_ALWAYS_ASSERT(amrex::almostEqual(outval, 3.60));
  tabfunc_par.deallocate();
  amrex::Print() << "TabulatedFunction Tests passed" << std::endl;

  // Test Unit Conversions
  AMREX_ALWAYS_ASSERT(m2c::Length(1.0) == 100.0);
  AMREX_ALWAYS_ASSERT(c2m::Length(1.0, 2) == 1.0e-4);
  AMREX_ALWAYS_ASSERT(m2c::Mass(c2m::Mass(1.0)) == 1.0);
  AMREX_ALWAYS_ASSERT(m2c::Mu(1.0) == m2c::Mass(1.0) / m2c::Length(1.0));
  amrex::Print() << "Unit Conversion Tests Passed" << std::endl;

  // Test locate function
  amrex::Array<amrex::Real, 5> sorted_data{0.0, 0.2, 0.5, 0.9, 1.0};
  int length = sorted_data.size();
  int idxlo = -1;
  amrex::Real x = -50.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 0);
  x = 1.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == length - 1);
  x = 0.25;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 1);
  amrex::Print() << "Locate Tests Passed" << std::endl;

  // Test rectangle circle intersection area
  // Quarter circle intersection
  amrex::Real x0 = 0.0;
  amrex::Real x1 = 1.0;
  amrex::Real y0 = 0.0;
  amrex::Real y1 = 1.0;
  amrex::Real cx = 1.0;
  amrex::Real cy = 1.0;
  amrex::Real r = 1.0;
  amrex::Real A1 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A1 - 0.25 * M_PI) < 1e-14);
  // No intersection
  cx = 3.0;
  cy = 3.0;
  amrex::Real A2 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A2) < 1e-14);
  // Square covered by circle
  r = 20.0;
  amrex::Real A3 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A3 - 1.0) < 1e-14);
  amrex::Print() << "Intersection Area Tests Passed" << std::endl;

#if AMREX_SPACEDIM == 2
  testSurfaceDecimation();
  amrex::Print() << "Surface Decimation Tests Passed" << std::endl;
#endif

  amrex::Finalize();
  return 0;
}
