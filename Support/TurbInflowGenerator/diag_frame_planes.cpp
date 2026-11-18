#include <main.H>

#include <sstream>

using namespace amrex;

namespace {
// Plotfile components varidx of one frame, filled into one MultiFab
// (component d holds plotfile component varidx[d]): one fillPatchFromPlt()
// per run of consecutive plotfile components.
MultiFab
fill_frame(
  pele::physics::pltfilemanager::PltFileManager& pltfile,
  const int save_level,
  const Vector<int>& varidx,
  const Box& domain,
  const RealBox& probDomain,
  const int coord_sys,
  const Array<int, AMREX_SPACEDIM>& perio)
{
  // One Fab, no ghost cells
  BoxArray bxs{domain};
  DistributionMapping dm{bxs, 1};
  const int nvar = static_cast<int>(varidx.size());
  MultiFab mf{bxs, dm, nvar, {0, 0, 0}};
  Geometry geom{domain, probDomain, coord_sys, perio};
  for (int d = 0; d < nvar;) {
    int n = 1;
    while (d + n < nvar && varidx[d + n] == varidx[d] + n) {
      ++n;
    }
    pltfile.fillPatchFromPlt(save_level, geom, varidx[d], d, n, mf);
    d += n;
  }
  return mf;
}

// Component d of frame k (plotfile component comp) as one plane, serialized
// as a FAB.  A plain function rather than a lambda so that the device lambda
// inside it is not nested in another lambda.
std::string
frame_plane(
  const MultiFab& frame,
  const int d,
  const int comp,
  const int k,
  const Box& domain,
  const IntVect& bg,
  const IntVect& dim_map,
  const IntVect& periodicity)
{
  Box tbx{{0, 0, 0}, {bg[0] + 2, bg[1] + 2, 0}};
  FArrayBox TMP(tbx, 1);
  TMP.setVal<RunOn::Device>(1.2e34);

  const auto& pltarr = frame[0].const_array(d);
  const auto& planearr = TMP.array();

  ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int kp) noexcept {
    const IntVect idx_plane{i, j, kp};
    const IntVect idx_plt =
      index_mapper(idx_plane, dim_map, domain, periodicity);
    planearr(idx_plane) = pltarr(idx_plt);
  });
  std::ostringstream os;
  TMP.writeOn(os);
  amrex::Print() << "for comp = " << comp << " and k = " << k
                 << ": minval = " << TMP.min<RunOn::Device>(0)
                 << ", maxval = " << TMP.max<RunOn::Device>(0) << std::endl;
  return os.str();
}
} // namespace

void
turbinflow_from_diag_frame_planes(std::ofstream& ifsd, std::ofstream& ifsh)
{

  ParmParse pp;

  int nf = pp.countval("ifiles");
  AMREX_ALWAYS_ASSERT(nf >= 4);
  Vector<std::string> ifiles(nf);
  pp.getarr("ifiles", ifiles, 0, nf);

  int save_level = 0;
  pp.query("level", save_level);

  Array<int, AMREX_SPACEDIM> perio{1, 1, 1};
  int np = pp.countval("periodicity");
  if (np == 0) {
    amrex::Print() << "periodcity not specified, assuming fully periodic domain"
                   << std::endl;
  } else if (np != AMREX_SPACEDIM) {
    amrex::Error("Must provide 3 values for periodicity");
  }
  pp.query("periodicity", perio);

  int normal;
  pp.get("normal", normal);
  amrex::Vector<std::string> varnames{3};
  if (normal == 0) {
    varnames[0] = "y_velocity";
    varnames[1] = "z_velocity";
    varnames[2] = "x_velocity";
  } else if (normal == 1) {
    varnames[0] = "x_velocity";
    varnames[1] = "z_velocity";
    varnames[2] = "y_velocity";
  } else if (normal == 2) {
    varnames[0] = "x_velocity";
    varnames[1] = "y_velocity";
    varnames[2] = "z_velocity";
  } else {
    amrex::Abort("invalid normal");
  }
  // Optional scalar fields, stored after the velocity
  const auto scalars = requested_scalars();
  varnames.insert(varnames.end(), scalars.begin(), scalars.end());

  Box domain;
  RealBox probDomain;
  int coord_sys;
  amrex::Vector<std::string> varlist;
  {
    pele::physics::pltfilemanager::PltFileManager firstfile(ifiles[0]);
    int Nlevels = firstfile.getNlev();
    AMREX_ALWAYS_ASSERT(save_level >= 0 && save_level < Nlevels);
    domain = firstfile.getGeom(save_level).Domain();
    probDomain = firstfile.getGeom(save_level).ProbDomain();
    coord_sys = firstfile.getGeom(save_level).CoordInt();
    varlist = firstfile.getVariableList();
  }

  // Find the variables we need in the plt file, before writing anything
  const auto varidx = find_plt_vars(varnames, varlist);
  const int nvar = static_cast<int>(varidx.size());

  //
  // Write the first part of the header.
  //
  IntVect bg;
  bg[0] = domain.length(0);
  bg[1] = domain.length(1);
  bg[2] = nf;
  ifsh << bg[0] + 3 << ' ' << bg[1] + 3 << ' ' << bg[2] << '\n';

  Vector<Real> planeSize(2), dx(2);
  planeSize[0] = probDomain.length(0);
  planeSize[1] = probDomain.length(1);

  dx[0] = planeSize[0] / bg[0];
  dx[1] = planeSize[1] / bg[1];
  // Full precision: TurbInflow derives the file's spacing from these, and
  // for a mesh-mapped file they must agree with the trailer's Xi bounds.
  ifsh << std::setprecision(17) << planeSize[0] + 2 * dx[0] << ' '
       << planeSize[1] + 2 * dx[1] << ' ' << nf << '\n';

  ifsh << perio[0] << ' ' << perio[1] << ' ' << 0 << '\n';

  Vector<Real> fileTimes(nf);
  IntVect periodicity{perio[0], perio[1], 0};
  IntVect dim_map{0, 1, 2};

  // DAT holds the nf planes of component 0, then those of component 1, and
  // so on.  Each plotfile is read once, all its components filled into one
  // MultiFab, and each plane written straight to its place: all planes have the
  // same box, so the same size in DAT, and plane k of component d is at dat0 +
  // (d * nf + k) * plane_bytes.
  const std::streamoff dat0 = ifsd.tellp();
  std::streamoff plane_bytes = -1;
  Vector<std::streamoff> offsets(static_cast<long>(nvar) * nf);
  for (int k = 0; k < nf; ++k) {
    std::cout << "Loading frame " << k << " ... " << std::endl;
    pele::physics::pltfilemanager::PltFileManager pltfile(ifiles[k]);
    fileTimes[k] = pltfile.getTime();
    int Nlevels = pltfile.getNlev();
    AMREX_ALWAYS_ASSERT(save_level >= 0 && save_level < Nlevels);
    const MultiFab frame = fill_frame(
      pltfile, save_level, varidx, domain, probDomain, coord_sys, perio);
    for (int d = 0; d < nvar; ++d) {
      const std::string plane =
        frame_plane(frame, d, varidx[d], k, domain, bg, dim_map, periodicity);
      if (plane_bytes < 0) {
        plane_bytes = static_cast<std::streamoff>(plane.size());
      }
      AMREX_ALWAYS_ASSERT(
        static_cast<std::streamoff>(plane.size()) == plane_bytes);
      const long i = static_cast<long>(d) * nf + k;
      offsets[i] = dat0 + i * plane_bytes;
      ifsd.seekp(offsets[i]);
      ifsd.write(plane.data(), static_cast<std::streamsize>(plane.size()));
    }
    std::cout << "done" << std::endl;
  }
  if (!ifsd.good()) {
    amrex::Abort("Failed to write the planes to DAT");
  }

  // Plane offsets of the velocity, then the plane times
  for (int i = 0; i < AMREX_SPACEDIM * nf; ++i) {
    ifsh << offsets[i] << '\n';
  }
  for (int k = 0; k < nf; ++k) {
    ifsh << std::setprecision(17) << fileTimes[k] << std::endl;
  }

  // Optional mesh-mapping trailer; must come after the plane times.
  write_meshmap_trailer(ifsh, probDomain, dim_map);

  // Optional scalar fields, whose plane offsets end the SCALARS_V1 trailer
  if (!scalars.empty()) {
    write_scalars_header(ifsh, scalars);
    for (long i = static_cast<long>(AMREX_SPACEDIM) * nf; i < offsets.size();
         ++i) {
      ifsh << offsets[i] << '\n';
    }
  }
}
