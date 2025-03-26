#include <stateinflow.H>

namespace pele {
namespace physics {
namespace stateinflow {
void
StateInflow::init(amrex::Geometry const& /*geom*/)
{
  amrex::ParmParse ppr;

  int n_sp = 0;
  n_sp = ppr.countval("stateinflows");
  amrex::Vector<std::string> sp_list;
  if (n_sp > 0) {
    sp.resize(n_sp);
    sp_list.resize(n_sp);
    for (int n = 0; n < n_sp; n++) {
      ppr.get("stateinflows", sp_list[n], n);
    }
  }

  for (int n = 0; n < n_sp; n++) {

    amrex::ParmParse pp("stateinflow." + sp_list[n]);
    if (pp.countval("state_file") > 0) {

      // Query data
      pp.query("state_file", sp[n].m_state_file);
      sp[n].dir = -1;
      pp.query("dir", sp[n].dir);
      AMREX_ASSERT_WITH_MESSAGE(
        sp[n].dir >= 0 && sp[n].dir < AMREX_SPACEDIM,
        "Injection direction is needed: 0, 1 or 2");
      std::string side;
      pp.query("side", side);
      if (side == "low") {
        sp[n].side = amrex::Orientation::low;
      } else if (side == "high") {
        sp[n].side = amrex::Orientation::high;
      } else {
        amrex::Abort("stateinflow.side can only be low or high");
      }
      pp.query("time_offset", sp[n].time_shift);
      pp.query("state_scale_loc", sp[n].state_scale_loc);
      //pp.query("turb_scale_vel", tp[n].turb_scale_vel); //Do we need vel scaling? probably not
      amrex::Print() << "Initializing stateinflow " << sp_list[n]
                     << " with file " << sp[n].m_state_file
                     << " (location coordinates in will be scaled by "
                     << sp[n].state_scale_loc << ") \n";
      //<< " and velocity out to be scaled by "
      //               << tp[n].turb_scale_vel << ") \n";

      // Get the statecenter on the injection face
      amrex::Vector<amrex::Real> state_center(AMREX_SPACEDIM - 1, 0);
      pp.getarr("state_center", state_center);
      AMREX_ASSERT_WITH_MESSAGE(
        state_center.size() == AMREX_SPACEDIM - 1,
        "state_center must have AMREX_SPACEDIM-1 elements");
      for (int idim = 0; idim < state_center.size(); ++idim) {
        state_center[idim] *= sp[n].state_scale_loc;
      }

      pp.query("state_nplane", sp[n].nplane);
      AMREX_ASSERT(sp[n].nplane > 0);
      pp.query("state_conv_vel", sp[n].state_conv_vel);
      AMREX_ASSERT(sp[n].state_conv_vel > 0);

      // Set other stuff
      std::string state_header = sp[n].m_state_file + "/HDR";
      std::ifstream is(state_header.c_str());
      if (!is.is_open()) {
        amrex::Abort("Unable to open input file " + state_header);
      }
      amrex::Array<int, AMREX_SPACEDIM> npts = {{0}};
      amrex::Array<amrex::Real, AMREX_SPACEDIM> probsize = {{0}};
      amrex::Array<int, AMREX_SPACEDIM> iper = {{0}};
      is >> npts[0] >> npts[1] >> npts[2];
      is >> probsize[0] >> probsize[1] >> probsize[2];
      is >> iper[0] >> iper[1] >>
        iper[2]; // Unused - we assume it is always fully periodic

      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        sp[n].dx[idim] = probsize[idim] / amrex::Real(npts[idim] - 1);
        sp[n].dxinv[idim] = 1.0 / sp[n].dx[idim];
      }

      // The following is relative to the injection face:
      // 0 and 1 are transverse directions, 2 is normal
      // one ghost point on each side, tangential to inflow face
      sp[n].pboxsize[0] = probsize[0] - 2.0 * sp[n].dx[0];
      sp[n].pboxsize[1] = probsize[1] - 2.0 * sp[n].dx[1];
      sp[n].pboxsize[2] = probsize[2];

      sp[n].npboxcells[0] = npts[0] - 3;
      sp[n].npboxcells[1] = npts[1] - 3;
      sp[n].npboxcells[2] = npts[2];

      // Center the state
      sp[n].pboxlo[0] = state_center[0] - 0.5 * sp[n].pboxsize[0];
      sp[n].pboxlo[1] = state_center[1] - 0.5 * sp[n].pboxsize[1];
      sp[n].pboxlo[2] = 0.;

      amrex::Box sbx(
        amrex::IntVect(AMREX_D_DECL(1, 1, 1)),
        amrex::IntVect(AMREX_D_DECL(npts[0], npts[1], sp[n].nplane)));

      sp[n].sdata = new amrex::FArrayBox(sbx, NVAR, amrex::The_Async_Arena());

      sp[n].kmax = npts[2];

      amrex::Real rdummy;
      if (sp[n].isswirltype) {
        for (int i = 0; i < sp[n].kmax; i++) {
          is >> rdummy; // Time for each plane - unused at the moment
        }
      }

      // Offset for each plane in Binary StateFile
      sp[n].m_offset.resize(sp[n].kmax * AMREX_SPACEDIM + sp[n].kmax*(NVAR-AMREX_SPACEDIM)); 
      sp[n].offset = sp[n].m_offset.data();
      sp[n].offset_size = sp[n].m_offset.size();
      for (int i = 0; i < sp[n].offset_size; i++) {
        is >> sp[n].offset[i];
      }
      is.close();
    }
    stateinflow_initialized = true;
  }
}

void
StateInflow::add_state(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const int dir,
  const amrex::Orientation::Side& side,
  int ncomp)
{
  AMREX_ASSERT(stateinflow_initialized);

  // Box on which we will access data
  amrex::Box bvalsBox = bx;
  int planeLoc =
    (side == amrex::Orientation::low ? geom.Domain().smallEnd()[dir] - 1
                                     : geom.Domain().bigEnd()[dir] + 1);
  bvalsBox.setSmall(dir, planeLoc);
  bvalsBox.setBig(dir, planeLoc);

  // Define box that we will fill with : need to be z-normal
  // Get transverse directions
  int tdir1 = (dir != 0) ? 0 : 1;
  int tdir2 = (dir != 0) ? ((dir == 2) ? 1 : 2) : 2;
  int tr1Lo = bvalsBox.smallEnd()[tdir1];
  int tr1Hi = bvalsBox.bigEnd()[tdir1];
  int tr2Lo = bvalsBox.smallEnd()[tdir2];
  int tr2Hi = bvalsBox.bigEnd()[tdir2];
  const amrex::IntVect lo(AMREX_D_DECL(tr1Lo, tr2Lo, planeLoc));
  const amrex::IntVect hi(AMREX_D_DECL(tr1Hi, tr2Hi, planeLoc));
  amrex::Box stateBox(lo, hi);
  amrex::FArrayBox s(stateBox, NVAR, amrex::The_Async_Arena());
  s.setVal<amrex::RunOn::Device>(0);
  
  // Add state from all the sp acting on this face
  for (int n = 0; n < sp.size(); n++) {

    if (sp[n].dir == dir && sp[n].side == side) {

      // 0 and 1 are the two transverse directions
      amrex::Vector<amrex::Real> x(stateBox.size()[0]), y(stateBox.size()[1]);
      for (int i = stateBox.smallEnd()[0]; i <= stateBox.bigEnd()[0]; ++i) {
        x[i - stateBox.smallEnd()[0]] =
          (geom.ProbLo()[tdir1] + (i + 0.5) * geom.CellSize(tdir1)) *
          sp[n].state_scale_loc;
      }
      for (int j = stateBox.smallEnd()[1]; j <= stateBox.bigEnd()[1]; ++j) {
        y[j - stateBox.smallEnd()[1]] =
          (geom.ProbLo()[tdir2] + (j + 0.5) * geom.CellSize(tdir2)) *
          sp[n].state_scale_loc;
      }

      // Get the state
      amrex::Real z =
        (time + sp[n].time_shift) * sp[n].state_conv_vel * sp[n].state_scale_loc;
      fill_state_plane(sp[n], x, y, z, s);
    }
  }
  // Moving it into data
  set_state(dir, tdir1, tdir2, s, data, dcomp, ncomp);  //permutes vels if needed

#if 0
  std::string junk = "TurbV_AftTP"+std::to_string(n)+"_D";
  std::ofstream os;
  os.precision(15);
  os.open(junk.c_str());
  data.writeOn(os);
  os.close();
  amrex::Abort();
#endif
}

void
StateInflow::set_state(
  int normDir,
  int transDir1,
  int transDir2,
  amrex::FArrayBox& s,
  amrex::FArrayBox& data,
  const int dcomp, int ncomp)
{
  // copy velocity fluctuations from plane into data
  const auto& box = s.box(); // z-normal plane
  const auto& s_in = s.array();
  const auto& s_out = data.array(dcomp);

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // From z-normal box index to data box index
      
    int idx[3] = {0};
    idx[transDir1] = i;
    idx[transDir2] = j;
    idx[normDir] = k;
    if (ncomp == AMREX_SPACEDIM || ncomp == NVAR) { //vels or everything
      s_out(idx[0], idx[1], idx[2], transDir1) =
	s_in(i, j, k, 0); // transverse velocity 1
      s_out(idx[0], idx[1], idx[2], transDir2) =
	s_in(i, j, k, 1); // transverse velocity 2
      s_out(idx[0], idx[1], idx[2], normDir) =
	s_in(i, j, k, 2); // normal velocity
    }
    if (ncomp == 1 || ncomp == NVAR || ncomp == NUM_SPECIES) { //density, temp, species or everything
      for (int n = 0; n<ncomp;n++) {
	s_out(i,j,k,n) = s_in(i,j,k,n);
	//std::cout << s_out(i,j,k,n) << std::endl;
      }
    }
    
        
    });
}

void
StateInflow::read_one_state_plane(StateParm& a_sp, int iplane, int k)
{
  // There are NVAR * kmax planes of FABs.
  // The first component are in the first kmax planes,
  // the second component in the next kmax planes, ....
  // Note also that both (*plane) and (*ncomp) start from
  // 1 not 0 since they're passed from Fortran.

  std::string state_data = a_sp.m_state_file + "/DAT";
  std::ifstream ifs(state_data.c_str());
  if (!ifs.is_open()) {
    amrex::Abort("Unable to open input file " + state_data);
  }

  amrex::Box dstBox = a_sp.sdata->box();
  dstBox.setSmall(2, iplane);
  dstBox.setBig(2, iplane);
  for (int n = 0; n < NVAR; ++n) {
    const long offset_idx = (k + 1) + (n * a_sp.kmax);
    AMREX_ASSERT_WITH_MESSAGE(
      offset_idx < a_sp.offset_size, "Bad state fab offset idx");

    const long start = a_sp.offset[offset_idx];
    ifs.seekg(start, std::ios::beg);

    if (!ifs.good()) {
      amrex::Abort("getplane(): seekg() failed");
    }
    amrex::FArrayBox tmp;
    tmp.readFrom(ifs);
    amrex::Box srcBox = tmp.box();
#if 0 //we good here
    const auto array = tmp.array();
    amrex::ParallelFor(srcBox, [=] AMREX_GPU_DEVICE(int i, int j, int l) noexcept {
	std::cout << "(" << i << ", " << j << ", " << l << ", " << n << ") = " << array(i,j,l) <<std::endl; 
      });
#endif    
    a_sp.sdata->copy<amrex::RunOn::Device>(tmp, srcBox, 0, dstBox, n, 1); 
  }
#if 0  //good here
  const auto& sdataArr = a_sp.sdata->array();
  amrex::ParallelFor(dstBox, [=] AMREX_GPU_DEVICE(int i, int j, int l) noexcept {
      for (int n = 0; n<NVAR; n++) {
	std::cout << "(" << i << ", " << j << ", " << l << ", " << n << ") = " << sdataArr(i,j,l,n) <<std::endl;
      }
      });
#endif
  ifs.close();
}

void
StateInflow::read_state_planes(StateParm& a_sp,amrex::Real z)
{
  int izlo = (int)(round(z * a_sp.dxinv[2])) - 1;
  int izhi = izlo + a_sp.nplane - 1;
  a_sp.szlo = static_cast<amrex::Real>(izlo) * a_sp.dx[2];
  a_sp.szhi = static_cast<amrex::Real>(izhi) * a_sp.dx[2];

#if 0
  amrex::AllPrint() << "read_state_planes filling " << izlo << " to " << izhi
                 << " covering " << a_sp.szlo + 0.5 * a_sp.dx[2]
                 << " to "       << a_sp.szhi - 0.5 * a_sp.dx[2] << " for z = " << z << std::endl;
#endif

  for (int iplane = 1; iplane <= a_sp.nplane; ++iplane) {
    int k = (izlo + iplane - 1) % (a_sp.npboxcells[2] - 2);
    read_one_state_plane(a_sp, iplane, k);
  }
#if 0 //probably can delete
  int myproc = amrex::ParallelDescriptor::MyProc();
  std::string junk = "StateData_proc"+std::to_string(myproc)+"_D";
  std::ofstream os;
  os.precision(15);
  os.open(junk.c_str());
  a_sp.sdata->writeOn(os);
  os.close();
  //amrex::Abort();
#endif
}

void
StateInflow::fill_state_plane(
  StateParm& a_sp,
  const amrex::Vector<amrex::Real>& x,
  const amrex::Vector<amrex::Real>& y,
  amrex::Real z,
  amrex::FArrayBox& s)
{
  if (
    (z < a_sp.szlo + 0.5 * a_sp.dx[2]) || (z > a_sp.szhi - 0.5 * a_sp.dx[2])) {
#if 0
    {
      amrex::AllPrint() << "Reading new data because z " << z << " is outside " << a_sp.szlo + 0.5 * a_sp.dx[2] << " and "
                     << a_sp.szhi - 0.5 * a_sp.dx[2] << std::endl;
    }
#endif
    read_state_planes(a_sp, z);
  }

  const auto& bx = s.box();
  const auto& vd = s.array();

  amrex::Gpu::DeviceVector<amrex::Real> x_dev(x.size());
  amrex::Gpu::DeviceVector<amrex::Real> y_dev(y.size());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, x.begin(), x.end(), x_dev.begin());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, y.begin(), y.end(), y_dev.begin());
  amrex::Real* xd = x_dev.data();
  amrex::Real* yd = y_dev.data();

  amrex::Real velScale = (a_sp.side == amrex::Orientation::high) ? -1 : 1;//-a_sp.turb_scale_vel
  const auto& npboxcells = a_sp.npboxcells;
  const auto& pboxlo = a_sp.pboxlo;
  const auto& szlo = a_sp.szlo;
  const auto& dxinv = a_sp.dxinv;
  const auto& sd = a_sp.sdata->array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real cx[3], cy[3], cz[3], ydata[3];
    amrex::Real zdata[3][3];

    amrex::Real zz =
      (z - szlo) * dxinv[2];        // How many dz away from the left side ?
    int k0 = (int)(std::round(zz)); // What's the closest point ?
    zz -= amrex::Real(k0);
    cz[0] = 0.5 * (zz - 1.0) * (zz - 2.0); // Weight of k0 - 1
    cz[1] = zz * (2.0 - zz);               // Weight of k0
    cz[2] = 0.5 * zz * (zz - 1.0);         // Weight of k0 + 1
    k0 += 1;                               // Index starting at 1

    for (int n = 0; n < NVAR; ++n) {
      amrex::Real xx = (xd[i - bx.smallEnd(0)] - pboxlo[0]) * dxinv[0];
      amrex::Real yy = (yd[j - bx.smallEnd(1)] - pboxlo[1]) * dxinv[1];
      int i0 = (int)(std::round(xx));
      int j0 = (int)(std::round(yy));
      xx -= amrex::Real(i0);
      yy -= amrex::Real(j0);
      cx[0] = 0.5 * (xx - 1.0) * (xx - 2.0);
      cy[0] = 0.5 * (yy - 1.0) * (yy - 2.0);
      cx[1] = xx * (2.0 - xx);
      cy[1] = yy * (2.0 - yy);
      cx[2] = 0.5 * xx * (xx - 1.0);
      cy[2] = 0.5 * yy * (yy - 1.0);

      if (i0 >= 0 && i0 < npboxcells[0] && j0 >= 0 && j0 < npboxcells[1]) {
        i0 += 2;
        j0 += 2;
        for (int ii = 0; ii <= 2; ++ii) {
          for (int jj = 0; jj <= 2; ++jj) {
            zdata[ii][jj] = cz[0] * sd(i0 + ii, j0 + jj, k0 - 1, n) +
                            cz[1] * sd(i0 + ii, j0 + jj, k0, n) +
                            cz[2] * sd(i0 + ii, j0 + jj, k0 + 1, n);
          }
        }
        for (int ii = 0; ii <= 2; ++ii) {
          ydata[ii] =
            cy[0] * zdata[ii][0] + cy[1] * zdata[ii][1] + cy[2] * zdata[ii][2];
        }
        vd(i, j, k, n) = cx[0] * ydata[0] + cx[1] * ydata[1] + cx[2] * ydata[2];
	if (n == VELX || n == VELY || n == VELZ) {
	  vd(i, j, k, n) *= velScale;
	}
      }
    }
  });
  
  amrex::Gpu::synchronize(); // Ensure that DeviceVector's don't leave scope
                             // early
#if 0
  const auto& bxc = s.box();
  const auto& vdc = s.array();
  amrex::ParallelFor(bxc, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if(k<=0) {
	for (int n =0; n<NVAR;n++) {
	  std::cout << vdc(i,j,k,n) << std::endl;
	}
      }
    });
#endif
  
}
} // namespace stateinflow
} // namespace physics
} // namespace pele
