.. highlight:: rst

.. _sec:Utility:

*******
Utility
*******

In addition to routines for evaluating chemical reactions, transport properties, and equation of state functions, PelePhysics includes other shared utilities that are utilized by both PeleC and PeleLM(eX). These utilities include support for:

* Premixed Flame (``PMF``) initialization from precomputed 1D flame profiles
* Turbulent inflows (``TurbInflow``) on domain boundaries from saved turbulence data
* Forced Turbulence (``TurbForcing``) for maintained HIT
* Mesh mapping (``MeshMap``) evaluation for solvers and utilities that use non-uniform, mapped grids
* Plt file management (``PltFileManager``)
* Output of runtime ``Diagnostics``
* Filtering of data (``Filter``), typically for LES
* Generic tabulated and neural network functions (``BlackBoxFunction``), used by the manifold models
* Basic ``Utilities``, including unit conversions

This section provides relevant notes on using these utilities across the Pele family of codes. There may be additional subtleties based on the code-specific implementation of these capabilities, in which case it will be necessary to refer to the documentation of the code of interest.

Premixed Flame Initialization
=============================

Pre-computed profiles from 1D freely propagating premixed flames are used to initialize a wrinkled `flamesheet <https://Pele-Suite.github.io/PeleLMeX/manual/html/Tutorials_FlameSheet.html>`_ in PeleLMeX, among other problems. Right now, this capability is not used in PeleC, but similar code that accomplishes the same task using data files of the same format is applied in PeleC. The code has two parts, a data container defined in ``PMFData.{H,cpp}`` that loads and stores data from the pre-computed profile, and a function defined in ``PMF.H`` that, when provided this data structure and the bounds of a cell of interest, returns temperature, velocity, and mole fractions from that location.

This code has two runtime parameters that may be set in the input file: ::

  pmf.datafile = pmf.dat
  pmf.do_cellAverage = 1

The first parameter specifies the path to the PMF data file. This file contains a two-line header followed by whitespace-delimited data columns in the order: position (cm), temperature (K), velocity (cm/s), density(g/cm3), species mole fractions. Sample files are provided in the relevant PeleLMeX (and PeleC) examples, and the procedure to generate these files with a provided script is described below. The second parameter specifies whether the PMF code does a finite volume-style integral over the queried cell (``pmf.do_cellAverage = 1``) or whether the code finds an interpolated value at the midpoint of the queried cell (``pmf.do_cellAverage = 0``)

Generating a PMF file
~~~~~~~~~~~~~~~~~~~~~

The script ``cantera_pmf_generator.py`` solves a 1D unstrained premixed flame using `Cantera <https://doi.org/10.5281/zenodo.6387882>`_ and saves it in the appropriate file format for use by the Pele codes. To use this script, first follow the :ref:`CEPTR instructions <sec_ceptr_software>` for setting up ``poetry``. Then return to the ``Utility/PMF/`` directory and run the script: ::

  poetry -C ../../../Support/ceptr/ run python cantera_pmf_generator.py -m dodecane_lu -f NC12H26 -d 0.01 -o ./

This example computes a dodecane/air flame using the ``dodecane_lu`` mechanism on a 0.01 m domain, and saves the resulting file to the present directory. It is also possible to pass all arguments via an input file, for example ``pmf_input.txt``: ::

  -m
  dodecane_lu 
  -f
  NC12H26
  -d
  0.01
  -o
  ./  

Then running the script with: ::
  
  poetry -C ../../../Support/ceptr/ run python cantera_pmf_generator.py @pmf_input.txt

You may choose any mechanism included with PelePhysics, as the Cantera ``.yaml`` format is included for all mechanisms. If you choose a reduced mechanism with QSS species, the 1D solution will be obtained using the skeletal version of the mechanism without QSS species being remove, as Cantera does not support QSS by default. In this case, the script will also save a data file with QSS species removed (ending in ``-qssa-removed.dat``) that is suitable for use with the QSS mechanisms in Pele. A range of additional conditions may be specified at the command line. To see the full set of options, use the help feature: ::

  poetry -C ../../../Support/ceptr/ run python cantera_pmf_generator.py --help

Note that when running Cantera, you may need to adjust the domain size to be appropriate for your conditions in order for the solver to converge. At times, you may also need to adjust the solver tolerances and other parameters that are specified within the script.

An additional script is provided to allow plotting of the PMF solutions. This script is used as follows (if no variable index is specified, temperature is plotted by default): ::

  poetry -C ../../../Support/ceptr/ run python plotPMF.py <pmf-file> <variable-index>

.. _sec_turbinflow:

Turbulent Inflows
=================

PelePhysics supports the capability of the flow solvers to have spatially and temporally varying inflow conditions based on precomputed planes of
turbulence data (3 components of velocity fluctuations). To reduce memory requirements, a fixed number of planes in the third dimension are read
at a time and then replaced when exhausted. This number of planes can be specified at runtime to balance I/O and memory requirements. Multiple
turbulent inflow patches can be applied on any domain boundary or on multiple domain boundaries. If multiple patches overlap on the same face, the
data from each overlapping patch are superimposed and added together.

This PelePhysics utility provides the machinery to load the data files and interpolate in space and time onto rectangular boundary patches,
which may or may not cover an entire boundary. Additional code within PeleC and PeleLMeX is required to drive this functionality, and the
documentation for the
relevant code should be consulted to fully understand the necessary steps. Typically, the TurbInflow capability provides velocity fluctuations,
and the mean inflow velocity must be provided through another means. For each code, an example test case for the capability is provided
in `Exec/RegTests/TurbInflow`. Inputs are available as part of this utility to rescale the data as needed. If differently shaped inlet patches
are required, this must be done by masking undesired parts of the patch on the PeleC or PeleLMeX side of the implementation.

Generating a turbulence file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relevant data files are generated using the tools in ``Support/TurbInflowGenerator``
and usage instructions are available in the :ref:`documentation <sec_turbfile>` on these tools.

Input file options
~~~~~~~~~~~~~~~~~~

The input file options that drive this utility (and thus are inherited by PeleC and PeleLMeX) are
provided below. ::

  turbinflows=low high                                    # Names of injections (can provide any number)

  turbinflow.low.turb_file      = TurbFileHIT/TurbTEST    # Path to directory created in previous step
  turbinflow.low.dir            = 1                       # Boundary normal direction (0,1, or 2) for patch
  turbinflow.low.side           = "low"                   # Boundary side (low or high) for patch
  turbinflow.low.turb_scale_loc = 633.151                 # Factor by which to scale the spatial coordinate between the data file and simulation
  turbinflow.low.turb_scale_vel = 1.0                     # Factor by which to scale the velocity between the data file and simulation
  turbinflow.low.turb_center    = 0.005 0.005             # Center point where turbulence patch will be applied
  turbinflow.low.turb_conv_vel  = 5.                      # Velocity to move through the 3rd dimension to simulate time evolution
  turbinflow.low.turb_nplane    = 32                      # Number of planes to read and store at a time
  turbinflow.low.time_offset    = 0.0                     # Offset in time for reading through the 3rd dimension
  turbinflow.low.verbose        = 0                       # Verbosity level
  turbinflow.low.extrap_nonperiodic = 0                   # Allow interpolation near edges of inflow patch where stencil may touch ghost cells
  turbinflow.low.tile_periodic  = 0                       # Cover the entire inflow face by periodically repeating/tiling the inflow patch
  turbinflow.low.interp_type    = quadratic               # Either quadratic (default) or linear (required if there are nonperiodic directions)
  turbinflow.low.time_periodic  = 0                       # Specify if the turbinflow data should be repeated in time, only valid for time varying turbinflow of type "istimeplanes"

  turbinflow.high.turb_file      = TurbFileHIT/TurbTEST   # All same as above, but for second injection patch
  turbinflow.high.dir            = 1
  turbinflow.high.side           = "high"
  turbinflow.high.turb_scale_loc = 633.151
  turbinflow.high.turb_scale_vel = 1.0
  turbinflow.high.turb_center    = 0.005 0.005
  turbinflow.high.turb_conv_vel  = 5.
  turbinflow.high.turb_nplane    = 32
  turbinflow.high.time_offset    = 0.0006
  turbinflow.high.verbose        = 2

TurbInflow Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To facilitate application of inflows generated from precursor Pele simulations, the data in the inflow files is now always interpreted as
cell-centered, which is a departure from legacy implementations designed to use fully-periodic data originating from solvers with nodal data.
However, some aspects of the implementation are still influenced by the legacy implementation. A main consideration in the current implementation
is that an inflow generated from a precursor simulation (using either of the ``periodic_plt`` or ``diag_frame_plane`` generation options) can
be directly applied as an inflow to a second simulation on the same grid without interpolation or data shifting.

Each TurbInflow plane contains a valid box corresponding to the physical domain used to generate the inflow plane,
plus one ghost cell on the low side and two ghost cells on the high side. But the dimensions listed
in the ``HDR`` correspond to the size of valid domain in the tangential directions plus two grid cell sizes. The normal direction dimension is just
the valid domain size. The behavior of the utility depends on the periodicity of the data. If the data is periodic in the normal direction, it is
treated as spatial data, the normal direction is traversed based on the specified ``turb_conv_vel``, and the inflow data can be recycled to allow
arbitrarily long simulations. If the data is not periodic in the normal direction, a list of time stamps for each plane is provided at the end of
the ``HDR`` file. The inflow is only valid from the first time stamp through, but not including, the 2nd last time stamp (the final plane
is necessary for interpolation).

For tangential directions, periodicity influences the type of interpolation that may be applied. For periodic tangential dimensions, quadratic
interpolation with a three point stencil is used by default, but linear interpolation is also available. The quadratic stencil is asymmetric,
using one point to the left and two to the right. Periodic directions may be tiled to cover the full face of the inflow using the `tile_periodic`
option, but by default they cover any cell of the simulation that has its center within the valid domain of the inflow. For nonperiodic
directions, linear interpolation is required. Because interpolation the outer 1/2 cell perimeter of a nonperiodic valid domain would
require an interpolation stencil that uses a ghost cell, and ghost cells may not be meaningfully populated (they are currently filled by
first order extrapolation), by default an error is raised if the TurbInflow utility is asked to populate a cell center in this region. This
error can be disabled with the ``extrap_nonperiodic`` option. In general, for nonperiodic directions, it is recommended to generate an inflow
that matches the grid from the finest level of the simulation using the inflow. The figure below demonstrates the structure of the turb
inflow plane data and interpolation approaches for periodic and nonperiodic tangential directions.

.. note:: The TurbInflow capability was not designed with embedded boundaries in mind. It can be applied for simulations using EB, but care should
          be take. Inflows should not be generated from simulations where EBs intersect the inflow plane.

Non-uniform and mapped grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inflow data is always uniformly spaced *in some coordinate*: the ``HDR``
carries only ``npts`` and ``probsize`` per direction, and the coordinate-to-index
conversion in the interpolator is a single affine expression. Which coordinate
that is depends on how the file was made. Synthetic (``turb_box``) data and data
extracted from a uniform-mesh precursor are uniform in physical position. Data
extracted with ``diag_frame_planes`` or ``periodic_plt`` from a precursor that
ran with a coordinate mapping (for example PeleLMeX's ``geometry.mesh_mapping``)
is uniform in that run's *computational* (:math:`\xi`) coordinate, because the
plane files and plotfiles carry the computational geometry with physical
velocity values.

Such a file is tagged by an optional ``MESHMAP_V2`` trailer at the end of the
``HDR``, placed after everything a legacy reader consumes (the plane times, or
the plane offsets when there are no times) so that older readers never see it.
It has one line per transverse direction, in the file's own transverse order::

  MESHMAP_V2
  <kind> <p> <p2> <p3> <q> <xi_lo> <xi_hi>
  <kind> <p> <p2> <p3> <q> <xi_lo> <xi_hi>

``kind``/``p``/``p2``/``p3``/``q`` are the per-axis payload of
``pele::physics::MeshMapEvaluator`` (``Source/Utility/MeshMap``, the same
descriptor PeleLMeX uses for its mapping: 0 identity, 1 constant, 2 exponential
stretch, 3 tanh stretch, 4 interior stretch), and ``xi_lo``/``xi_hi`` are the
precursor's computational-domain bounds along that axis. The generators write
the trailer when their input carries the precursor's ``geometry.mesh_mapping``
block (copy those lines verbatim; see ``Support/TurbInflowGenerator/README.md``),
and write nothing otherwise. A file *without* the trailer is, by declaration,
uniform in physical position. (An earlier ``MESHMAP_V1`` form, ``kind p q xi_lo
xi_hi``, is still read.)

``TurbInflow`` handles both kinds of file with one rule: the physical position
of each target cell is converted to the file's own coordinate, and the
interpolation then proceeds unchanged. For a physically-uniform file the
conversion is the identity; for a mapped file it is the inverse of the file's
map, so the interpolation happens in the precursor's :math:`\xi`. This is the
standard mapped-grid argument -- for a smooth map the interpolant's error is
the same order in :math:`\Delta\xi` as it would be in :math:`\Delta x` -- and it
covers every combination of uniform or mapped file and uniform or mapped target.
When the target's map and computational grid coincide with the file's, the
conversion lands on file cell centres and the file is reproduced exactly (the
same-grid guarantee above, now also for mapped meshes). With ``tile_periodic``
the position is wrapped into the file's physical period before the inverse.

For a mapped file ``turb_center`` is optional: by default the file is placed
where the precursor had it. If given, it is in the file's :math:`\xi` units,
must lie within ``[xi_lo, xi_hi]``, and the equivalent physical position is
printed. The availability of this sampling support is advertised by the
``PELEPHYSICS_TURBINFLOW_SAMPLES_MESHMAP`` macro (the trailer reader alone by
``PELEPHYSICS_TURBINFLOW_HAS_MESHMAP_HDR``); ``TurbInflow::file_has_map()``
and ``file_map()`` expose the file's descriptor so a solver can report whether
its own map matches.

A run *consuming* the data may itself be on a non-uniform grid. The turbulence
file is indexed by physical position, so a solver whose AMReX grid is a uniform
computational grid carrying a coordinate mapping must not use the
``amrex::Geometry`` overload of ``TurbInflow::add_turb()``, which would sample
the file at computational rather than physical coordinates. Such solvers pass
the physical cell-centre positions of the injection face explicitly, using the
overload taking ``x_phys`` / ``y_phys`` vectors; the ordering of the two
transverse directions is given by ``TurbInflow::transverseDirs()``. The presence
of that overload is advertised by the ``PELEPHYSICS_TURBINFLOW_HAS_COORD_ADDTURB``
macro.

Whichever entry point is used, the physically meaningful constraint is the ratio
of the target grid's local spacing on the injection face to the file's spacing.
``TurbInflow::file_transverse_dx()`` reports the latter in case units (the mean
for a mapped file, whose physical spacing varies; ``file_transverse_dx_range()``
gives its extremes; both are also printed at ``verbose > 0``) so that a solver
can check the ratio: substantially above one and the file cannot fill the scales
the grid resolves; substantially below one and the injected field is aliased
onto the grid.

.. note:: Diagnostics that slice the precursor solution -- notably
          ``DiagFramePlane`` -- locate the requested ``center`` using the AMReX
          geometry, which for a mapped run is the uniform computational grid, so
          by default ``center`` is a *computational* coordinate. Set
          ``center_is_physical = 1`` to give a physical position instead; the
          diagnostic reads the run's ``geometry.mesh_mapping`` block and
          inverts the map (the conversion is printed, and recorded in the
          ``PlaneData`` records the requested and converted centers, but the full
          mapping descriptor enters the turbulence file only through the generator's
          input, as described above.

.. figure:: ./Visualization/TurbInflowData.png

.. _sec_turbforce:

Forced Turbulence
=================

PelePhysics supports the capability for maintained homogeneous isotropic turbulence.  A source term in the momentum equation injects energy at
large scales, which naturally cascades to form (maintained) homogeneous isotropic turbulence.

Input file options
~~~~~~~~~~~~~~~~~~

The input file options that drive this utility (and thus are inherited by PeleC and PeleLMeX) are
provided below. ::

  turbforce.v                   = 0                       # Verbosity level
  turbforce.urms                = [Must be User Defined]  # Target urms
  turbforce.time_offset         = 0.0                     # Offset of the forcing function.
  turbforce.hack_lz             = 0                       # Allow periodic reproduction in z.
  turbforce.force_scale_fudge   = 1.0                     # Used for fine scale tuning of the forcing function.

  turbforce.ff_factor           = 4                       # Fast force coarsening factor.
  turbforce.nmodes              = 4                       # Largest mode to force.
  turbforce.forcing_epsilon     = 0.1                     # Reduce amplitude of modes used for breaking symmetry.
  turbforce.spectrum_type       = 2                       # Shape of the spectrum of the forcing.
  turbforce.moderate_zero_modes = 1                       # Reduce the impact of any zero mode.
  turbforce.mode_start          = 0                       # Starting mode

For ff_factor, nmodes, forcing_epsilon, spectrum_type and moderate_zero_modes, it is best to leave these as defaults unless you are confident on
the consiquences.

Forced Turbulence Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The strategy is to inject energy into the large scales using a source term in the momentum equation; the turbulence then develops naturally through the Richardson-Kolmogorov cascade.  The original description is in `Aspden et al. <https://doi.org/10.2140/camcos.2008.3.103>`_; a few modifications (e.g. tweaking the amplitudes of the modes through the random seed, which reduced the temporal variation in u_rms) and explicitly using a divergence-free form (see `Almgren et al. <https://doi.org/10.1137/110829386>`_). We typically run in a cube, or 4:1 box; larger boxes can achieved using the hack_lz option, which uses a periodic reproduction of the forcing in the z direction.  As the source term is a superposition of a substantial number of fourier modes, the direct evaluation is computationally expensive.  The resulting field is smooth (as we're forcing the large scales), so for computational efficiency, the source term is evaluated on a coarser box (by ff_factor, usually 4), and interpolated onto the required resolution; this is substantially faster, and indistinguishable from the direct approach.


.. _sec_meshmap:

Mesh Mapping
============

AMReX grids are uniform, but a solver may lay a coordinate mapping :math:`X = \Phi(\xi)` over the uniform grid so that the
physical grid is stretched; in this case the AMReX geometry describes the *computational* (:math:`\xi`) coordinate rather than
the physical one. The ``MeshMap`` utility provides ``pele::physics::MeshMapEvaluator``, a single definition of the supported
mappings that is shared by the solver and by every utility that must agree with the solver on where a cell physically is. This
is the same descriptor that PeleLMeX uses for its ``geometry.mesh_mapping`` capability (PeleLMeX aliases this type), and within
PelePhysics it is used by ``TurbInflow`` (see :ref:`the turbulent inflow section <sec_turbinflow>`), by ``DiagFramePlane``
(see :ref:`sec_diagnostics`), and by the tools in ``Support/TurbInflowGenerator``. The utility is header-only and consists of
two files: ``MeshMapEvaluator.H``, which defines the evaluator, and ``MeshMapEvaluatorInputs.H``, which builds an evaluator
from the input file.

The mappings are single-block, diagonal, and axis-aligned: the physical coordinate along a given direction depends only on the
computational coordinate along that same direction. Defining :math:`s = (\xi - \xi_{lo})/L` as the fraction of the axis in
computational space, where :math:`\xi_{lo}` and :math:`L` are the lower bound and length of the computational domain along that
direction, the available kinds of mapping are:

* ``Identity`` (0): no mapping, :math:`X = \xi`. This is the default when no mapping is specified.
* ``Constant`` (1): a constant stretch factor :math:`f` per direction, :math:`X = \xi_{lo} + f (\xi - \xi_{lo})`.
* ``ExpStretch`` (2): one-sided exponential stretching along a single direction that clusters cells at a wall on either the low or the high
  end of the domain. For a wall on the low end, :math:`(X - \xi_{lo})/L = (e^{\beta s} - 1)/(e^{\beta} - 1)`; the high end case is the mirror image.
* ``TanhStretch`` (3): two-sided stretching that clusters cells at both ends of the domain,
  :math:`(X - \xi_{lo})/L = \frac{1}{2} \left( 1 + \tanh(\beta (2 s - 1)) / \tanh \beta \right)`, with one :math:`\beta` per direction.
* ``InteriorStretch`` (4): clusters cells around a point in the interior of the domain, with separate coarsening rates
  :math:`\beta_{lo}` and :math:`\beta_{hi}` towards the low and high ends of the domain. The local stretch factor is proportional to
  :math:`\cosh^2(\beta \eta)`, where :math:`\eta` is the normalized distance from the clustering point on the relevant side.

All of the stretched mappings other than ``Constant`` preserve the length of the domain, and all reduce to the identity as
:math:`\beta \rightarrow 0` (stretching parameters with a magnitude below :math:`10^{-8}` are treated as zero).

Input file options
~~~~~~~~~~~~~~~~~~

The function ``make_evaluator_from_inputs()`` builds a ``MeshMapEvaluator`` from the same input file block that PeleLMeX reads, so
the lines from a PeleLMeX input file can be copied verbatim into the input file of a utility. The mapping is selected with
``geometry.mesh_mapping``, and the parameters for each mapping are then read using the name of the mapping as a prefix. Only the
options for the selected mapping are relevant. If ``geometry.mesh_mapping`` is not present, an ``Identity`` evaluator is returned. ::

  geometry.mesh_mapping          = TanhStretchMap         # One of ConstantMap, ExpStretchMap, TanhStretchMap, or InteriorStretchMap

  ConstantMap.scaling_factor     = 1.0 1.0 1.0            # Stretch factor in each direction (must be positive)

  ExpStretchMap.direction        = 1                      # Direction to be stretched (0, 1, or 2)
  ExpStretchMap.wall_end         = lo                     # End of the domain where cells are clustered (lo or hi)
  ExpStretchMap.beta             = 0.0                    # Stretching parameter (must be non-negative)

  TanhStretchMap.beta            = 0.0 0.0 0.0            # Stretching parameter in each direction (must be non-negative)

  InteriorStretchMap.beta_lo     = 0.0 0.0 0.0            # Coarsening rate towards the low end in each direction (must be non-negative)
  InteriorStretchMap.beta_hi     = 0.0 0.0 0.0            # Coarsening rate towards the high end in each direction (must be non-negative)
  InteriorStretchMap.center_frac = 0.5 0.5 0.5            # Location of the clustering point as a fraction of the domain length
  InteriorStretchMap.center      = 0.0 0.0 0.0            # [OPT] Physical location of the clustering point, overrides center_frac

The values shown above are the defaults. For the ``InteriorStretchMap``, the clustering point must lie strictly within the domain in
any direction that is stretched; the ``ExpStretchMap`` should be used to cluster cells at a boundary.

MeshMap Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``MeshMapEvaluator`` is a trivially-copyable struct that holds the kind of mapping (``m_kind``) and a small set of
per-direction parameters (``m_p``, ``m_p2``, ``m_p3``, and ``m_q``, the meaning of which depends on the kind of mapping). It is
therefore passed by value into GPU kernels, and all of its member functions can be called on the device. The main member
functions are:

* ``x_phys_from_xi``: the physical coordinate that corresponds to a computational coordinate along a given direction
* ``xi_from_x_phys``: the inverse, which is used to convert a physical position (for example, the location of a sampling plane)
  into the coordinate that indexes the AMReX grid. The ``Constant``, ``ExpStretch``, and ``TanhStretch`` mappings are inverted in closed
  form, while the ``InteriorStretch`` mapping is inverted with a safeguarded Newton iteration. The physical position is clamped to the
  extent of the domain, so a caller that needs a periodic image must wrap the position into the domain before the call.
* ``x_phys_cc`` and ``x_phys_face``: the physical coordinate of a cell center or a cell face from its integer index
* ``dx_phys_cc``: the physical width of a cell, computed as the difference between the positions of its two faces
* ``domain_size_phys``: the physical size of the domain along a given direction

Each of these is available in a form that takes an ``amrex::GeometryData`` (the geometry of the computational grid), and most are also
available in a form that takes the bounds of the computational domain directly, for callers that hold no AMReX geometry, such as
a turbulence file reader describing the domain of the precursor simulation. The static functions in the struct (``exp_fac``,
``exp_offset_norm``, ``tanh_fac``, ``tanh_offset_norm``, ``interior_fac``, ``interior_offset_norm``, and their inverses) are the
underlying stretch factors and normalized offsets for each family of mapping, from which the functions above are constructed.

.. _sec_pltfilemanager:

Plt File Management
===================

The ``PltFileManager`` class, defined in the ``pele::physics::pltfilemanager`` namespace, reads the data from an AMReX plotfile and
interpolates it onto the grids of a running simulation. It is used by the routines that allow the Pele codes to be initialized
based on data from plt files, including cases where the resolution or the number of levels of the simulation differs from that of
the plt file, and by the tools in ``Support/TurbInflowGenerator``. This utility has no input file options of its own: the
options that trigger initialization from a plt file and select the variables to be used are defined in PeleC and PeleLMeX, and the
documentation for the relevant code should be consulted for these.

The plt file is specified by passing its path to the constructor. Only the metadata is loaded at this point, which consists of
the ``Header`` file of the plt file and the ``Cell_H`` file of each level. The following information is then available through the
member functions of the class:

* ``getVariableList()``: the names of the variables in the plt file
* ``getNlev()``: the number of levels
* ``getTime()`` and ``getNsteps()``: the time and the number of time steps at which the plt file was written
* ``getGrid(lev)`` and ``getGeom(lev)``: the ``BoxArray`` and ``Geometry`` on each level
* ``getRefRatio(lev)``: the refinement ratio between each level and the next finer one

The function ``fillPatchFromPlt()`` fills a ``MultiFab`` defined on one level of the simulation with data from the plt file: ::

   pele::physics::pltfilemanager::PltFileManager pltData(pltfile);
   ...
   // Fill ncomp components of a_mf, starting from dataComp, using the
   // components of the plt file starting from pltComp
   pltData.fillPatchFromPlt(lev, geom[lev], pltComp, dataComp, ncomp, a_mf);

The data itself is only read the first time this function is called, and it is kept in memory until the ``PltFileManager`` is destroyed.
The mapping between the variables in the plt file and the components of the solution must be established by the calling code
using the list of variable names. An optional final argument, ``interp_type``, selects the spatial interpolation: cell
conservative linear interpolation (``1``, the default) or piecewise constant interpolation (any other value).

PltFileManager Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The way in which a level is filled depends on its resolution relative to that of the levels in the plt file, determined from the
ratio of the number of cells across the domain:

* If level 0 of the simulation has the same resolution as level 0 of the plt file, the data from level 0 of the plt file is
  copied directly.
* If the level being filled is finer than level 0 of the plt file (which is always expected for the levels above level 0 of the
  simulation), the entire level is first filled by interpolating from level 0 of the plt file to ensure that there is data
  everywhere. The finer levels of the plt file are then visited from coarsest to finest, and the regions that each one covers
  are overwritten: a plt file level that is coarser than the level being filled is interpolated to the required resolution,
  and when a plt file level with the same resolution is reached its data is copied directly and the process stops.
* If level 0 of the simulation is coarser than level 0 of the plt file, the code aborts, because averaging down of the
  plt file data is not supported.

The following points should also be noted:

* The number of spatial dimensions of the plt file must match that of the executable, and the problem domain of the simulation must be
  contained within that of the plt file.
* Only the standard AMReX plotfile layout is supported, with the cell-centered data for each level stored in ``Level_N/Cell``.
* All variables in the plt file are loaded on all levels, even if only a subset is requested, which should be kept in mind
  when considering memory requirements for large plt files.
* Only the valid cells of the ``MultiFab`` are filled. Ghost cells must be filled afterwards by the calling code.
* The plt file does not carry information about boundary conditions. Where the interpolation stencil extends outside the
  domain, first order extrapolation is used in non-periodic directions; the dummy boundary function ``FillExtDirDummy`` defined in
  ``PltFileManagerBCFill.H`` is provided to the AMReX fill patch routines to satisfy their interface and does not
  modify the data.
* Embedded boundaries are not accounted for in the interpolation.

.. _sec_diagnostics:

Diagnostics
===========

Analysing the data a-posteriori can become extremely cumbersome when dealing with extreme datasets.
PelePhysics supports shared capability for diagnostics available at runtime during PeleC and PeleLMeX simulation
and more are under development.
Currently, the list of diagnostic contains:

* ``DiagFramePlane`` : extract a plane aligned in the 'x','y' or 'z' direction across the AMR hierarchy, writing
  a 2D plotfile compatible with Amrvis, Paraview or yt. Only available for 3D simulations.
* ``DiagPDF`` : extract the PDF of a given variable and write it to an ASCII file.
* ``DiagConditional`` : extract statistics (average and standard deviation, integral or sum) of a
  set of variables conditioned on the value of given variable and write it to an ASCII file.

When using `DiagPDF` or `DiagConditional`, it is possible to narrow down the diagnostic to a region of interest
by specifying a set of filters, defining a range of interest for a variable. Note also that for these two diagnostics,
fine-covered regions are masked. An arbitrary number of diagnostics can be named in a list with the ``diagnostics``
keyword, then each diagnostic name in the list becomes a keyword allowing the type and specific options for
that diagnostic to be specified.
The following provide examples for each diagnostic in PeleLMeX (in PeleC, all diagnostics would be prefixed with
`pelec` instead of `peleLM`:

::

    #--------------------------DIAGNOSTICS------------------------

    peleLM.diagnostics = xnormP condT pdfTest

    peleLM.xnormP.type = DiagFramePlane                             # Diagnostic type
    peleLM.xnormP.file = xNorm5mm                                   # Output file prefix
    peleLM.xnormP.normal = 0                                        # Plane normal (0, 1 or 2 for x, y or z)
    peleLM.xnormP.center = 0.005                                    # Coordinate in the normal direction (grid coordinate; Xi under a mesh mapping)
    peleLM.xnormP.center_is_physical = 0                            # [OPT, DEF=0] take `center` as a physical position and convert it through geometry.mesh_mapping
    peleLM.xnormP.int    = 5                                        # Frequency (as step #) for performing the diagnostic
    peleLM.xnormP.interpolation = Linear                            # [OPT, DEF=Quadratic] Interpolation type : Linear or Quadratic
    peleLM.xnormP.field_names = x_velocity mag_vort density         # List of variables outputted to the 2D pltfile
    peleLM.xnormP.n_files = 2                                       # [OPT, DEF="min(256,NProcs)"] Number of files to write per level
    peleLM.xnormP.dump_ghost_if_OOB = 1                             # [OPT, DEF=false] if the specified coordinate is out-of-bounds, a plane of ghost cells in that direction will be dumped (for debugging purposes). If false, an error is raised if the requested plane is OOB.
    peleLM.xNormP.dump_flat_3D_plotfile                             # instead of a 2D plotfile, dump a 2D plotfile with 1 cell in the z-direction

    peleLM.condT.type = DiagConditional                             # Diagnostic type
    peleLM.condT.file = condTest                                    # Output file prefix
    peleLM.condT.int  = 5                                           # Frequency (as step #) for performing the diagnostic
    peleLM.condT.filters = xHigh stoich                             # [OPT, DEF=None] List of filters
    peleLM.condT.xHigh.field_name = x                               # Filter field
    peleLM.condT.xHigh.value_greater = 0.006                        # Filter definition : value_greater, value_less, value_inrange
    peleLM.condT.stoich.field_name = mixture_fraction               # Filter field
    peleLM.condT.stoich.value_inrange = 0.053 0.055                 # Filter definition : value_greater, value_less, value_inrange
    peleLM.condT.conditional_type = Average                         # Conditional type : Average, Integral or Sum
    peleLM.condT.nBins = 50                                         # Number of bins for the conditioning variable
    peleLM.condT.condition_field_name = temp                        # Conditioning variable name
    peleLM.condT.field_names = HeatRelease I_R(CH4) I_R(H2)         # List of variables to be treated

    peleLM.pdfTest.type = DiagPDF                                   # Diagnostic type
    peleLM.pdfTest.file = PDFTest                                   # Output file prefix
    peleLM.pdfTest.int  = 5                                         # Frequency (as step #) for performing the diagnostic
    peleLM.pdfTest.filters = innerFlame                             # [OPT, DEF=None] List of filters
    peleLM.pdfTest.innerFlame.field_name = temp                     # Filter field
    peleLM.pdfTest.innerFlame.value_inrange = 450.0 1500.0          # Filter definition : value_greater, value_less, value_inrange
    peleLM.pdfTest.nBins = 50                                       # Number of bins for the PDF
    peleLM.pdfTest.normalized = 1                                   # [OPT, DEF=1] PDF is normalized (i.e. integral is unity) ?
    peleLM.pdfTest.volume_weighted = 1                              # [OPT, DEF=1] Computation of the PDF is volume weighted ?
    peleLM.pdfTest.range = 0.0 2.0                                  # [OPT, DEF=data min/max] Specify the range of the PDF
    peleLM.pdfTest.field_name = x_velocity                          # Variable of interest


Filter
======

A utility for filtering data stored in AMReX data structures. When initializing the ``Filter`` class, the filter type
and filter width to grid ratio are specified. A variety of filter types are supported:

* ``type = 0``: no filtering
* ``type = 1``: standard box filter
* ``type = 2``: standard Gaussian filter

We have also implemented a set of filters defined in Sagaut & Grohens (1999) Int. J. Num. Meth. Fluids:

* ``type = 3``: 3 point box filter approximation (Eq. 26)
* ``type = 4``: 5 point box filter approximation (Eq. 27)
* ``type = 5``: 3 point box filter optimized approximation (Table 1)
* ``type = 6``: 5 point box filter optimized approximation (Table 1)
* ``type = 7``: 3 point Gaussian filter approximation
* ``type = 8``: 5 point Gaussian filter approximation (Eq. 29)
* ``type = 9``: 3 point Gaussian filter optimized approximation (Table 1)
* ``type = 10``: 5 point Gaussian filter optimized approximation (Table 1)

.. warning:: This utility is not aware of EB or domain boundaries. If the filter stencil extends across these boundaries,
             the boundary cells are treated as if they are fluid cells.
             It is up to the user to ensure an adequate number of ghost cells in the arrays are appropriately populated,
             using the ``get_filter_ngrow()`` member function of the class to determine the required number of ghost cells.


Developing
~~~~~~~~~~

The weights for these filters are set in ``Filter.cpp``. To add a
filter type, one needs to add an enum to the ``filter_types`` and
define a corresponding ``set_NAME_weights`` function to be called at
initialization.

The application of a filter can be done on a Fab or MultiFab. The loop nesting
ordering was chosen to be performant on existing HPC architectures and
discussed in PeleC milestone reports. An example call to the filtering operation is

::

   les_filter = Filter(les_filter_type, les_filter_fgr);
   ...
   les_filter.apply_filter(bxtmp, flux[i], filtered_flux[i], Density, NUM_STATE);

The user must ensure that the correct number of grow cells is present in the Fab or MultiFab.

.. _sec_blackboxfunction:

BlackBoxFunction
================

The ``BlackBoxFunction`` utility provides a generic interface for evaluating a function that maps a set of inputs to a set
of named outputs based on data loaded from a file at runtime. Two implementations are available: multilinear interpolation
in a structured table (``TabulatedFunction``) and inference with a neural network (``NeuralNetFunction``). Within PelePhysics,
this capability is used by the ``Manifold`` equation of state and transport models, where the inputs are the manifold
parameters and the outputs are the thermochemical and transport quantities that are needed by the flow solvers. Tools for
generating tables and neural networks in the appropriate formats are available in the `CMLM <https://github.com/NREL/cmlm>`_
repository, and a few sample files that are used for testing are provided in ``Support/CMLM``.

All implementations derive from the abstract ``BlackBoxFunction`` class and provide the same GPU-callable functions,
where ``indata`` is the array of input values:

* ``get_value(ivar, indata, out)``: evaluate a single output variable with index ``ivar``
* ``get_values(nvar, ivar, indata, out)``: evaluate a list of ``nvar`` output variables with indices ``ivar``; zero is returned
  for any negative index
* ``get_all_values(indata, out)``: evaluate all output variables
* ``get_derivs(ivar, indata, derivs)``: evaluate the derivatives of an output variable with respect to each of the inputs
* ``model()``: return the type of the function (``BlackBoxModel::TABLE`` or ``BlackBoxModel::NEURAL_NET``)

The index of an output variable is obtained from its name on the host with the ``get_var_index()`` function. By default
this function aborts if the variable is not found, but it can optionally return a negative index instead.

Input file options
~~~~~~~~~~~~~~~~~~

The data for a function is loaded and stored by a parameter class, either ``TabulatedFunctionParams`` or
``NeuralNetFunctionParams``. The input file options are read with a prefix that is specified by the code making use of the
function, to which ``table`` or ``network`` is appended. If no prefix is specified, the options are simply read as ``table.*``
or ``network.*``. For example, the ``Manifold`` equation of state uses the prefix ``manifold``, for which the options are: ::

  manifold.table.filename        = peletable_2dim_10grid  # Path to the table file [Must be User Defined]
  manifold.table.v               = 2                      # Verbosity level

  manifold.network.filename      = fgm_net.pnn            # Path to the network file, must have the .pnn extension [Must be User Defined]
  manifold.network.info_filename = fgm_net_info.txt       # Path to the network metadata file [Must be User Defined]
  manifold.network.v             = 2                      # Verbosity level
  manifold.network.cmlm_net      = 0                      # Network has the layer structure generated by CMLM (see below)

Only the options that correspond to the type of function in use are read. With a verbosity level greater than 0, the name of each
file is printed when it is loaded, and with a verbosity level greater than 1 (the default), a summary of the contents is also printed,
including the names of the inputs and outputs. The option that selects between the table and the network (``manifold.model`` in the
case of the ``Manifold`` equation of state) belongs to the code making use of the function, as described below.

Tabulated functions
~~~~~~~~~~~~~~~~~~~

The table file is a binary file that contains, in order:

#. The number of dimensions (inputs) of the table (1 integer)
#. The name of each dimension (64 characters each)
#. The number of grid points in each dimension (1 integer each); at least two grid points are required in every dimension
#. The number of output variables (1 integer)
#. The grid point locations for each dimension in turn (reals)
#. The name of the model (64 characters)
#. The name of each output variable (64 characters each)
#. The data for each output variable in turn (reals), with the index of the first dimension varying fastest

The integers and reals are read in the native binary representation of the executable, so the reals must match the precision
with which the code is compiled.

The grid for each dimension is defined independently and need not be uniformly spaced. Values are obtained by multilinear
interpolation, and inputs outside the bounds of the table are clipped to the bounds in each dimension. Derivatives are computed by
finite differences across the cell of the table that contains the input location, and are therefore piecewise constant in the direction of
differentiation. The number of dimensions of the table may be provided at compile time as a template argument, ``TabulatedFunction<N>``,
which allows the interpolation to be fully inlined. With the default template argument, the number of dimensions is determined
at runtime and tables of up to 5 dimensions are supported.

Neural network functions
~~~~~~~~~~~~~~~~~~~~~~~~

The neural network implementation is self-contained and does not rely on any external machine learning library. A network
is described by two files. The first is a binary file with the ``.pnn`` extension that contains the network itself. It
starts with four integers: the length of the strings in the file, the size in bytes of the floating point data, the size in bytes
of the integer data, and the number of inputs to the network. This is followed by a record for each layer consisting of the
name of the layer type, a format string that specifies the type and dimensions of each array of parameters for the layer,
and the parameters themselves. Three types of layer are supported, named following the PyTorch conventions:

* ``Linear``: a fully connected layer defined by a weight matrix and a bias vector
* ``LeakyReLU``: a leaky rectified linear unit activation function defined by its negative slope
* ``BatchNorm1d``: a batch normalization layer in inference mode, defined by a scaling vector and a bias vector that are applied
  element by element

The second file is a text file containing metadata in a ``key = value`` format. The ``=`` must be surrounded by spaces, multiple
values are separated by spaces, a line ending in a ``\`` (separated by a space from the last value) is continued on the next line,
and lines beginning with ``#`` followed by a space are comments.
The following keys are required: ::

  # Name of the neural network model
  model_name = premixed
  # Number of input dimensions
  ndim = 2
  # Number of output dimensions
  nvar = 16
  # Number of manifold parameters
  nmanpar = 2
  # Names of input variables
  dimnames = xi0 xi1
  # Names of output variables
  varnames = Y-CO2 Y-CO Y-OH SRC_H2O SRC_H2 SRC_CO2 SRC_CO SRC_CH4 SRC_O2 T RHO lnRHO DIFF VISC \
             SRC_xi0 SRC_xi1

Any other keys in the file are ignored by this utility, although they may be used elsewhere, as in the complete example
``Support/CMLM/fgm_net_info.txt``.

.. note:: At present, the only network architecture that can be evaluated is the one generated by CMLM, which consists of
          seven layers: ``Linear``, ``LeakyReLU``, ``BatchNorm1d``, ``Linear``, ``LeakyReLU``, ``BatchNorm1d``, ``Linear``.
          The ``cmlm_net`` option must be set to 1 to use such a network, as it is in all the examples provided
          with PelePhysics. The code aborts for any other network.

Derivatives of the network outputs are computed using first order finite differences with a relative perturbation
of :math:`10^{-5}` applied to each input in turn (an absolute perturbation of :math:`10^{-5}` is used for an input that is exactly zero), so the cost of ``get_derivs()`` is
one evaluation of the network per input in addition to the base evaluation. The number of inputs must be less than 100.

BlackBoxFunction Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameter classes are registered with a factory under the names ``Table`` and ``NeuralNet``, which allows the type of function to be
selected with a string from the input file. Their ``initialize()`` function reads the data files on the host into memory
allocated from the AMReX pinned arena, and copies the data structure that holds the sizes and the pointers to this data
(``TabulatedFunctionData`` or ``NeuralNetFunctionData``, both of which derive from ``BlackBoxFunctionData``) to the device.
The function objects themselves are lightweight and do not own any data, so they are constructed from a pointer to this data structure
wherever they are needed, including within GPU kernels. The ``deallocate()`` function must be called to free the memory once the function is
no longer needed. The following is an example based on the unit test in ``Testing/Exec/UtilityUnitTest``, which reads the table
specified with ``test.table.filename``: ::

   // Load the table based on the input file options test.table.*
   auto tabfunc_par = pele::physics::PeleParamsGeneric<
     pele::physics::BlackBoxFunctionData>::create("Table");
   tabfunc_par->host_only_parm().parm_parse_prefix = "test";
   tabfunc_par->initialize();

   // Look up the index of an output variable on the host
   const int idxT = pele::physics::get_var_index("T", &tabfunc_par->host_parm());

   // Evaluate (host_parm() on the host, device_parm() in a GPU kernel)
   pele::physics::TabulatedFunction<> tabfunc(
     static_cast<const pele::physics::TabulatedFunctionData*>(
       &tabfunc_par->host_parm()));
   amrex::Real outval;
   const amrex::GpuArray<amrex::Real, 2> invals{1.0, 1.0};
   tabfunc.get_value(idxT, invals.data(), outval);

   tabfunc_par->deallocate();

The ``BlackBoxFunctionFactory`` struct, defined in ``BlackBoxFunctionFactory.H``, constructs a function object from a pointer to a ``BlackBoxFunctionData`` and
provides access to it through ``get_func()``. Its template argument determines whether the type of function is fixed at compile time
or selected at runtime:

* ``BlackBoxFunctionFactory<TabulatedFunction<N>>``: a table with ``N`` dimensions
* ``BlackBoxFunctionFactory<NeuralNetFunction>``: a neural network
* ``BlackBoxFunctionFactory<BlackBoxFunction>``: either a table or a neural network, depending on the data loaded at runtime.
  This relies on virtual function calls and a table with a number of dimensions determined at runtime, and is therefore
  slower than the other two options.

For the first two options, the code aborts if the type of the data that was loaded at runtime does not match the type selected at compile time.

Utilities
=============================
The ``utilities`` namespace contains some functions that may be useful to users when writing problem-specific source code.

This includes unit conversions, which are particularly useful for PeleLM(eX) users working with mixed unit systems. The following aliases, defined in a file header, enable straightforward conversions between MKS and CGS units for use in equation of state (``eos``) function calls:

::

   namespace m2c = pele::physics::utilities::mks2cgs;
   namespace c2m = pele::physics::utilities::cgs2mks;

For example, users can call ``eos.PYT2R`` as follows:

::

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real P_mean = 101325.0_rt; // Pressure in Pa (MKS)
   amrex::Real massfrac[NUM_SPECIES]; // Mass fractions tracked by PeleLM(eX)
   amrex::Real Temp = 300.0_rt; // Temp in K
   amrex::Real rho_cgs = 0.0_rt; // Density in g/cm^3 (CGS)

   // Calculate density in CGS units, then convert to MKS
   eos.PYT2R(m2c::P(P_mean), massfrac, Temp, rho_cgs);
   amrex::Real rho = c2m::Rho(rho_cgs); // Convert eos density to MKS

The ``utilities`` namespace also contains some other useful functions: ``locate`` to find the nearby indices of a value in an array for interpolation
and ``rectangle_circle_intersection_area`` to find the area of the intersections between rectangles and circles.
