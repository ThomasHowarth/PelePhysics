.. highlight:: rst

.. _sec:radiation:

*********
Radiation
*********

PelePhysics includes a radiative heat transfer module, located in ``Source/Radiation``, that computes the radiative source term of the energy equation for a participating gas and soot mixture. The module is header-only and all of its classes live in the ``PeleRad`` namespace. It consists of:

* A set of solvers for the P1 approximation of the radiative transfer equation, built on the AMReX linear solvers (``POneSingle``, ``POneSingleEB``, ``POneMulti``, ``POneMultiEB`` and ``POneMultiLevbyLev``)
* A gray Planck-mean absorption model for the gas phase and soot (``PlanckMean.H`` and ``SpectralModels.H``)
* A ``Radiation`` class that couples the two and is driven by the flow solver (``PeleLMRad.H``)

This section describes the model and its implementation, and the options that are inherited by the Pele codes that use it.

.. note:: For the flow solvers in the Pele suite, the radiation module is presently only used by PeleLMeX.

Radiation Equations
===================

The P1 approximation reduces the radiative transfer equation to a single elliptic equation for the incident radiation :math:`G`. For a gray, non-scattering medium with absorption coefficient :math:`\kappa`, the equation solved is

.. math::
   \kappa G - \nabla \cdot \left( \frac{1}{3 \kappa} \nabla G \right) = 4 \kappa \sigma T^4,

where :math:`\sigma` is the Stefan-Boltzmann constant and :math:`T` is the temperature. Once :math:`G` is known, the radiative source term of the energy equation is

.. math::
   \dot{q}_{\rm{rad}} = \kappa G - 4 \kappa \sigma T^4,

which is negative where the medium emits more than it absorbs.

The equation is discretized with the AMReX cell-centered linear operator ``MLABecLaplacian`` (``MLEBABecLap`` when embedded boundaries are used), which solves :math:`(A \alpha - B \nabla \cdot \beta \nabla) \phi = f`. The radiation module sets :math:`A = 1`, :math:`B = 1/3`, :math:`\alpha = \kappa`, :math:`\beta = 1/\kappa` and :math:`f = 4 \kappa \sigma T^4`. The cell-centered values of :math:`\beta` are averaged to the cell faces, and the system is solved with the AMReX geometric multigrid solver ``MLMG``.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

The boundary condition on each face of the domain is one of the AMReX linear solver types ``Dirichlet``, ``Neumann``, ``Robin`` or ``Periodic``, selected at runtime (see :ref:`sec:radiation:inputs`). Robin boundaries have the form

.. math::
   a G + b \frac{\partial G}{\partial n} = f,

where the coefficients :math:`a`, :math:`b` and :math:`f` are stored in the ghost cells outside the domain. The ``Radiation`` class of ``PeleLMRad.H`` fills them with :math:`a = -\kappa_b`, :math:`b = -2/3` and :math:`f = 0`, where :math:`\kappa_b` is the absorption coefficient of the adjacent interior cell.

With embedded boundaries, the condition applied on the EB surface is selected by ``ebbc_type``. If ``ebbc_type = 1``, a Dirichlet condition is applied on the EB with the values taken from the current solution field; for any other value no condition is set on the EB, which corresponds to a homogeneous Neumann condition.

Solver classes
~~~~~~~~~~~~~~

All solver classes are constructed from an ``MLMGParam`` object, the geometry, grids and distribution mapping, and references to the solution, right-hand side, :math:`\alpha`, :math:`\beta` and Robin coefficient ``MultiFab`` data, and all provide a ``solve()`` function. The solution ``MultiFab`` is used as the initial guess and, in its ghost cells, as the boundary value for Dirichlet boundaries.

* ``POneSingle``: single-level solver. It also provides ``calcRadSource`` to evaluate the radiative source term.
* ``POneSingleEB``: single-level solver with embedded boundaries.
* ``POneMulti``: multi-level solver that performs a composite solve over all AMR levels. The bottom solver is BiCGStab, or Hypre if ``use_hypre`` is set.
* ``POneMultiEB``: multi-level composite solver with embedded boundaries.
* ``POneMultiLevbyLev``: multi-level solver that solves each AMR level separately, from the coarsest to the finest, using the solution on the next coarser level as the coarse/fine boundary condition. This solver is exercised by the unit tests but is not used by the ``Radiation`` class.

Absorption Coefficient
======================

The absorption coefficient is computed with a gray Planck-mean model. The Planck-mean absorption coefficients of the individual species, :math:`\kappa_{p,i}(T)`, are tabulated as a function of temperature and are read from data files at the start of the simulation by the ``PlanckMean`` class. The functions that evaluate the absorption coefficient from these tables are located in the ``PeleRad::RadProp`` namespace in ``SpectralModels.H``.

The gas phase contribution is evaluated by ``getRadPropGas`` and includes ``CO2``, ``H2O`` and ``CO``:

.. math::
   \kappa_{\rm{gas}} = \frac{p}{10^5} \left( x_{\rm{CO_2}} \kappa_{p,\rm{CO_2}} + x_{\rm{H_2O}} \kappa_{p,\rm{H_2O}} + x_{\rm{CO}} \kappa_{p,\rm{CO}} \right),

where :math:`p` is the pressure in Pa and :math:`x_i` is the amount of species :math:`i` supplied by the calling code. With these units, :math:`\kappa_{\rm{gas}}` is in :math:`{\rm cm}^{-1}`. If all three :math:`x_i` are below :math:`10^{-5}`, :math:`\kappa_{\rm{gas}}` is set to :math:`10^{-3}` :math:`{\rm cm}^{-1}`.

If soot is present, its contribution is added by ``getRadPropSoot``:

.. math::
   \kappa = \kappa_{\rm{gas}} + f_v \kappa_{p,\rm{soot}},

where :math:`f_v` is the soot volume fraction.

The tabulated coefficients are interpolated linearly in temperature. Below 300 K and above 2800 K, the values at the ends of the table are used.

Planck-mean data files
~~~~~~~~~~~~~~~~~~~~~~

``PlanckMean`` reads six files from the directory given by the ``kppath`` input: ``kpl_co2.dat``, ``kpl_h2o.dat``, ``kpl_co.dat``, ``kpl_ch4.dat``, ``kpl_c2h4.dat`` and ``kpl_soot.dat``. A copy of these files is provided in ``Testing/Exec/Radiation/kpDB``.

Each file contains two columns, the temperature in K and the Planck-mean absorption coefficient, and must have exactly 126 rows with temperatures from 300 K to 2800 K in steps of 20 K. The code aborts if a file is missing, or if its temperatures or number of rows do not match this layout.

.. note:: All six files must be present, but the tables for ``CH4`` and ``C2H4`` are not used in the present evaluation of the absorption coefficient.

Coupling with the Flow Solvers
==============================

The interface between the flow solver and the P1 solvers is the ``PeleRad::Radiation`` class defined in ``PeleLMRad.H``. It owns the ``PlanckMean`` data, the P1 solver (``POneMulti``, or ``POneMultiEB`` if AMReX is built with EB support) and, for each AMR level, the ``MultiFab`` data of the solution, right-hand side, coefficients and Robin boundary coefficients. The flow solver uses it as follows:

#. Construct the ``Radiation`` object with the geometry, grids and distribution mapping of all levels, a ``RadComps`` structure holding the indices of the radiating species in the mechanism, and an ``amrex::ParmParse`` object from which the solver options are read. With EB, the vector of EB factories is passed as well. The constructor loads the Planck-mean data files.

#. For each level and each ``MFIter``, call ``updateSpecProp`` with the arrays of ``CO2``, ``H2O`` and ``CO``, the temperature, the pressure and, if soot is enabled, the soot volume fraction. This computes :math:`\kappa`, the coefficients :math:`\alpha` and :math:`\beta`, the emission term :math:`4 \kappa \sigma T^4` and the Robin boundary coefficients.

#. Call ``evaluateRad`` to solve for :math:`G` on all levels.

#. For each level and each ``MFIter``, call ``calcRadSource`` to add :math:`\kappa G - 4 \kappa \sigma T^4` to the energy source term.

The fields are also accessible for output through the ``G()``, ``kappa()`` and ``emis()`` functions.

``PeleLMRad.H`` expects its inputs in SI units (temperature in K, pressure in Pa) and returns the source term in :math:`{\rm W/m}^3`. Internally, the absorption coefficient is converted to :math:`{\rm m}^{-1}` and limited to a minimum of 0.01 :math:`{\rm m}^{-1}` in the interior of the domain; a minimum of 1 :math:`{\rm m}^{-1}` is used for the absorption coefficient :math:`\kappa_b` in the Robin boundary coefficients.

The soot contribution is included if ``PELELM_USE_SOOT`` is defined at compile time, in which case ``updateSpecProp`` takes the soot volume fraction as an additional argument.

.. note:: ``PeleCRad.H`` contains an earlier single-level interface in CGS units intended for PeleC. It has not been kept up to date with the solver classes: it does not match the present ``POneSingle`` constructor and does not read the ``kppath`` input. It defines the same class name and include guard as ``PeleLMRad.H``, so only one of the two headers can be included.

Radiation Flags and Inputs
==========================

* With GNUmake, the radiation sources are added to the build by including ``Source/Radiation/Make.package``. In PeleLMeX, this is done by specifying ``USE_RADIATION = TRUE`` in the ``GNUmakefile``, and the radiation solve is turned on in the input file using ``peleLM.do_rad_solve = 1``.

* The chemistry model specified with ``Chemistry_Model`` should contain the radiating species ``CO2``, ``H2O`` and ``CO``. The indices of these species are provided to the radiation module by the flow solver through the ``RadComps`` structure.

* To include the soot contribution, the soot model must also be enabled (see :ref:`Soot`).

* If ``use_hypre`` is set, AMReX must be built with Hypre support.

.. _sec:radiation:inputs:

Input file options
~~~~~~~~~~~~~~~~~~

The solver options are read by the ``MLMGParam`` structure from the ``amrex::ParmParse`` object supplied by the calling code. PeleLMeX uses the ``pelerad`` prefix, which is also used in the listing below. The values shown are those set by the default constructor of ``MLMGParam``. ::

  pelerad.kppath               = [Must be User Defined]  # Path to the directory containing the Planck-mean data files
  pelerad.lo_bc                = [Must be User Defined]  # Boundary condition on the low side in each direction:
                                                         # Dirichlet, Neumann, Robin or Periodic
  pelerad.hi_bc                = [Must be User Defined]  # Boundary condition on the high side in each direction
  pelerad.ebbc_type            = 2                       # EB boundary condition: 1 for Dirichlet, otherwise homogeneous Neumann

  pelerad.verbose              = 0                       # Verbosity of the multigrid solver
  pelerad.bottom_verbose       = 0                       # Verbosity of the bottom solver
  pelerad.max_iter             = 50                      # Maximum number of multigrid iterations
  pelerad.max_fmg_iter         = [Must be User Defined]  # Maximum number of F-cycle iterations
  pelerad.reltol               = 1.e-4                   # Relative tolerance of the multigrid solver
  pelerad.abstol               = 1.e-4                   # Absolute tolerance of the multigrid solver
  pelerad.bottom_reltol        = 1.e-6                   # Relative tolerance of the bottom solver
  pelerad.bottom_abstol        = 1.e-6                   # Absolute tolerance of the bottom solver
  pelerad.linop_maxorder       = 3                       # Maximum order of the linear operator stencil at boundaries
  pelerad.max_coarsening_level = 20                      # Maximum number of multigrid coarsening levels
  pelerad.agglomeration        = 0                       # Use agglomeration in the multigrid solver
  pelerad.consolidation        = 1                       # Use consolidation in the multigrid solver
  pelerad.use_hypre            = 0                       # Use Hypre as the bottom solver

.. note:: The values above are only set by the default constructor of ``MLMGParam``. When the structure is constructed from a ``ParmParse`` object, as is done by the ``Radiation`` class, options that are not present in the input file are left uninitialized. All of the options above should therefore be specified explicitly in the input file.

``lo_bc`` and ``hi_bc`` take one entry per spatial dimension, for example ::

  pelerad.lo_bc = Robin Periodic Periodic
  pelerad.hi_bc = Robin Periodic Periodic

The options ``max_bottom_iter``, ``agg_grid_size``, ``con_grid_size``, ``composite_solve`` and ``fine_level_solve_only`` are also read from the input file but are not used by the solver classes. The single-level solvers (``POneSingle`` and ``POneSingleEB``) do not use ``linop_maxorder``, ``bottom_reltol``, ``bottom_abstol``, ``agglomeration`` and ``consolidation``.

Testing
=======

A set of tests of the P1 solvers is provided in ``Testing/Exec/Radiation``. The ``CMakeLists.txt`` file in this directory is written to be added to the CMake build of a flow solver that contains PelePhysics as a submodule, such as PeleLMeX, and registers the tests with CTest. The tests read their options without a prefix, from the input files located in ``Testing/Exec/Radiation/inputs``. In addition to the solver options described above, the tests read the mesh options of the ``AMRParam`` structure (``max_level``, ``ref_ratio``, ``n_cell``, ``max_grid_size`` and ``plot_file_name``).

* ``POneSingle`` and ``POneMulti`` solve a problem with a manufactured trigonometric solution on a single level and on multiple AMR levels, respectively, and check the L1 norm of the error with respect to the exact solution. ``POneMulti`` is only built in 3D.
* ``POneSingleAF`` and ``POneMultiAF`` solve for the radiation field of prescribed temperature, composition and soot volume fraction fields, with the absorption coefficient evaluated from the Planck-mean data files in ``kpDB``.
* ``POneSingleEB`` and ``POneMultiEB`` solve a manufactured problem inside a rotated box defined with embedded boundaries. These are built in place of the other tests when EB is enabled.

The ``POneMulti`` and ``POneMultiAF`` tests are additionally run with ``composite_solve = 0``, in which case they use the ``POneMultiLevbyLev`` solver.
