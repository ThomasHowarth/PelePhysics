.. highlight:: rst

.. _sec:testing:

*******
Testing
*******

PelePhysics includes several short stand-alone programs that exercise its various capabilities without requiring one of the Pele flow solvers. These programs are located in the ``Testing/Exec`` directory:

* ``EosEval``: evaluation of equation of state functions
* ``TranEval``: evaluation of transport coefficients
* ``ReactEval``: integration of chemical reactions on a grid of cells
* ``IgnitionDelay``: ignition delay time of a 0D reactor
* ``Jacobian``: comparison of the analytical chemical Jacobian with finite differences
* ``PMFReader``: initialization of a field from a premixed flame (PMF) data file
* ``Radiation``: P1 radiation solver tests
* ``CeptrTest``: comparison of a 0D reactor against Cantera
* ``UtilityUnitTest``: unit tests of the basic utilities

These programs are also what is run by the continuous integration (CI) of the repository, as described in :ref:`sec:testing:ci`. In addition, the Python package used to generate the mechanism files, CEPTR, has its own set of tests (:ref:`sec:testing:ceptr`).

.. _sec:testing:building:

Building and running the tests
==============================

With the exception of ``Radiation``, each directory contains a ``GNUmakefile`` and is built with GNUmake, following the procedure described for `ReactEval` in :ref:`sec:BuildingRunning`: ::

  cd ${PELE_PHYSICS_HOME}/Testing/Exec/EosEval
  make TPL
  make -j 4
  ./Pele2d.gnu.ex inputs.2d

The name of the executable depends on the build options (dimension, compiler, and so on). Each ``GNUmakefile`` sets the usual AMReX options (``DIM``, ``COMP``, ``DEBUG``, ``USE_MPI``, ``USE_OMP``, ``USE_CUDA``, ``USE_HIP``, ...), followed by the PelePhysics model choices, and then includes ``Testing/Exec/Make.PelePhysics``, which contains the build rules shared by all of the tests. The main PelePhysics options are:

* ``Eos_Model``: the equation of state, one of ``GammaLaw``, ``Fuego``, ``Soave-Redlich-Kwong``, or ``Manifold``
* ``Chemistry_Model``: the name of a directory in ``${PELE_PHYSICS_HOME}/Mechanisms``. ``Null`` must be used with the ``GammaLaw`` EOS.
* ``Transport_Model``: the transport model, one of ``Constant``, ``Sutherland``, ``Simple``, or ``Manifold``
* ``Manifold_Dim``: the number of manifold dimensions, when the ``Manifold`` EOS is used
* ``PELE_USE_KLU`` and ``PELE_USE_MAGMA``: enable the optional SuiteSparse and MAGMA dependencies

Any of these options may be changed by editing the ``GNUmakefile`` or by specifying them on the command line, which is how a single test directory is used to cover several models: ::

  make -j 4 Eos_Model=Soave-Redlich-Kwong Chemistry_Model=LiDryer Transport_Model=Simple

.. note:: The ``make TPL`` step and the build must be done with the same set of compile-time options. Use ``make realclean`` before rebuilding a test with different models.

All of the programs accept AMReX-style inputs, either from the input file given as the first argument or appended on the command line (``./Pele2d.gnu.ex inputs.2d gridsize=64``).

EosEval
=======

This test evaluates EOS functions on a single-level, periodic, uniform grid. The temperature, density, and mass fractions are initialized with a smooth variation across the domain. The program then computes the specific heat at constant pressure, the species production rates (set to zero for the ``GammaLaw`` EOS), and the temperature recovered from the internal energy and the mass fractions. The test passes if the program runs to completion.

The default build is 2D and uses the ``Fuego`` EOS with the ``LiDryer`` mechanism. Setting ``FUEGO_GAS = FALSE`` in the ``GNUmakefile`` selects the ``GammaLaw`` EOS with ``Null`` chemistry, and setting ``USE_MANIFOLD = TRUE`` selects the ``Manifold`` EOS and transport.

The relevant input options are: ::

  gridsize      = 128      # Number of cells in each direction
  max_size      = 32       # Maximum box size
  amr.do_plot   = 1        # Write a plotfile with the results
  amr.plot_file = plt      # Plotfile prefix (written as plt00001)

The plotfile contains the density, mass fractions, specific heat, production rates, temperature, and internal energy. Three input files are provided:

* ``inputs.2d``: an empty file, so the defaults listed above are used. This file is used for the ``Fuego``, ``GammaLaw``, and ``Soave-Redlich-Kwong`` EOS.
* ``inputs.2d.table2d``: for the ``Manifold`` EOS, using the tabulated manifold ``Support/CMLM/peletable_2dim_10grid``
* ``inputs.2d.network``: for the ``Manifold`` EOS, using the neural network manifold ``Support/CMLM/fgm_net.pnn``

TranEval
========

This test evaluates the transport coefficients on a single-level, periodic, uniform grid: the mixture-averaged species diffusion coefficients, the dynamic viscosity, the bulk viscosity, and the thermal conductivity. The results are written to a plotfile with variables ``D0``, ``D1``, ..., ``mu``, ``xi``, and ``lam``. Details on the transport models can be found in :ref:`sec:transport`.

The default build is 2D and uses the ``Fuego`` EOS with the ``air`` mechanism and ``Simple`` transport. The ``GNUmakefile`` also sets ``USE_LCT_TRANSPORT = TRUE``, which compiles the multicomponent thermal conductivity models of ``Simple`` transport. One input file is provided for each transport model:

* ``inputs.2d_Constant``, ``inputs.2d_Sutherland``: values of the model parameters for the ``Constant`` and ``Sutherland`` models
* ``inputs.2d_Simple``: an empty file, as ``Simple`` transport requires no inputs
* ``inputs.2d_Manifold``, ``inputs.2d_Network``: for a build with the ``Manifold`` EOS and transport, with transport properties obtained from a separate table or from the same neural network as the EOS, respectively
* ``inputs.2d_Simple_Matrix``: for ``Simple`` transport, additionally computes the multicomponent diffusion matrix (see :ref:`sec:transport:multicomponent`) and verifies it
* ``inputs.2d_Simple_LCT``: for ``Simple`` transport, additionally computes the multicomponent thermal conductivity and the thermal diffusion ratios (see :ref:`sec:transport:conductivity`) and verifies them

The options specific to this test are: ::

  tran_eval.npts               = 128    # Number of cells in each direction
  tran_eval.matrix_diffusion   = 0      # Compute and check the multicomponent diffusion matrix
  tran_eval.write_matrix       = 1      # Add the diffusion matrix to the plotfile
  tran_eval.conductivity_check = 0      # Compute and check the multicomponent conductivity
  max_size                     = 32     # Maximum box size
  amr.plot_file                = plt    # Plotfile prefix (written as plt00001)

When ``tran_eval.matrix_diffusion = 1`` or ``tran_eval.conductivity_check = 1``, the result of each individual check is printed along with ``PASS`` or ``FAIL``, and the program aborts if any check exceeds its tolerance. The tolerances and the number of cells that are compared against a reference dense solve can be adjusted with the ``tran_eval.matrix_*`` and ``tran_eval.lct_*`` options; see ``main.cpp`` for the complete list and the default values. The conductivity checks are only performed if the code is compiled with ``USE_LCT_TRANSPORT = TRUE`` and one of the multicomponent models is selected with ``transport.conductivity_model``. Without these options, the test passes if the program runs to completion.

ReactEval
=========

This test integrates the chemical reactions in every cell of a single-level, periodic grid, using any of the reactors available in PelePhysics. A fuel/oxygen/nitrogen mixture at atmospheric pressure is initialized with a temperature that varies between 2000 K and 2500 K along the second direction, so that the cells react at different rates. The default build is 2D and uses the ``Fuego`` EOS with the ``drm19`` mechanism.

A detailed walk-through of this test with CVODE, including the meaning of the ``ode`` and ``cvode`` options and the expected output, is given in :ref:`sec:subsReactEvalCvode`, and the options of the CVODE reactor are described in :ref:`sec:subsPPOptions`. The inputs read by the test program itself are: ::

  fuel_name        = CH4           # Fuel species [Required]: H2, CH4, NC12H26, or IC8H18
  chem_integrator  = ReactorCvode  # Reactor [Required]: ReactorCvode, ReactorArkode, ReactorRK64, ReactorBDF, or ReactorNull
  ncells           = 8 64 8        # Number of cells in each direction
  max_grid_size    = 16            # Maximum box size
  plotfile         = plt           # Plotfile prefix. If omitted, no plotfiles are written
  reactFormat      = Array4        # Reactor interface that is called: Array4 or 1dArray
  initFromFile     = 0             # Initialize the data from a plotfile instead
  initFile         = <plotfile>    # Plotfile used when initFromFile = 1

  ode.dt           = 1.e-5         # Time over which the reactions are integrated
  ode.ndt          = 1             # Number of reactor calls used to cover ode.dt
  ode.reactor_type = 1             # 1: internal energy (as in PeleC), 2: enthalpy (as in PeleLMeX)
  ode.ode_ncells   = 1             # Number of cells integrated per call, with reactFormat = 1dArray
  ode.rtol         = 1.e-10        # Relative tolerance
  ode.atol         = 1.e-10        # Absolute tolerance
  ode.use_typ_vals = 0             # Use typical values to set the absolute tolerances

The fuel species must be present in the mechanism the test is compiled with. When ``plotfile`` is specified, the initial state is written to ``plt00000`` and the state after integration to ``plt00001``, with the species stored as ``SPEC0``, ``SPEC1``, ... and the temperature as ``TEMP``. The test passes if the program runs to completion.

The following input files are provided:

* ``inputs.3d_Array4``: uses the ``Array4`` interface of the reactors, in which a box of cells is passed to the reactor
* ``inputs.3d_1dArray``: uses the ``1dArray`` interface, in which the data is flattened and passed to the reactor ``ode.ode_ncells`` cells at a time
* ``inputs.3d-regt_GPU``: a short integration on a larger number of cells
* ``inputs.3d_ChkFile``: initializes the state and the external forcing terms from a plotfile. This requires ``ode.reactor_type = 2`` and a plotfile containing the ``rhoY(<spec>)``, ``F_rhoY(<spec>)``, ``rhoH``, ``Temp``, and ``F_rhoH`` variables, which is not provided in the repository.

IgnitionDelay
=============

This test computes the ignition delay time of a 0D reactor, defined as the time at which the temperature exceeds its initial value by 400 K. The calculation is done in two passes, driven by the ``exec_ignDelay.sh`` script: a first run provides an estimate of the ignition delay, which is used to generate a second input file (``inputs/inputs.0d_refine``) with an integration time and a time increment adapted to that estimate.

The default build uses the ``Fuego`` EOS with the ``dodecane_lu`` mechanism. Note that the calculation of the initial composition from the equivalence ratio in ``GPU_misc.H`` assumes that the fuel is n-dodecane, and needs to be modified for another fuel. The test is run as follows (``numpy`` is required by the Python scripts): ::

  cd ${PELE_PHYSICS_HOME}/Testing/Exec/IgnitionDelay
  make TPL
  make -j 4
  bash exec_ignDelay.sh firstpass
  python check_ignDelay.py

The argument of ``exec_ignDelay.sh`` selects the input file ``inputs/inputs.0d_<argument>``. The script uses the first executable named ``Pele*.ex`` found in the directory, so only one executable should be present. The temperature history is written to ``PPreaction.txt`` and the resulting ignition delay, along with its uncertainty (the time increment of the second pass), is written to the file ``log``.

In addition to the inputs of ``ReactEval``, the initial conditions are specified with: ::

  hr.t0             = 600         # Initial temperature [K]
  hr.press          = 1013250.0   # Pressure [dyn/cm^2]
  hr.equiv_ratio    = 1.0         # Equivalence ratio
  hr.outputFolderHR = ./          # Directory where PPreaction.txt is written

The ``check_ignDelay.py`` script provides the pass/fail signal: it exits with an error if the computed ignition delay differs from a reference value by more than a multiple of the uncertainty. The reference value and the multiple are optional arguments, and default to ``0.0763`` and ``10``, which correspond to ``inputs/inputs.0d_firstpass``.

A second input file, ``inputs/inputs.0d_manifold``, repeats the calculation with a tabulated manifold model (``Support/CMLM/autoignition.ctb``). In this case the code must be compiled with ``Eos_Model=Manifold Transport_Model=Manifold Manifold_Dim=1 Chemistry_Model=Null``, and the test is run with ``bash exec_ignDelay.sh manifold``.

Jacobian
========

This test verifies the analytical chemical Jacobian of a mechanism. Species concentrations are set to (reproducible) random values at a temperature of 1000 K, and the derivatives of the species production rates with respect to the species concentrations obtained from the analytical Jacobian (``aJacobian``) are compared with a centered finite difference of the production rates (``productionRate``). The program takes no input file and prints no output: it returns a non-zero exit code if the difference between the two Jacobians is too large.

The default build is 3D and uses the ``dodecane_lu_qss`` mechanism, as the main purpose of this test is to verify the symbolic Jacobian generated by CEPTR for mechanisms with QSS species (see :ref:`the QSS documentation <sec_qss>`). In the CI, the mechanism is regenerated with the symbolic Jacobian before the test is built: ::

  cd ${PELE_PHYSICS_HOME}/Support/ceptr
  MECH=${PELE_PHYSICS_HOME}/Mechanisms/dodecane_lu_qss
  poetry run qssa -f ${MECH}/skeletal.yaml -n ${MECH}/non_qssa_list.yaml
  poetry run convert -f ${MECH}/qssa.yaml --qss_format_input ${MECH}/qssa_input_non_optimized.toml --qss_symbolic_jacobian

  cd ${PELE_PHYSICS_HOME}/Testing/Exec/Jacobian
  make TPL DEBUG=TRUE
  make -j 4 DEBUG=TRUE
  ./Pele3d.gnu.DEBUG.ex

PMFReader
=========

This test exercises the premixed flame initialization utility described in :ref:`sec:Utility`. A 1D flame profile is read from the provided data file, ``drm19_pmf.dat``, and used to fill the velocity, density, species, enthalpy, and temperature on a 3D grid, with the flame normal to the last direction. The result is written to the plotfile ``pltInitData``. The test passes if the program runs to completion.

The default build is 3D and uses the ``drm19`` mechanism, which must match the species contained in the data file. The provided input file, ``inputs.3d``, contains: ::

  geometry.is_periodic = 0 0 0              # For each dir, 0: non-perio, 1: periodic
  geometry.coord_sys   = 0                  # 0 => cart, 1 => RZ
  geometry.prob_lo     = 0.0    0.00  0.0   # x_lo y_lo (z_lo)
  geometry.prob_hi     = 0.016  0.016 0.032 # x_hi y_hi (z_hi)

  amr.n_cell         = 32 32 64             # Number of cells in each direction
  amr.max_grid_size  = 32                   # Maximum box size

  pmf.datafile       = drm19_pmf.dat        # PMF data file
  pmf.do_cellAverage = 1                    # Average the profile over each cell
  pmf.v              = 1                    # Verbosity

  pele.standoff      = -0.022               # Shift applied to the flame position

Radiation
=========

This directory contains tests of the P1 radiation solvers located in ``Source/Radiation``. Unlike the other tests, there is no ``GNUmakefile``. The tests are defined in a ``CMakeLists.txt`` file that is meant to be added to the CMake build of a code that includes PelePhysics as a submodule, as is done in PeleLMeX. Each test is then registered with CTest, and run with the input file of the same name in the ``inputs`` directory.

The tests that are built depend on the options of the parent project:

* ``POneSingle``: single-level solve on a non-periodic domain, compared to an exact solution
* ``POneMulti``: multi-level solve on a non-periodic domain, compared to an exact solution (3D only)
* ``POneSingleAF``, ``POneMultiAF``: single-level and multi-level solves in which the absorption coefficient is computed from fields of temperature, pressure, species (CO2, H2O, and CO), and soot volume fraction, using the Planck mean absorption coefficient data files in the ``kpDB`` directory
* ``POneSingleEB``, ``POneMultiEB``: single-level and multi-level solves on a periodic domain with an embedded boundary, compared to an exact solution. These replace the tests above when EB is enabled.

The multi-level tests ``POneMulti`` and ``POneMultiAF`` are each run a second time (``POneMulti_LevbyLev`` and ``POneMultiAF_LevbyLev``) with ``composite_solve = 0``, in which case the levels are solved one at a time rather than with a composite solve.

The tests that have an exact solution print the normalized L1 norm of the error and return a non-zero exit code if it exceeds a tolerance. The ``AF`` tests pass if the program runs to completion.

The input files contain the options for the grid and for the AMReX MLMG linear solver, for example: ::

  max_level     = 2                     # Maximum level of refinement
  ref_ratio     = 2                     # Refinement ratio
  n_cell        = 32                    # Number of cells in each direction on the base level
  max_grid_size = 32                    # Maximum box size

  composite_solve = 1                   # Composite solve (1) or level-by-level solve (0)
  max_iter        = 200                 # Maximum number of MLMG iterations
  reltol          = 1.0e-5              # Relative tolerance
  abstol          = 1.0e-5              # Absolute tolerance
  linop_maxorder  = 2                   # Maximum order of the linear operator

  lo_bc = Robin Robin Neumann           # Boundary conditions on the low side of the domain
  hi_bc = Robin Robin Neumann           # Boundary conditions on the high side of the domain

CeptrTest
=========

This test compares a constant-volume 0D reactor computed with PelePhysics against the same reactor computed with `Cantera <https://cantera.org>`_, in order to verify the mechanism files generated by CEPTR. It consists of three parts:

* ``PelePhysicsSimulation``: a program similar to ``ReactEval`` that integrates a single cell and, after each reactor call, writes the time, temperature, density, and viscosity to ``results.txt``. The default build is 3D and uses the ``Fuego`` EOS with the ``dodecane_lu`` mechanism and ``Simple`` transport. The provided input file is ``inputs-cvode-0d.inp``.
* ``CanteraSimulation/homogeneous-reactor.py``: a Python script that runs the equivalent reactor in Cantera and also writes ``results.txt``. The mechanism is specified with the ``-f`` option, e.g., ``-f ${PELE_PHYSICS_HOME}/Mechanisms/dodecane_lu/mechanism.yaml``.
* ``compare.py``: a Python script that takes the Cantera (``-c``) and PelePhysics (``-p``) results files, converts the PelePhysics results from CGS to MKS units, plots both in ``plots.pdf``, and raises an error if the temperature or the density differ by more than a relative tolerance.

The Python scripts require ``cantera``, ``pandas``, and ``matplotlib``. The plots are rendered with LaTeX, which must also be available.

.. note:: The initial conditions (temperature, pressure, and composition) and the integration time are specified separately in the PelePhysics program and in the Cantera script, and must be kept consistent if either is modified.

UtilityUnitTest
===============

This program contains unit tests of the functions in ``Source/Utility``: interpolation in a tabulated function, the unit conversions between MKS and CGS, the ``locate`` function, and the rectangle-circle intersection area. The default build is 2D and uses the ``GammaLaw`` EOS with ``Null`` chemistry. It is run with the provided file named ``input``, which specifies the location of the table used in the first test: ::

  ./Pele2d.gnu.ex input

Each group of tests is checked with assertions, and a message such as ``Unit Conversion Tests Passed`` is printed when the group succeeds. A failed assertion aborts the program.

Regression tests
================

The ``Testing/Regression`` directory contains a configuration file, ``PelePhysics-tests.ini``, for the `regression testing tools <https://github.com/AMReX-Codes/regression_testing>`_ distributed by the AMReX developers, along with a ``README.md`` that describes how to set up and run the tests. These tools pull the latest version of AMReX and PelePhysics in a scratch area, build and run each test defined in the configuration file, compare the resulting plotfile against a previously stored benchmark, and assemble the results in a web page. In summary, once the repositories are cloned and the paths in the configuration file are set, the benchmarks are generated with: ::

  ./regtest.py --make_benchmarks "<a useful comment>" PelePhysics-tests.ini

and the tests are subsequently run with: ::

  ./regtest.py PelePhysics-tests.ini

.. note:: The test definitions in ``PelePhysics-tests.ini`` refer to directories and input files that are no longer present in the repository (e.g., ``Testing/Exec/ReactEval_C``), and need to be updated before this file can be used. The tests that are run automatically are those of the CI described below.

.. _sec:testing:ci:

Continuous integration
======================

Every pull request to the ``development`` branch triggers a set of GitHub Actions workflows, defined in ``.github/workflows``. The main workflow (``ci.yml``) performs the following:

* Formatting: verifies that the C++ files in ``Source``, ``Testing``, and ``Mechanisms`` conform to ``clang-format``
* Codespell: checks the repository for common misspellings
* CEPTR: checks the formatting of the Python code, runs the CEPTR tests (:ref:`sec:testing:ceptr`), and converts the ``air`` and ``LiDryer`` mechanisms
* Transport-EOS-Reactions: builds the programs of ``Testing/Exec`` with GNUmake, in 2D, for several compilers and GPU backends (GNU, LLVM, CUDA, HIP, and SYCL). The executables are run for the GNU and LLVM builds only; the GPU builds are only compiled.

The last job covers the following:

* ``UtilityUnitTest``
* ``TranEval`` for the ``Constant``, ``Simple``, ``Sutherland``, and ``Manifold`` transport models, including the multicomponent diffusion and conductivity checks for ``Simple`` transport
* ``EosEval`` for the ``Fuego``, ``GammaLaw``, ``Soave-Redlich-Kwong``, and ``Manifold`` EOS
* ``ReactEval`` for the ``drm19`` and ``LiDryer`` mechanisms, with the ``ReactorCvode``, ``ReactorArkode``, and ``ReactorRK64`` reactors and both the ``Array4`` and ``1dArray`` interfaces
* ``IgnitionDelay`` with the ``dodecane_lu`` mechanism and with the tabulated manifold
* ``Jacobian`` with the ``dodecane_lu_qss`` mechanism

The ``Manifold`` cases are not built with SYCL, and ``IgnitionDelay`` and ``Jacobian`` are not built for the GPU backends. The ``PMFReader``, ``CeptrTest``, and ``Radiation`` tests are not part of this workflow.

A second workflow (``downstream.yml``) verifies that changes to PelePhysics do not break the codes that depend on it: PeleC and PeleLMeX are built against the modified version of PelePhysics, a short case of each is run, and the PelePhysics sources are checked with ``clang-tidy`` through the CMake builds of these codes. Finally, ``docs.yml`` builds this documentation, with warnings treated as errors.

.. _sec:testing:ceptr:

CEPTR tests
===========

The tests of CEPTR are located in ``Support/ceptr/tests`` and are run with ``pytest``. After setting up ``poetry`` following the :ref:`CEPTR instructions <sec_ceptr_software>`: ::

  cd ${PELE_PHYSICS_HOME}/Support/ceptr
  poetry run pytest

These tests generate the mechanism files for a few of the mechanisms included in PelePhysics (``air``, ``LiDryer``, and ``dodecane_lu``) and verify some of the functions used to write the equilibrium constants.
