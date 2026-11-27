.. highlight:: rst

.. _sec:reactions:

*********
Reactions
*********

The ``Source/Reactions`` directory of `PelePhysics` contains the chemistry integrators, referred to as `reactors`, that advance the
species and energy of a set of cells over a time interval, subject to chemical reactions and to constant external source terms provided
by the flow solver. The following reactors are available:

* ``ReactorCvode``, an implicit integrator based on the CVODE package of `SUNDIALS`, with a selection of linear solvers
* ``ReactorArkode``, an explicit Runge-Kutta integrator with adaptive time stepping based on the ARKODE package of `SUNDIALS`
* ``ReactorRK64``, a self-contained explicit six-stage Runge-Kutta integrator with adaptive time stepping
* ``ReactorBDF``, a self-contained implicit integrator with a fixed number of substeps
* ``ReactorNull``, which only applies the external source terms and does not integrate any chemistry

The choice between these reactors is made at run time. The chemical mechanism, in contrast, is selected at compile time
through the ``Chemistry_Model`` parameter in the ``GNUmakefile`` (see :ref:`sec:ceptr`).

This section describes the interface shared by all reactors and the input file options of each of them. A more detailed discussion of the
CVODE-based reactor, with examples, is provided in the following sections of this chapter.

Selecting a reactor
===================

All reactors derive from the ``pele::physics::reactions::ReactorBase`` class, which relies on the ``Factory`` class of
``Source/Factory.H`` to register each reactor under a keyword. An application creates a reactor from this keyword: ::

    std::string chem_integrator = "ReactorCvode";
    auto reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(reactor_type, ncells);

The accepted keywords are ``ReactorCvode``, ``ReactorArkode``, ``ReactorRK64``, ``ReactorBDF`` and ``ReactorNull``. If any other keyword is provided,
the list of registered reactors is printed and the code aborts.

`PelePhysics` does not read the keyword itself: the name of the corresponding input file option is defined by the application.
The test cases distributed with `PelePhysics` (``Testing/Exec/ReactEval`` and ``Testing/Exec/IgnitionDelay``) use ::

    chem_integrator = "ReactorCvode"

while `PeleC` and `PeleLMeX` use ``pelec.chem_integrator`` and ``peleLM.chem_integrator``, respectively.

.. note:: When the ``Manifold`` equation of state is used, only ``ReactorRK64`` and ``ReactorNull`` are supported. The code aborts if the equation of state parameters are passed to any other reactor.

The reactor interface
=====================

The ``init`` function takes the reactor type as its first argument, which defines the form of the energy equation:

* ``1``: the reactor is formulated in terms of internal energy, as needed by `PeleC`
* ``2``: the reactor is formulated in terms of enthalpy, as needed by `PeleLMeX`

Any other value leads to an abort. The test cases distributed with `PelePhysics` read this value from the ``ode.reactor_type`` input file option.
In both cases, the state integrated in each cell is made up of the species partial densities :math:`\rho Y_k` and the temperature :math:`T`,
i.e. ``NUM_SPECIES + 1`` unknowns, and the temperature equation is obtained from the conservation of energy.
The corresponding equations are described in :ref:`sec:subsDiffReacts`.

Two versions of the ``react`` function are provided by each reactor:

* a version operating on the ``Array4`` data of a box, which takes the species partial densities, the species external source terms, the temperature, the energy
  (:math:`\rho E` or :math:`\rho H` depending on the reactor type), the energy external source term, an array used to return an estimate of the integration cost in each cell, and an integer mask;
* a version operating on one-dimensional arrays, where the data of ``ncells`` cells has already been gathered by the caller.

Both versions take the time interval ``dt_react`` to integrate over and the current ``time``. The external source terms are held constant during the integration.
The ``flatten`` and ``unflatten`` functions are provided to convert between the two data layouts; the temperature is recomputed from the energy and composition in both of these functions.

The cost estimate returned by the ``Array4`` version is the number of right-hand side evaluations for ``ReactorCvode`` and ``ReactorArkode``, the number of substeps taken for
``ReactorRK64`` and ``ReactorBDF``, and zero for ``ReactorNull``.

.. note:: The mask is only used by ``ReactorCvode`` in CPU builds, where cells with a mask value of ``-1`` are not integrated. It is ignored by the other reactors, and by ``ReactorCvode`` in GPU builds.

In CPU builds, ``ReactorCvode`` integrates the cells of a box one at a time, each cell being an independent system of ``NUM_SPECIES + 1`` equations.
In GPU builds, and with ``ReactorArkode`` in all cases, all the cells of the box are gathered into a single system that is advanced with a common time step.
``ReactorRK64`` and ``ReactorBDF`` always integrate each cell independently.

Options shared by the reactors
==============================

The run time options of the reactors are read with the ``ode`` prefix. The following options are shared by several reactors: ::

    ode.verbose             = 0       # Verbosity level (all reactors)
    ode.rtol                = 1e-6    # Relative tolerance (ReactorCvode, ReactorArkode)
    ode.atol                = 1e-10   # Absolute tolerance (ReactorCvode, ReactorArkode, ReactorRK64)
    ode.clean_init_massfrac = 0       # Clip the initial species in flatten/unflatten (all reactors but ReactorNull)
    ode.atomic_reductions   = 1       # Use atomic reductions in the SUNDIALS GPU vectors (ReactorCvode, ReactorArkode), 0 by default with HIP

When ``ode.clean_init_massfrac`` is activated, the species partial densities passed to the integrator by the ``flatten`` function are clipped between zero and the density,
and only the increment computed by the integrator is added back to the unclipped state in the ``unflatten`` function.

With ``ode.verbose`` larger than zero, the reactors print their options during initialization; with a value larger than one, ``ReactorCvode`` and ``ReactorArkode`` also print the solver statistics after each integration.

For ``ReactorCvode`` and ``ReactorArkode``, the absolute tolerance is a vector with one entry per unknown. By default, all entries are equal to ``ode.atol``.
An application can also provide `typical values` of the species partial densities and of the temperature through the ``set_typ_vals_ode`` function of the reactor,
which takes a vector of size ``NUM_SPECIES + 1`` with the temperature in the last entry. In this case, the absolute tolerance of each unknown is its typical value multiplied by ``ode.atol``,
which then acts as a relative factor. The typical values are not used by the other reactors.

ReactorCvode
============

This reactor uses the variable-order, variable-step Backward Differentiation Formula (BDF) method of CVODE, and is appropriate for stiff chemistry.
An introduction to the method is given in the section `A brief introduction to CVODE` of this chapter. In addition to the shared options listed above, the following options are available: ::

    ode.max_nls_iters   = 3       # Maximum number of nonlinear solver iterations per step
    ode.max_fp_accel    = 2       # Number of acceleration vectors of the fixed-point nonlinear solver
    ode.print_profiling = 0       # Print the SUNDIALS profiler output (requires SUNDIALS built with profiling)

    cvode.solve_type    = [Must be User Defined]  # Linear solver, see below
    cvode.precond_type  = [No default]            # Preconditioner, required on CPU with cvode.solve_type = precGMRES
    cvode.max_order     = 2       # Maximum order of the BDF method
    cvode.max_substeps  = 10000   # Maximum number of internal steps per call (only applied in GPU builds)

The ``cvode.solve_type`` option controls the solution of the nonlinear system arising at each implicit step. The available choices
are summarized in Table :numref:`tab:CvodeSolveTypes`, where AJ refers to the analytical Jacobian provided by the mechanism files and FD to the finite difference Jacobian approximation of CVODE.

.. _tab:CvodeSolveTypes:

.. table:: Options for ``cvode.solve_type``
    :align: center

    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``solve_type``     | Method                                   | CPU                       | GPU                            |
    +====================+==========================================+===========================+================================+
    | ``dense_direct``   | Dense direct solver, FD Jacobian         | Yes                       | No                             |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``denseAJ_direct`` | Dense direct solver, AJ                  | Yes                       | No                             |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``sparse_direct``  | Sparse direct solver, AJ                 | KLU,                      | cuSolver batched QR,           |
    |                    |                                          | ``PELE_USE_KLU = TRUE``   | CUDA only                      |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``custom_direct``  | Mechanism-specific Gauss-Jordan direct   | Yes                       | No, see below                  |
    |                    | solver, AJ                               |                           |                                |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``magma_direct``   | MAGMA batched dense direct solver, AJ    | No                        | CUDA or HIP,                   |
    |                    |                                          |                           | ``PELE_USE_MAGMA = TRUE``      |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``GMRES``          | Matrix-free GMRES, no preconditioner     | Yes                       | Yes                            |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``precGMRES``      | Matrix-free GMRES, left preconditioner   | Yes, see                  | CUDA only, see                 |
    |                    |                                          | ``cvode.precond_type``    | ``cvode.precond_type``         |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``fixed_point``    | Fixed-point nonlinear solver, no linear  | Yes                       | Yes                            |
    |                    | solver                                   |                           |                                |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+
    | ``diagnostic``     | Print the sparsity pattern of the        | Yes                       | No                             |
    |                    | Jacobians and abort                      |                           |                                |
    +--------------------+------------------------------------------+---------------------------+--------------------------------+

When ``cvode.solve_type = precGMRES``, the preconditioner is built from a simplified analytical Jacobian provided by the mechanism files, and its storage and solution method are selected with ``cvode.precond_type``:

* In CPU builds: ``dense_simplified_AJacobian``, ``sparse_simplified_AJacobian`` (requires ``PELE_USE_KLU = TRUE``) or ``custom_simplified_AJacobian``. The code aborts if none of these is provided.
* In GPU builds: ``cuSparse_simplified_AJacobian``, which is only available with CUDA.

There is no default value for ``cvode.solve_type``: the code aborts with the list of valid options if it is not provided.

.. note:: ``cvode.solve_type = custom_direct`` is accepted in CUDA builds, but the corresponding linear solve is currently disabled in ``ReactorCvodeCustomLinSolver.cpp``. This option should thus only be used in CPU builds.

.. note:: In GPU builds, the unknowns of the cells gathered in one CVODE instance can be ordered in two different ways. By default, the unknowns are grouped by species, which is faster but only compatible with ``cvode.solve_type = GMRES``. All the other GPU options of Table :numref:`tab:CvodeSolveTypes` require the unknowns to be grouped by cell, which is activated by adding ``PELE_CVODE_FORCE_YCORDER = TRUE`` to the ``GNUmakefile``, and is automatically the case when ``PELE_USE_MAGMA = TRUE``. CPU builds always group the unknowns by cell.

.. note:: With HIP, the analytical Jacobian routines of the mechanism files are only compiled if ``PELE_COMPILE_AJACOBIAN = TRUE`` is added to the ``GNUmakefile``.

The optional libraries are activated with the following ``GNUmakefile`` options, and are downloaded and built together with `SUNDIALS` when executing ``make TPL`` (see :ref:`sec:BuildingRunning`): ::

    PELE_USE_KLU   = TRUE    # SuiteSparse KLU, sparse solver and preconditioner in CPU builds
    PELE_USE_MAGMA = TRUE    # MAGMA, dense batched solver with CUDA or HIP

ReactorArkode
=============

This reactor uses the explicit Runge-Kutta methods of ARKODE with adaptive time stepping, and is therefore intended for chemistry that is not stiff
over the time interval considered. In addition to the shared options listed above, the following options are available: ::

    ode.use_erkstep   = 1     # 1: use the ERKStep module of ARKODE, 0: use the ARKStep module with an explicit method
    ode.rk_method     = 40    # Explicit Runge-Kutta method, see below
    ode.rk_controller = 0     # Time step controller: 0 = PID, 1 = PI, 2 = I, 3 = explicit Gustafsson

The ``ode.rk_method`` option selects one of the explicit Butcher tables of ARKODE:

.. table:: Options for ``ode.rk_method``
    :align: center

    +-------------------+-----------------------------------+
    | ``ode.rk_method`` | ARKODE Butcher table              |
    +===================+===================================+
    | 20                | ``ARKODE_HEUN_EULER_2_1_2``       |
    +-------------------+-----------------------------------+
    | 30                | ``ARKODE_BOGACKI_SHAMPINE_4_2_3`` |
    +-------------------+-----------------------------------+
    | 31                | ``ARKODE_ARK324L2SA_ERK_4_2_3``   |
    +-------------------+-----------------------------------+
    | 40 (default)      | ``ARKODE_ZONNEVELD_5_3_4``        |
    +-------------------+-----------------------------------+
    | 41                | ``ARKODE_ARK436L2SA_ERK_6_3_4``   |
    +-------------------+-----------------------------------+
    | 42                | ``ARKODE_SAYFY_ABURUB_6_3_4``     |
    +-------------------+-----------------------------------+
    | 43                | ``ARKODE_ARK437L2SA_ERK_7_3_4``   |
    +-------------------+-----------------------------------+
    | 50                | ``ARKODE_CASH_KARP_6_4_5``        |
    +-------------------+-----------------------------------+
    | 51                | ``ARKODE_FEHLBERG_6_4_5``         |
    +-------------------+-----------------------------------+
    | 52                | ``ARKODE_DORMAND_PRINCE_7_4_5``   |
    +-------------------+-----------------------------------+
    | 53                | ``ARKODE_ARK548L2SA_ERK_8_4_5``   |
    +-------------------+-----------------------------------+
    | 54                | ``ARKODE_ARK548L2SAb_ERK_8_4_5``  |
    +-------------------+-----------------------------------+
    | 60                | ``ARKODE_VERNER_8_5_6``           |
    +-------------------+-----------------------------------+
    | 80                | ``ARKODE_FEHLBERG_13_7_8``        |
    +-------------------+-----------------------------------+

Any other value of ``ode.rk_method`` or ``ode.rk_controller`` falls back to the default. The reader is referred to the `SUNDIALS documentation <https://sundials.readthedocs.io>`_ for a description of these methods and controllers.

ReactorRK64
===========

This reactor is an explicit, six-stage Runge-Kutta integrator with an error estimate, implemented directly in `PelePhysics`. It does not depend on `SUNDIALS` vectors or solvers,
and each cell is integrated independently with its own adaptive time step. As for ``ReactorArkode``, it is intended for chemistry that is not stiff.
In addition to ``ode.verbose``, ``ode.atol`` and ``ode.clean_init_massfrac``, the following options are available: ::

    ode.rk64_nsubsteps_guess = 10     # Initial time step is dt_react / rk64_nsubsteps_guess
    ode.rk64_nsubsteps_min   = 5      # Maximum time step is dt_react / rk64_nsubsteps_min
    ode.rk64_nsubsteps_max   = 500    # Minimum time step is dt_react / rk64_nsubsteps_max

After each step, the largest absolute value of the error estimate over all unknowns is compared to ``ode.atol``, and the time step is increased or decreased accordingly, within the bounds defined above.
The step is always accepted. The relative tolerance ``ode.rtol`` and the typical values are not used by this reactor.

ReactorBDF
==========

This reactor is an implicit integrator implemented directly in `PelePhysics`, which does not depend on `SUNDIALS`. Each cell is integrated independently with a fixed number of
substeps of equal size. At each substep, the nonlinear system is solved with Newton iterations, using the analytical Jacobian provided by the mechanism files and a GMRES linear solver.
In addition to ``ode.verbose`` and ``ode.clean_init_massfrac``, the following options are available: ::

    ode.bdf_scheme          = 3       # 0 = BDF1 (backward Euler), 1 = trapezoidal, 2 = BDF2, 3 = BDF3
    ode.bdf_nsubsteps       = 1       # Number of substeps per call
    ode.bdf_nonlinear_iters = 1       # Maximum number of Newton iterations per substep
    ode.bdf_nonlinear_tol   = 1e-6    # Tolerance on the norm of the nonlinear residual
    ode.bdf_gmres_precond   = 2       # GMRES preconditioner: 0 = none, 1 = Jacobi, 2 = symmetric Gauss-Seidel
    ode.bdf_gmres_restarts  = 10      # Number of GMRES restarts
    ode.bdf_gmres_kspiters  = NUM_SPECIES + 1   # Number of Krylov iterations per restart
    ode.bdf_gmres_tol       = 1e-10   # GMRES tolerance

The BDF2 and BDF3 schemes need the solution at previous substeps, which is not available at the beginning of each call: the trapezoidal scheme is used instead for the first substep with BDF2, and for the first two substeps with BDF3.
In particular, with the default values of ``ode.bdf_scheme`` and ``ode.bdf_nsubsteps``, a single trapezoidal step is performed.
There is no error control: the tolerances ``ode.rtol`` and ``ode.atol`` are not used, and the accuracy is controlled by the user through ``ode.bdf_nsubsteps``.

.. note:: In GPU builds, this reactor is only available for mechanisms with less than 55 species.

ReactorNull
===========

This reactor does not integrate the chemistry. The ``Array4`` version of the ``react`` function adds the external source terms multiplied by ``dt_react`` to the species partial densities and to the energy, and
recomputes the temperature; the one-dimensional version leaves the state untouched. It is intended for non-reacting calculations. The only option read by this reactor is ``ode.verbose``.
