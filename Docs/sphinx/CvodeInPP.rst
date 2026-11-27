.. highlight:: rst

.. role:: cpp(code)
   :language: c++

CVODE implementation in `PelePhysics`
=====================================

.. _sec:subsDiffReacts:

The different reactors
----------------------

Throughout this document, what we call a `reactor` is in fact a zero-dimensional model, 
the simplest representation of a chemically reacting system. Depending upon the choice of state variables 
driving the system, several different types of reactor can be considered; 
and the "correct" choice is case dependent. In general, the state variables for a reactor model are

- The reactor mass
- The reactor volume
- The energy of the system
- The mass fractions for each species

The most common type of reactor is the `constant-volume` (CV) reactor, which is the one used to advance the chemistry 
within `PeleC`. This reactor type is equivalent to a rigid vessel with fixed volume but variable pressure. 
In `PelePhysics`, the constant-volume constraint is ensured by keeping the density :math:`\rho` fixed 
-since there is no change of mass; and the indirect choice of energy in the CV reactor implementation is the total energy 
:math:`E`. :math:`E`'s evolution in our case is solely due to a constant external source term :math:`\dot{E}_{ext}`, which accounts 
for the effects of advection and convection in the Spectral Deferred Correction (SDC) scheme that all `Pele` codes use (see the `PeleLM <https://pelelm.readthedocs.io/en/latest/index.html>`_ documentation for example). 
In that sense, the CV reactor is an abstraction and is not a true closed vessel.

Note that CVODE still integrates the mass fractions (:math:`\rho Y`) together with energy for stability reasons, 
but a change of variable is applied to effectively transport the temperature :math:`T` via

.. math::

    \rho C_v \frac{\partial T}{\partial t} = \rho\dot{E}_{ext}  - \sum_k e_k {\dot{\omega}_k}^M

where the :math:`e_k` are the species internal energy and :math:`{\dot{\omega}_k}^M` is the species :math:`k` mass production rate. 

In a second implementation, that we will label `constant-volume-enthalpy` (CVH), the mass-weighted total enthalpy :math:`\rho H` is used and 
conserved along with :math:`\rho`. This reactor type is also an abstraction. Here also, :math:`\rho H` 
evolves according to an external source term :math:`\dot{\rho H}_{ext}`, and in CVODE, the mass fractions (:math:`\rho Y`) and 
temperature :math:`T` are integrated according to

.. math::

    \rho C_p \frac{\partial T}{\partial t} = \rho\dot{H}_{ext}  - \sum_k h_k  {\dot{\omega}_k}^M

where the :math:`h_k` are the species enthalpy. 

.. _sec:subsubValidCVreact:

Validation of the CV reactor implementation in CVODE (with CANTERA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`CANTERA <https://cantera.org/>`_ is an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes. 
It is a very robust and fast tool written in C++ that is also based on CVODE to perform the chemistry integration. 
CANTERA is well recognized in the combustion community, and by comparing our results to CANTERA reactor simulations, 
we will be able to validate our implementation. 

Note that only the CV reactor model described above can be validated, since as we discussed before, 
the CVH reactor model is an abstraction needed for our Low-Mach PeleLM chemistry integration. Also, to have a real CV reactor, 
the external source terms for the energy and species equations in `PelePhysics` have been set to 0 (see :ref:`sec:subsDiffReacts`).

The parameters chosen to initialize the simulation in both CANTERA and `PelePhysics` are described in 
Table :numref:`fig:ReactEvalCVODETable`. The kinetic mechanism used for hydrogen combustion is available in `PelePhysics`. 
Note that small sub-steps are explicitly taken until the final time is reached, 
but CVODE's internal machinery can subdivides the :math:`dt` even further. 
For the purpose of validation, the direct dense solver of CVODE is selected 
in `PelePhysics` (see section :ref:`sec:subsPPOptions`).

.. _fig:ReactEvalCVODETable:

.. table:: Parameters for initializing simulation
    :align: center

    +------------+-----------------+-------------+----------------+-------------+----------------+-----------------+
    | Mechanism  |     Mixture     |  Initial T  |  Initial phi   |   Pressure  |       dt       |    Final time   |
    +------------+-----------------+-------------+----------------+-------------+----------------+-----------------+
    |  Li Dryer  |      H2/O2      |   1500 K    |      0.8       |  101325 Pa  |     1.0e-8s    |     3.0e-6s     |
    +------------+-----------------+-------------+----------------+-------------+----------------+-----------------+

Results are plotted in Fig :numref:`fig:ReactEvalCVODE` and :numref:`fig:ReactEvalCVODESpecs`. for the :math:`H_2/O_2` mixture. 
All curves are indistinguishable, so the relative error of all major quantities is also plotted in Fig. :numref:`fig:ReactEvalCVODEErrss`. 
Note that :math:`H_2` and :math:`O_2` relative errors have similar features, and that relative errors observed 
for :math:`H` and :math:`H_2O` are representative of those exhibited by, respectively, intermediates and products.

.. |a| image:: ./Visualization/Main.001.png
     :width: 100%

.. |b| image:: ./Visualization/Specs.png
     :width: 100%

.. |c| image:: ./Visualization/ERRs.png
     :width: 100%

.. _fig:ReactEvalCVODE:

.. table:: Evolution of temperature, pressure and enthalpy in a CV reactor, computed with the LiDryer mechanism. Black: CANTERA, red: PelePhysics.
     :align: center

     +-----+
     | |a| |
     +-----+
..
    .. figure:: ./Visualization/Main.001.png
     :width: 100%
     :name: fig-ReactEvalCVODE
     :alt: Evolution of temperature, pressure and enthalpy in a CV reactor, computed with the LiDryer mechanism. Black: CANTERA, red: PelePhysics.


.. _fig:ReactEvalCVODESpecs:

.. table:: Evolution of major species in a CV reactor, computed with the LiDryer mechanism. Black: CANTERA, red: PelePhysics. 
     :align: center

     +-----+
     | |b| |
     +-----+

.. _fig:ReactEvalCVODEErrss:

.. table:: Relative errors on the temperature, pressure, enthalpy and major species in a CV reactor, computed with the LiDryer mechanism. 
     :align: center

     +-----+
     | |c| |
     +-----+


Overall, considering the many CVODE controlling parameters, results are deemed acceptable and that 
concludes the validation of the reactors implemented in `PelePhysics`.



.. _sec:subsPPOptions:

Activating the different CVODE solver options via the input files
-----------------------------------------------------------------
**Note that at this point, it is believed that the user has properly installed CVODE as well as the SuiteSparse package. If not, refer to** :ref:`sec:GetCVODE`.

Choosing between CVODE and the other ODE integrators available in `PelePhysics` is done at run time, via a keyword in the input file (see :ref:`sec:reactions`).
The type of reactor and specifics of the numerical algorithm are also selected via keywords in the input file. There is a subtlety though: 
when any sparsity feature is required, the choice should also be made at compile time since external libraries will be required; 
and if the compilation is not performed properly, subsequent options via keywords in the input file will lead to an error. 
This is discussed in more depth in what follows.

.. _subsubs:GNUtype:

The GNUmakefile
^^^^^^^^^^^^^^^

`SUNDIALS`, and thus CVODE, is a required dependency of `PelePhysics`: it is always built, and no ``GNUmakefile`` flag is needed to activate it. 
The user only has to choose between a number of different methods to integrate the linear system arising during the implicit solve. 
Add the following line if sparsity features are required: ::

    PELE_USE_KLU = TRUE

Note that this is a `PelePhysics` flag, so it will automatically be recognized in the `Pele` codes. `SUNDIALS` and, if requested, `SuiteSparse` are then built by executing ``make TPL``, 
as described in :ref:`sec:BuildingRunning`.
    
All of the flags discussed in this subsection are used in ``$PELE_PHYSICS_HOME/ThirdParty/Make.ThirdParty``.


The input file
^^^^^^^^^^^^^^

The input file is made up of specific blocks containing keywords that apply to specific areas of the integrationof the problem at hand. 
The suffix associated with each block of keywords should help the user in determining which keywords 
are needed in his case, depending on the integrator selected.
If CVODE is selected, for example, keywords starting with ``cvode.*`` are relevant. 
The general ``ode.*`` keywords are shared by all ODE integrators and thus are also relevant for CVODE:

- ``ode.reactor_type`` enable to switch from a CV reactor (``=1``) to a CVH reactor (``=2``).
- ``cvode.solve_type`` controls the CVODE linear integration method: choose ``dense_direct`` to enable the dense direct linear solver, 
  ``denseAJ_direct`` for the dense direct linear solver with an Analytical Jacobian, 
  ``sparse_direct`` for the sparse direct linear solver (if the KLU library has been linked), ``GMRES`` for the non-preconditioned Krylov iterative solver 
  and ``precGMRES`` for the preconditioned Krylov iterative solver. 
- ``cvode.precond_type`` selects the preconditioner when ``cvode.solve_type = precGMRES``: 

  - With ``dense_simplified_AJacobian``, the preconditioner is based on a simplified Analytical Jacobian in a dense format. 
  
  - With ``sparse_simplified_AJacobian``, **and** if the KLU library is linked, the preconditioned solve is done in a sparse format. 

The complete list of options is given in :ref:`sec:reactions`.

.. note:: Earlier versions of `PelePhysics` relied on integer values of ``cvode.solve_type`` (``1``, ``5`` and ``99``) combined with an ``ode.analytical_jacobian`` keyword. These are no longer recognized by the code.


.. _sec:subsReactEvalCvode:

The ReactEval test case with CVODE in details
---------------------------------------------

This test case is located in ``$PELE_PHYSICS_HOME/Testing/Exec/ReactEval``. 
The domain considered here is a :math:`2x1024x2` box (the number of cells is set via the input file), where the initial temperature is different in each :math:`(i,j,k)-` cell, according to a :math:`y-` evolving sinusoidal profile, see Fig. :numref:`fig:ErrH2`:

.. math::

    T(i,j,k) =  T_l + (T_h-T_l)\frac{y(i,j,k)}{L} + dTsin\left(2\pi\frac{y(i,j,k)}{P}\right) 

The different parameters involved are summarized in Table :numref:`tab::ParamReactEvalCvode`. The initial pressure is 1 atm. The initial composition is the same in every cell, and is a mixture of 0.1 :math:`C_nH_m`, 0.2 :math:`O_2` and 0.7 :math:`N_2` in mass fractions. 

Various fuels and kinetic mechanisms can be employed. For the purpose of this tutorial, two common fuels will be considered: **methane** (n=1 and m=4) and **n-dodecane** (n=12 and m=26), modelled via the **drm** and **dodecane_wang** kinetic schemes, respectively. Both mechanisms are available in `PelePhysics`.

The following focuses on the :math:`CH_4`/:math:`O_2` example, but performances for both mechanisms and initial composition will be reported in the results section.


.. _tab::ParamReactEvalCvode:

.. table:: Parameters used to initialize T in the ReactEval test case
    :align: center

    +------------+-----------------+-------------+----------------+-------------+
    | Tl         |     Th          |  dT         |  L             |   P         |
    +------------+-----------------+-------------+----------------+-------------+
    |  2000 K    |      2500 K     |   100 K     |      1024      |  L/4        |
    +------------+-----------------+-------------+----------------+-------------+


.. _fig:ErrH2:

.. figure:: ./Visualization/Case_ReactEvalCvode.001.png
     :width: 50%
     :align: center
     :name: fig-ReactEvalCVODE
     :target: ./Visualization/Case_ReactEvalCvode.001.png
     :alt: The ReactEval test case

     The ReactEval test case

The GNUmakefile
^^^^^^^^^^^^^^^

For this example, the chemistry model should be set to ``drm19``. Starting from the ``GNUmakefile`` distributed with the test case, 
the ``DIM``, ``USE_MPI`` and ``TINY_PROFILE`` flags are modified, and the full file reads as follows:

.. code-block::

    # define the location of the PELE_PHYSICS top directory
    PELE_PHYSICS_HOME    ?= ../../..
    
    # AMReX
    DIM          = 3
    PRECISION    = DOUBLE
    PROFILE      = FALSE
    VERBOSE      = FALSE
    DEBUG        = FALSE
    TINY_PROFILE = TRUE
    
    # Compiler
    COMP       = gnu
    USE_MPI    = TRUE
    USE_OMP    = FALSE
    USE_CUDA   = FALSE
    USE_HIP    = FALSE
    USE_SYCL   = FALSE
    
    # PelePhysics TPLs
    PELE_USE_KLU = FALSE
    PELE_USE_MAGMA = FALSE
    
    # this flag activates the subcycling mode in the D/Cvode routines
    DEFINES  += -DMOD_REACTOR
    
    Eos_Model       = Fuego
    Chemistry_Model = drm19
    Transport_Model = Simple
    
    Bpack   := ./Make.package
    Blocs   := .

    include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics         

Note that the ``TINY_PROFILE`` flag has been activated to obtain statistics on the run. This is an `AMREX` option.

The input file
^^^^^^^^^^^^^^

The run parameters that can be controlled via the input file for this example are as follows: ::

    # Number of cells
    ncells = 2 1024 2
    # Max size of problem
    max_grid_size = 2
    # Choose name of output pltfile
    plotfile = plt
    # Fuel species
    fuel_name = CH4
    # ODE integrator
    chem_integrator = "ReactorCvode"

    #ODE solver options
    # REACTOR mode
    ode.dt = 1.e-05  
    ode.ndt = 10
    # Reactor formalism: 1=full e, 2=full h
    ode.reactor_type = 1
    # Tolerances for ODE solve
    ode.rtol = 1e-9
    ode.atol = 1e-9
    #CVODE SPECIFICS
    # Choose the linear solver
    cvode.solve_type = dense_direct

so in this example, a **CV reactor model is chosen** to integrate each cell, and the **dense direct solve without analytical Jacobian** is activated. 
Each cell is then integrated for a total of :math:`1.e-05` seconds, with 10 external time steps. 
This means that the actual :math:`dt` is :math:`1.e-06s`, which is more than what is typically used in the `PeleC` code, 
but consistent with what is used in `PeleLM`. Note that the fuel is explicitly specified to be methane.
By default, the number of cells integrated simultaneously by one CVODE instance is 1 [#Foot1]_, but the `AMREX` block-integration proceeds by blocks of :math:`2x2x2`.


Results
^^^^^^^

It took 52.61s to integrate the 4096 cells of this box, with 4 MPI processes and no OMP process. 
The resulting temperature evolution for all cells in the y-direction is displayed in Fig. :numref:`fig:ReacEvalCv`.


.. _fig:ReacEvalCv:

.. figure:: ./Visualization/ReactEvalCv.001.png
     :width: 100%
     :align: center
     :name: fig-ReactEvalCv
     :alt: Evolution of temperature in the 2x1024x2 example box, using a CV reactor and a dense direct solve, and computed with the DRM mechanism. Black: $t=0$, red: $t=1e-05s$

     Evolution of temperature in the 2x1024x2 example box, using a CV reactor and a dense direct solve, and computed with the DRM mechanism. Black: t=0s, red: t=1e-05s


To go further: ReactEval with CVODE and the KLU library
-------------------------------------------------------

The GNUmakefile
^^^^^^^^^^^^^^^

Only one line of the ``GNUmakefile`` needs to be modified compared to the previous example, and the dependencies need to be built again with ``make TPL``.

.. code-block::

    ...
    # PelePhysics TPLs
    PELE_USE_KLU = TRUE
    ...


The input file
^^^^^^^^^^^^^^

For the KLU library to be of use, a solver utilizing sparsity features should 
be selected. We modify the input file as follows:

.. code-block::

    ...
    #ODE solver options 
    ...
    #CVODE SPECIFICS
    # Choose the linear solver
    cvode.solve_type = precGMRES
    # Choose the preconditioner
    cvode.precond_type = sparse_simplified_AJacobian
    
So that now, a preconditioned iterative Krylov solver is selected, where the preconditioner is specified in a sparse format.

Results
^^^^^^^

This run now takes 1m34s to run. As expected from the dense Jacobian of the system obtained when using the small DRM mechanism 
(the fill in pattern is :math:`>90 \%`), using an iterative solver does not enable to reach speed-ups over the simple dense direct 
solve. **NOTE**, and this is important, that this tendency will revert when sufficiently small time steps are used. 
For example, if instead of :math:`1e-6s` we took time steps of :math:`1e-8s` (consistent with `PeleC` time steps), then using 
the iterative GMRES solver would have provided significant time savings. This is because the smaller the time step the 
closer the system matrix is from the identity matrix and the GMRES iterations become really easy to complete.

This example illustrates that choosing the "best" and "most efficient" algorithm is far from being a trivial task, 
and will depend upon many factors. Table :numref:`tab:RunsReactEvalCvode` provides a summary of the CPU run time in solving the 
ReactEval example with a subset of the various available CVODE linear solvers. As can be seen from the numbers, using an AJ is much more efficient than relying upon CVODE's built-in difference quotients. Using a sparse solver does not appear to provide additional time savings.

.. note:: The run times reported in this chapter were obtained with an earlier version of `PelePhysics`. The keywords listed in the tables have been updated to their current equivalent.

.. _tab:RunsReactEvalCvode:

.. table:: Summary of ReactEval runs with various algorithms (methane/air)
    :align: center

    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  Solver              |  Direct        |  Direct          |  Direct         |  Iter.         |  Iter.                        |
    |                      |  Dense         |  Dense AJ        |  Sparse AJ      |  not Precond.  |  Precond. (S)                 |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  KLU                 |  OFF           |  OFF             |  ON             |  OFF           |  ON                           |
    +======================+================+==================+=================+================+===============================+
    |  ode.reactor_type    |  1             |  1               |  1              |  1             |  1                            |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  cvode.solve_type    |  dense_direct  |  denseAJ_direct  |  sparse_direct  |  GMRES         |  precGMRES                    |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  cvode.precond_type  |  N/A           |  N/A             |  N/A            |  N/A           |  sparse_simplified_AJacobian  |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  Run time            |  52.61s        |  44.87s          |  48.64s         |  1m42s         |  1m34s                        |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+


The same series of tests are performed for a mixture of n-dodecane and air (see :ref:`sec:subsReactEvalCvode`), the configuration being otherwise the same as in the methane/air case. Results are summarized in Table :numref:`tab:RunsReactEvalCvodeDOD`. The overall tendencies remain similar. Note that the non-preconditioned GMRES solver becomes very inefficient for this larger system. Here also, the direct sparse solve --which relies upon the KLU library, does not seem to provide additional time savings. The fill-in pattern is :math:`70 \%`.

.. _tab:RunsReactEvalCvodeDOD:

.. table:: Summary of ReactEval runs with various algorithms (n-dodecane/air)
    :align: center

    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  Solver              |  Direct        |  Direct          |  Direct         |  Iter.         |  Iter.                        |
    |                      |  Dense         |  Dense AJ        |  Sparse AJ      |  not Precond.  |  Precond. (S)                 |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  KLU                 |  OFF           |  OFF             |  ON             |  OFF           |  ON                           |
    +======================+================+==================+=================+================+===============================+
    |  ode.reactor_type    |  1             |  1               |  1              |  1             |  1                            |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  cvode.solve_type    |  dense_direct  |  denseAJ_direct  |  sparse_direct  |  GMRES         |  precGMRES                    |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  cvode.precond_type  |  N/A           |  N/A             |  N/A            |  N/A           |  sparse_simplified_AJacobian  |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+
    |  Run time            |  6m25s         |  5m33s           |  6m32s          |  21m44s        |  10m14s                       |
    +----------------------+----------------+------------------+-----------------+----------------+-------------------------------+



Current Limitations
-------------------

Note that currently, all sparse operations rely on an Analytical Jacobian. This AJ is provided via the chemistry routines dumped by the CEPTR code. Those routines are generated in a pre-processing step, when the sparsity pattern of the AJ is still unknown. As such, all entries of the AJ are computed at all times, and when a sparsity solver is chosen, the AJ is in fact "sparsified" to take advantage of the sparse linear algebra. The "sparsification" process involves a series of loop in the cpp that takes a significant amount of the CPU time most of the time. However, it is always good to verify that this is the case. `AMREX`'s ``TINY_PROFILER`` features is a handy tool to do so.

.. _sec:subssubsTricks:

Tricks and hacks, stuff to know
-------------------------------

DVODE, the Fortran integrator used in earlier versions of `PelePhysics`, included a `hack` enabling the user to reuse the Jacobian instead of reevaluating it from scratch. 
This option was triggered by setting the ``extern_probin_module`` flag ``new_Jacobian_each_cell`` to ``0`` in the ``probin`` file.

A similar feature is currently not available in CVODE, although it would be possible to modify the ``CVodeReInit`` function 
to reinitialize only a subset of counters. 
The maximum number of steps between two Jacobian evaluations is set to 100 in `PelePhysics` via the CVODE function ``CVodeSetJacEvalFrequency``.

How does CVODE compare with DVODE ?
-----------------------------------

.. note:: DVODE and the Fortran routines interfacing it are no longer distributed with `PelePhysics`: the ``USE_SUNDIALS_PP`` flag and the ``probin`` file mentioned in this section no longer exist. The comparison is kept as a reference.

Depending on whether the famous Jacobian `hack` is activated or not in DVODE, 
the completion time of the run can be decreased significantly. The same test case as that described in the previous section can also be integrated with DVODE. 
For that purpose, the FORTRAN routines implementing the DVODE integration were interfaced with C++ via a FORTRAN header. The run is thus identical to ReactEval with CVODE.
Only the ``GNUmakefile`` needed to be modified:

.. code-block:: c++

    ...
    #######################
    # ODE solver OPTIONS: DVODE (default) / SUNDIALS / RK explicit
    #######################
    # Activates use of SUNDIALS: CVODE (default) / ARKODE
    USE_SUNDIALS_PP = FALSE
    ...
    
    #######################
    ...

and, as explained in section :ref:`sec:subssubsTricks`, the famous AJ `hack` could be activated via the ``probin`` file.

Two runs are performed, activating the hack or not. Times are reported in Table :numref:`tab:CVODEvsDVODE`.

.. _tab:CVODEvsDVODE:

.. table:: Summary of a CVODE vs a DVODE chemistry integration on the same test case
    :align: center

    +--------------------+----------------+----------+------------------+
    |  Solver            |  Direct        |  Direct  |  Direct          |
    |                    |  Dense         |  Dense   |  Dense + `hack`  |
    +--------------------+----------------+----------+------------------+
    |  KLU               |  OFF           |  OFF     |  OFF             |
    +--------------------+----------------+----------+------------------+
    |  Integrator        |  CVODE         |  DVODE   |  DVODE           |
    +====================+================+==========+==================+
    |  ode.reactor_type  |  1             |  1       |  1               |
    +--------------------+----------------+----------+------------------+
    |  cvode.solve_type  |  dense_direct  |  N/A     |  N/A             |
    +--------------------+----------------+----------+------------------+
    |  Run time          |  52.61s        |  53.21s  |  52.83s          |
    +--------------------+----------------+----------+------------------+


In this case, the hack does not seem to provide significant time savings. Note also that CVODE is usually slightly more efficient than DVODE, consistently with findings of other studies available in the literature -- although in this case all options give comparable results.


CVODE implementation in `PelePhysics` on GPU
============================================

Requirements and input files
----------------------------

**To use CVODE on a GPU, Sundials should be built with GPU support. This is done automatically when executing** ``make TPL`` **with the GPU options activated in the** ``GNUmakefile`` **. Refer to** :ref:`sec:BuildingRunning`.

Note that the SuiteSparse package does not support GPU architecture and is thus no longer required. Sparse linear algebra operations, when needed, are performed with the help of CUDA's `cuSolver <https://docs.nvidia.com/cuda/cusolver/index.html>`_.

The GNUmakefile
^^^^^^^^^^^^^^^

To run on GPUs, `AMREX` should be build with CUDA enabled. To do so, add this line to the ``GNUmakefile``: ::

    USE_CUDA   = TRUE

This should activate the CUDA features of CVODE in `PelePhysics` too.


The input file
^^^^^^^^^^^^^^

In the input file, the same main keywords control the algorithm (``ode.reactor_type``, ``cvode.solve_type``, ``cvode.precond_type``). However, note that there are less linear solver options available. 

- Both preconditioned or non-preconditioned GMRES options are available (``cvode.solve_type = precGMRES`` or ``cvode.solve_type = GMRES``). The preconditioned version requires ``cvode.precond_type = cuSparse_simplified_AJacobian``.
- The user has the choice between different direct solvers. 

  - Sundials offers one option (the SUNLinSol_cuSolverSp_batchQR) relying upon the cuSolver to perform batched sparse QR factorizations. This version is enabled via ``cvode.solve_type = sparse_direct``. 
  - Another version is selected via ``cvode.solve_type = custom_direct``. This version relies upon a pre-computed Gauss-Jordan `Solver <https://github.com/accelerated-odes/gauss-jordan-solver>`_. Note that the corresponding solve is currently disabled in the source code.
  - A batched dense solver from the `MAGMA` library is enabled via ``cvode.solve_type = magma_direct``, if ``PELE_USE_MAGMA = TRUE`` is added to the ``GNUmakefile``.

All of these options but the non-preconditioned GMRES require a specific ordering of the unknowns, see :ref:`sec:reactions`.

Grouping cells together
-----------------------

To take full advantage of the GPU power, many intensive operations of similar nature should be performed in parallel. In `PelePhysics`, this is achieved by grouping many cells together, and integrating each one in separate threads within one CVODE instance. Indeed, the flow of operations to solve one set of ODEs is very similar from one cell to the next, and one could expect limited thread divergence from this approach. Fig. :numref:`fig:GroupingCells` summarizes the idea. Note that the Jacobian of the group of cells is block-sparse, and any chosen integration method should take advantage of this.


.. _fig:GroupingCells:

.. figure:: ./Visualization/GroupingOfCells.png
     :width: 40%
     :align: center
     :name: fig-GroupingCells
     :alt: Grouping of cells 

     n cells are solved together in one CVODE instance. The big-matrix is block-sparse.

In the current implementation, the number of cells that are grouped together is equal to the number of cells contained in the box under investigation within a MultiFab iteration.

The ReactEval test case on GPU in details
-----------------------------------------

A series of tests are performed on the GPU for a mixture of methane and air, with the intent of evaluationg the performance of the chemistry solvers. 
The test case, configuration and initial conditions are similar to that described in :ref:`sec:subsReactEvalCvode`. The mechanism employed is the **drm**. 

The GNUmakefile
^^^^^^^^^^^^^^^

The ``GNUmakefile`` is the same as that described in :ref:`sec:subsReactEvalCvode`, with the following modifications: ::

    USE_MPI    = FALSE
    USE_CUDA   = TRUE


The input file
^^^^^^^^^^^^^^

The Results
^^^^^^^^^^^

Results are summarized in Table :numref:`tab:RunsReactEvalCvodeDRMGPU`. 


.. _tab:RunsReactEvalCvodeDRMGPU:

.. table:: Summary of ReactEval runs on GPU with various algorithms (methane/air)
    :align: center

    +----------------------+-----------------+-----------------+----------------+---------------------------------+
    |  Solver              |  Direct         |  Direct         |  Iter.         |  Iter.                          |
    |                      |  Sparse I       |  Sparse II      |  not Precond.  |  Precond. (S)                   |
    +======================+=================+=================+================+=================================+
    |  reactor_type        |  1              |  1              |  1             |  1                              |
    +----------------------+-----------------+-----------------+----------------+---------------------------------+
    |  cvode.solve_type    |  custom_direct  |  sparse_direct  |  GMRES         |  precGMRES                      |
    +----------------------+-----------------+-----------------+----------------+---------------------------------+
    |  cvode.precond_type  |  N/A            |  N/A            |  N/A           |  cuSparse_simplified_AJacobian  |
    +----------------------+-----------------+-----------------+----------------+---------------------------------+
    |  Run time            |  13s            |  20s            |  19s           |  36s                            |
    +----------------------+-----------------+-----------------+----------------+---------------------------------+



Current Limitations
-------------------

The current GPU implementation of CVODE relies on the launch of many kernels from the host. As such, a CVODE instance does not live *directly* on the GPU; rather, the user is in charge of identifying and delegating computation-intensive part of the RHS, Jacobian evaluation, etc.
The current implementation thus suffers from the cost of data movement, and parallelization is limited due to required device synchronizations within CVODE.

.. [#Foot1] NOTE that only one cell at a time should be integrated with CVODE right now. The vectorized version on CPU is still WIP and not properly implemented for all linear solvers so that no computational gain should be expected from solving several cells at a time.
