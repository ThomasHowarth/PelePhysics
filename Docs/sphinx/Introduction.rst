.. highlight:: rst

.. _sec:subsWD:

************
Introduction
************

`PelePhysics` is a repository of physics databases and implementation code for use with the `Pele suite of codes <https://Pele-Suite.github.io>`_,
primarily:

- `PeleC <https://Pele-Suite.github.io/PeleC/>`_, an adaptive-mesh compressible hydrodynamics code for reacting flows.
- `PeleLMeX <https://Pele-Suite.github.io/PeleLMeX/>`_, an adaptive-mesh low-Mach number hydrodynamics code for reacting flows.

Both codes, and `PelePhysics` itself, rely on `AMReX <https://amrex-codes.github.io/amrex>`_, a software framework that provides the data structures
and enables massive parallelization.

.. |a| image:: ./Visualization/PeleSuite.png

.. table::
   :align: center

   +-----+
   | |a| |
   +-----+

`PelePhysics` is primarily intended as a library: it does not solve for a flow by itself, but provides the thermodynamic, transport, chemistry and
multi-physics models that the flow solvers need, such that these models are implemented once and shared across the suite. In particular, the choice
of chemistry and transport models as well as the associated functions and capabilities are managed in `PelePhysics`. A number of stand-alone programs are also
provided to test the models independently of the flow solvers, see :ref:`sec:BuildingRunning`.

This section gives an overview of the content of the repository and of the main design choices of the code. Relevant references to the other sections of this
documentation are given for an in-depth discussion.

Organization of the repository
==============================

The documentation is organized to follow the layout of the repository, which contains the following directories:

* ``Source``: the C++ implementation of the physics modules and utilities described below.
* ``Mechanisms``: the chemical mechanisms shipped with `PelePhysics`. For each mechanism, the directory contains the `CEPTR`-generated source code
  (``mechanism.H`` and ``mechanism.cpp``) as well as the Cantera and, where available, CHEMKIN format data files it was generated from.
* ``Support``: stand-alone tools to aid in `Pele` simulation workflows, see :ref:`sec:Support`.
* ``Testing``: stand-alone programs for exercising various aspects of the physics modules (``Testing/Exec``) and the configuration of the regression test suite (``Testing/Regression``).
* ``Submodules``: the two required dependencies, `AMReX` and `SUNDIALS`, included as git submodules, see :ref:`sec:GetCVODE`.
* ``ThirdParty``: the makefiles used to build `SUNDIALS` and to download and build the optional dependencies (`SuiteSparse` and `MAGMA`). The resulting libraries are installed in this directory.
* ``Docs``: the source files of this documentation.

Physics modules
---------------

The physics modules are found in the ``Source`` directory:

* ``Eos``: the equations of state, including multi-species ideal gas and real gas options, and the associated thermodynamic functions. See :ref:`sec:eos` and :ref:`sec:thermodynamics`.
* ``Transport``: the routines to evaluate the transport properties of multi-species mixtures. See :ref:`sec:transport`.
* ``Reactions``: the routines for the substepped integration of stiff chemical systems, including with CVODE. See :ref:`sec:chemistry`.
* ``Soot``: an implementation of the Hybrid Method of Moments soot model, formerly part of `PeleMP <https://github.com/Pele-Suite/PeleMP>`_. See :ref:`Soot`.
* ``Spray``: a Lagrangian spray droplet library, formerly part of `PeleMP <https://github.com/Pele-Suite/PeleMP>`_. See :ref:`Spray`.
* ``Radiation``: a radiative heat transfer model based on the spherical harmonics method, formerly `PeleRad <https://github.com/Pele-Suite/PeleRad>`_.

Utilities
---------

The ``Source/Utility`` directory contains several modules for tasks that are shared across the `Pele` codes, see :ref:`sec:Utility`:

* ``PMF``: initialization of a simulation from a one-dimensional premixed flame solution.
* ``TurbInflow``: turbulent inflow boundary conditions from precomputed turbulence files, see :ref:`sec_turbinflow`.
* ``TurbForcing``: forcing to generate and maintain homogeneous isotropic turbulence, see :ref:`sec_turbforce`.
* ``MeshMap``: evaluation of the grid mappings used with non-uniform and mapped grids.
* ``PltFileManager``: data structures to handle the data read from plot files.
* ``Diagnostics``: runtime diagnostics, see :ref:`sec_diagnostics`.
* ``Filter``: spatial filtering of data.
* ``BlackBoxFunction``: generic functions evaluated from tabulated data or neural networks, which are used by the manifold models.
* ``Utilities``: basic functions, including unit conversions.

Support codes
-------------

The ``Support`` directory contains the following tools:

* ``ceptr``: Chemistry Evaluation for Pele Through Recasting (`CEPTR`), a tool that converts chemical mechanisms from the `Cantera <https://cantera.org>`_ yaml format to
  efficient source code files for use in the `Pele` codes. With `CEPTR`, a unique set of chemistry files is generated for a specific kinetic scheme.
  The mechanisms found in the ``Mechanisms`` directory were all generated in this manner, and instructions to generate your own chemistry files are provided in
  :ref:`sec:ceptr`. `CEPTR` can also generate reduced mechanisms based on the quasi-steady state assumption, see :ref:`sec:QSSAssumption` and :ref:`sec:tutorials`.
* ``TurbInflowGenerator``: generation of the turbulence files used by the turbulent inflow utility, see :ref:`sec_turbfile`.
* ``MechanismPAH``: a tool to add a PAH module to an existing mechanism, as needed for the soot model.
* ``liquidProp``: scripts to fit the liquid properties required by the spray model.
* ``CMLM``: sample table and neural network data files for the manifold models.

Code design
===========

Language and namespaces
-----------------------

`PelePhysics` is written in C++, and all the source code is contained in the ``pele::physics`` namespace. Each module has its own nested namespace, such as
``pele::physics::eos``, ``pele::physics::transport``, ``pele::physics::reactions`` or ``pele::physics::utilities``.

Compile-time selection of the physics models
--------------------------------------------

The equation of state, the transport model and the chemical mechanism are selected at compile time. When using GNUmake, this is done by setting the following
parameters in the ``GNUmakefile``: ::

    Eos_Model       := Fuego
    Transport_Model := Simple
    Chemistry_Model := drm19

where ``Chemistry_Model`` is the name of one of the directories of ``Mechanisms``. The values of ``Eos_Model`` and ``Transport_Model`` are converted into
preprocessor definitions, which are in turn used to define the ``pele::physics::EosType`` and ``pele::physics::TransportType`` types. These two types are
combined in ``PelePhysics.H`` into a single type: ::

    using PhysicsType = PelePhysics<EosType, TransportType>;

from which the equation of state and transport objects are obtained: ::

    auto eos = pele::physics::PhysicsType::eos();
    auto trans = pele::physics::PhysicsType::transport();

Because the models are known at compile time, the functions of these objects can be inlined and called from within the loops over the cells, on CPU as well as on GPU.

Not all the combinations of models are meaningful. The valid combinations are listed in ``PelePhysicsConstraints.H``, and attempting to use an invalid one
leads to a compilation error with the message ``Invalid physics combination attempted``:

* The ``GammaLaw`` equation of state can only be used with ``Constant`` transport, and requires ``Chemistry_Model`` to be ``Null``.
* The ``Soave-Redlich-Kwong`` equation of state can only be used with ``Simple`` transport.
* The ``Manifold`` equation of state can only be used with ``Manifold`` or ``Constant`` transport, and ``Manifold`` transport can only be used with the ``Manifold`` equation of state.

The soot and spray modules are also activated at compile time, see :ref:`Soot` and :ref:`Spray`.

Run-time selection of the chemistry integrator
----------------------------------------------

In contrast with the models above, the integrator used to advance the chemical system is selected at run time. The integrators all derive from the ``ReactorBase`` class and
are registered by name using the factory defined in ``Factory.H``, such that an integrator can be created from a keyword read in the input file: ::

    auto reactor = pele::physics::reactions::ReactorBase::create("ReactorCvode");

The available integrators are ``ReactorCvode``, ``ReactorArkode``, ``ReactorRK64``, ``ReactorBDF`` and ``ReactorNull``, the latter being used when the reactions are not integrated.
The name of the input file keyword is defined by the code using `PelePhysics`; the stand-alone test programs use ``chem_integrator``. The options specific
to the integrators, such as the tolerances or the linear solver used by CVODE, are also specified in the input file. See :ref:`sec:chemistry` for additional details.

The same factory is used for the run-time selection of the functions of the ``BlackBoxFunction`` utility.

Performance portability
-----------------------

`PelePhysics` relies on `AMReX` for performance portability. The functions evaluating the thermodynamic, transport and chemistry properties of a single cell are marked such that
they can be compiled for the host as well as for the device, and are called by the flow solvers from within `AMReX` parallel loops. The same source code can therefore be used on CPU,
with MPI and OpenMP, and on GPU, with CUDA, HIP or SYCL. For the chemistry integration on GPU, many cells are grouped together and integrated at once.

Units
-----

The routines generated by `CEPTR` and the equation of state and transport functions operate in CGS units, with the temperature in K.
Because `PeleLMeX` uses MKS units, conversion functions between the two unit systems are provided in the ``pele::physics::utilities`` namespace, see :ref:`sec:Utility`.

Build systems
-------------

`PelePhysics` is built as part of the code that uses it rather than as a stand-alone library. For GNUmake, ``Testing/Exec/Make.PelePhysics`` gathers the source files,
the include paths and the definitions corresponding to the selected models, and is included by the ``GNUmakefile`` of each of the stand-alone test programs. The dependencies
are built separately with the ``make TPL`` command, as described in :ref:`sec:BuildingRunning`. When `PeleC` or `PeleLMeX` are built with CMake, the `PelePhysics` source files are compiled
as part of the build of these codes; `PelePhysics` does not provide a top-level ``CMakeLists.txt``.
