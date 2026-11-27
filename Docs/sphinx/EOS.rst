.. highlight:: rst

.. _sec:eos:

*****************
Equation of State
*****************

PelePhysics allows the user to use different equation of state (EOS) as the constitutive equation and close the compressible Navier-Stokes system of equations. All the routines needed to fully define an EOS are implemented through PelePhysics module. Available models include:

* A simple ``GammaLaw`` model for a single component perfect gas
* An ideal gas mixture model (similar to the CHEMKIN-II approach) labeled ``Fuego``
* The ``Soave-Redlich-Kwong`` cubic equation of state, labeled ``SRK``
* A ``Manifold`` model, where the thermochemical state is obtained from a pretabulated or neural network manifold rather than from transported species

The EOS implementations can be found in ``PelePhysics/Source/Eos``. The choice between these EOS models is made at compile time. When using GNUmake, this is done by setting the ``Eos_Model`` parameter in the ``GNUmakefile`` to ``GammaLaw``, ``Fuego``, ``SRK`` or ``Manifold``.

The following sections first describe the interface that is shared by all the EOS models, and then each of the models. The implementation of Soave-Redlich-Kwong, a non-ideal cubic EOS, for a general mixture of species is described in the most detail. For an advanced user interested in implementing a new EOS this chapter should provide a good starting point.

.. note::  For the flow solvers in the Pele suite, the SRK EOS is presently only supported in PeleC, and not PeleLM(eX).

.. _sec:eos:interface:

EOS interface
=============

Each EOS is a ``struct`` in the ``pele::physics::eos`` namespace (``GammaLaw``, ``Fuego``, ``SRK`` and ``Manifold``), defined in the header file of the same name in ``PelePhysics/Source/Eos``. The compile-time selection defines one of the macros ``USE_GAMMALAW_EOS``, ``USE_FUEGO_EOS``, ``USE_SRK_EOS`` or ``USE_MANIFOLD_EOS``, from which ``EOS.H`` defines the alias ``pele::physics::EosType``. The alias ``pele::physics::PhysicsType``, defined in ``PelePhysics.H``, combines the selected EOS and transport models, and is the usual way to obtain an EOS object: ::

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real P = 0.0;
  eos.RTY2P(rho, T, Y, P);

The EOS objects are lightweight and are meant to be created where they are needed, including inside GPU kernels.

The member functions are named according to their inputs and outputs, as ``<inputs>2<outputs>``, with the arguments passed in the same order as they appear in the name. The symbols used in the names are:

* ``R``: density
* ``T``: temperature
* ``P``: pressure
* ``Y``: species mass fractions, ``X``: species mole fractions, ``C``: species molar concentrations
* ``E``: mixture specific internal energy, ``Ei``: species specific internal energies
* ``H``: mixture specific enthalpy, ``Hi``: species specific enthalpies
* ``Cv``, ``Cp``: mixture specific heats at constant volume and pressure, ``Cpi``: species specific heats at constant pressure
* ``Cs``: speed of sound
* ``G``: ratio of specific heats (for SRK, the equivalent quantity that relates the speed of sound to :math:`p / \rho`)
* ``WBAR``: mean molecular weight
* ``WDOT``: species mass production rates
* ``JAC``: Jacobian of the chemical system
* ``Z``: elemental mass fractions
* ``RY``: species partial densities

For instance, ``RTY2P`` computes the pressure from the density, temperature and mass fractions, ``REY2T`` computes the temperature from the density, internal energy and mass fractions, and ``RY2RRinvY`` computes the density, its inverse and the mass fractions from the partial densities. The functions ``molecular_weight`` and ``inv_molecular_weight`` return the species molecular weights and their inverses. All quantities are in CGS units.

The same set of functions is declared by every EOS, so that code calling the EOS compiles for any of the models. However, not every function is meaningful for every model: where this is the case, the function exists but calls ``amrex::Error`` when executed. For example, the functions that do not take the density as an input (``TY2Cp``, ``TY2H``, ``T2Hi``, ...) are not available for SRK, and ``RTY2JAC`` is only implemented for Fuego.

A few free functions templated on the EOS type are also provided in ``EOS.H``: ``speciesNames`` returns the names of the transported species, and ``atomic_weightsCHON`` and ``element_compositionCHON`` return the atomic weights and the species elemental compositions for the elements C, H, O and N (these two are not supported for GammaLaw and Manifold).

EOS parameters
--------------

EOS models that need data at run time receive it through the ``EosParm<EOSType>`` struct defined in ``EosParams.H``. This struct is empty for Fuego and SRK. For GammaLaw it contains the value of :math:`\gamma`, and for Manifold it contains the manifold data and the indices of the manifold quantities used by the EOS (see :ref:`sec:eos:manifold`). The parameters are managed by the ``PeleParams`` container, which reads the relevant input file options when ``initialize()`` is called and provides a pointer to the parameters that can be used on the GPU: ::

  pele::physics::PeleParams<
    pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>> eos_parms;
  eos_parms.initialize();
  auto const* leosparm = eos_parms.device_parm();
  ...
  auto eos = pele::physics::PhysicsType::eos(leosparm);
  ...
  eos_parms.deallocate();

The Manifold EOS must be constructed with a pointer to its parameters in this way. The other models can also be constructed without any argument, in which case GammaLaw uses the default value of :math:`\gamma`.

GammaLaw
========

This EOS corresponds to a single component perfect gas, i.e. a gas that follows :math:`p = \rho \hat{R} T` and :math:`e = c_v T`, with constant :math:`c_v`, which together imply that :math:`p = (\gamma - 1)\rho e`, where :math:`\gamma = c_p / c_v`. The values of the relevant physical properties (:math:`\gamma = 1.4`, :math:`MW=28.97` g/mol) are chosen to correspond to air at standard conditions and are set in ``PelePhysics/Source/PhysicsConstants.H``.

When the GammaLaw EOS is used, the ``Chemistry_Model`` should be set to ``Null``.

The value of :math:`\gamma` can be changed in the input file using: ::

  eos.gamma = 1.4

A warning is issued if this value differs from the compile-time default, as the run-time specification of :math:`\gamma` may not be supported everywhere in the codes that use PelePhysics.

Fuego
=====

This is a multi-component ideal gas EOS to be used for multi-component (reacting) calculations with any ``Chemistry_Model`` besides ``Null``. The gas mixture follows :math:`p = \rho \hat{R} T`, but with :math:`h(T) \equiv \sum Y_i h_i(T)` and :math:`h_i(T)` computed based on `NASA polynomials <https://ntrs.nasa.gov/citations/20020085330>`_ that are included in the ``Chemistry_Model`` (see :ref:`sec:thermodynamics`). All the functions of this EOS are wrappers around the routines of the mechanism files generated by CEPTR.

Soave-Redlich-Kwong (SRK)
=========================

The cubic model is built on top of the ideal gas models. It should be used for multi-component (reacting) calculations where the pressure is sufficiently high that real gas effects are important; there is a significant computational cost penalty relative ideal gas calculations. Note that the function to compute the analytical Jacobian of the chemical system using the SRK EOS has not been implemented, so this EOS is not compatible with reaction integrators that rely on computing the analytical Jacobian. The additional parameters needed for each species (attraction and repulsion terms, critical temperature and acentric factor dependence) are computed by CEPTR when the mechanism files are generated and are provided by the ``GET_CRITPARAMS_SRK`` routine of the mechanism. They are based on critical properties tabulated in CEPTR when these are available for a species, or else are inferred from the transport data of the mechanism, as described below.

SRK EOS as a function of Pressure (p), Temperature(T), and :math:`\tau` (specific volume) is given by

.. math::
   p = R T \sum \frac{Y_k}{W_k} \frac{1}{\tau - b_m} - \frac{a_m}{\tau(\tau + b_m)}

where :math:`Y_k` are species mass fractions, :math:`R` is the universal gas constant, and
:math:`b_m` and :math:`a_m` are mixture repulsion and attraction terms, respectively.

Mixing rules
------------

For a mixture of species, the following mixing rules are used to compute :math:`b_m` and :math:`a_m`.

.. math::
   a_m = \sum_{ij} Y_i Y_j \alpha_i \alpha_j \;\;\;  b_m = \sum_k Y_k b_k

where :math:`b_i` and :math:`a_i` for each species is defined using critical pressure and temperature.

.. math::
   a_i(T) = 0.42748 \frac{\left(R T_{c,i} \right)^2}{W_i^2 p_{c,i}} \bar{a}_i \left(T/T_{c,i}\right) \;\;\;
   b_i = 0.08664 \frac{R T_{c,i}}{W_i p_{c,i}}

where

.. math::
   \bar{a}_i (T/T_{c,i}) = \left(1 + \mathcal{A} \left[ f\left( \omega_i \right) \left(1-\sqrt{T/T_{c,i}} \right ) \right] \right)^2

where :math:`\omega_i` are the accentric factors and

.. math::
   f\left( \omega_i \right) = 0.48508 + 1.5517 \omega_i - 0.151613 \omega_{i}^2

For chemically unstable species such as radicals, critical temperatures and pressures are not available.
For species where critical properties are not available, we use the Lennard-Jones potential for that species to construct attractive and repulsive coefficients, and the acentric factor is taken to be zero.

.. math::
   T_{c,i} = 1.316 \frac{\epsilon_i}{k_b} \;\;\;  a_i(T_{c,i}) = 5.55 \frac{\epsilon_i \sigma_i^3}{m_i^2} \;\;\;
   \mathrm{and} \;\;\; b_i = 0.855 \frac{\sigma_i^3}{m_i}

where :math:`\sigma_i`, :math:`\epsilon_i` are the Lennard-Jones potential molecular diameter and well-depth, respectively,
:math:`m_i` the molecular mass, and :math:`k_b` is Boltzmann's constant.

In terms of implementation, a member function called ``MixingRuleAmBm`` can be found in the SRK EOS implementation (``SRK.H``). The following code block shows this function, which receives species mass fractions and temperature as input. The outputs of this function are :math:`b_m` and :math:`a_m`. The species data ``Fomega``, ``sqrtOneOverTc``, ``sqrtAsti`` and ``Bi`` are the quantities :math:`f(\omega_i)`, :math:`1/\sqrt{T_{c,i}}`, :math:`\sqrt{a_i(T_{c,i})}` and :math:`b_i` obtained from the mechanism.

.. code-block:: c++

   void MixingRuleAmBm(
     const amrex::Real T,
     const amrex::Real Y[NUM_SPECIES],
     amrex::Real& am,
     amrex::Real& bm)
   {
     am = 0.0;
     bm = 0.0;
     const amrex::Real sqrtT = std::sqrt(T);
     amrex::Real amloc[NUM_SPECIES];

     for (int ii = 0; ii < NUM_SPECIES; ii++) {
       amloc[ii] =
         (1.0 + Fomega[ii] * (1.0 - sqrtT * sqrtOneOverTc[ii])) * sqrtAsti[ii];
       bm += Y[ii] * Bi[ii];
     }

     for (int ii = 0; ii < NUM_SPECIES; ii++) {
       for (int jj = 0; jj < NUM_SPECIES; jj++) {
         am += Y[ii] * Y[jj] * amloc[ii] * amloc[jj];
       }
     }
   }

Thermodynamic Properties
------------------------

Most of the thermodynamic properties can be calculated from the equation of state and involve derivatives of various thermodynamic quantities and of EOS parameters. In the following, some of these thermodynamic properties for SRK and the corresponding member functions of the ``SRK`` struct are presented.

Specific heat
^^^^^^^^^^^^^

For computing mixture specific heat at constant volume and pressure, the ideal gas contribution and the departure from the ideal gas are computed. Specific heat at constant volume can be computed using the following

.. math::
   c_v = \left( \frac{\partial e_m}{\partial T}\right)_{\tau,Y}

For SRK EOS, the formula for :math:`c_v` reduces to

.. math::
   c_v = c_v^{id} + T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau})

where :math:`c_v^{id}` is the ideal gas specific heat at constant volume. Mixture specific heat at constant volume is implemented through the function ``RTY2Cv``

.. code-block:: c++

   void RTY2Cv(
     const amrex::Real R,
     const amrex::Real T,
     const amrex::Real Y[NUM_SPECIES],
     amrex::Real& Cv)
   {
     amrex::Real am, bm, d2AmdT2;
     MixingRuleAmBm(T, Y, am, bm);
     // Second derivative of the EOS am w.r.t. temperature
     Calc_d2AmdT2(T, Y, d2AmdT2);
     // Ideal gas specific heat at constant volume
     CKCVBS(T, Y, Cv);
     // Real gas specific heat at constant volume
     Cv += T * d2AmdT2 * (1.0 / bm) * log1p(bm * R);
   }

Specific heat at constant pressure is given by

.. math::

   c_p = \left( \frac{\partial h_m}{\partial T}\right)_{p,Y}   \;\; \\
   c_p =  \frac{\partial h_m}{\partial T} - \frac {\frac{\partial h}{\partial \tau}} {\frac{\partial p}{\partial \tau}} \frac{\partial p}{\partial T}

where all the derivatives in the above expression for SRK EOS are given by

.. math::

   \frac{\partial p}{\partial T} = \sum Y_k / W_k  \frac{R}{\tau-b_m} - \frac{\partial a_m}{\partial T} \frac{1}{\tau(\tau +b_m)} \\
   \frac{\partial p}{\partial \tau} = -\sum Y_k / W_k  \frac{R T}{(\tau-b_m)^2} + \frac{a_m (2 \tau + b_m)}{[\tau(\tau +b_m)]^2} \\
   \frac{\partial h_m}{\partial \tau} = -\left(T \frac{\partial a_m}{\partial T}  - a_m \right) \frac{1}{\tau(\tau+b_m)} + \frac{a_m}{(\tau+b_m)^2} -\sum Y_k / W_k  \frac{R T b_m}{(\tau-b_m)^2}  \\
   \frac{\partial h_m}{\partial T} = c_p^{id} +T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau}) - \frac{\partial a_m}{\partial T} \frac{1}{\tau+b_m} +\sum Y_k / W_k  \frac{R b_m}{\tau-b_m}

These expressions are evaluated in the function ``RTY2Cp``, where the ideal gas specific heat at constant pressure :math:`c_p^{id}` is obtained from the ``CKCPBS`` routine of the mechanism and the temperature derivatives of :math:`a_m` are computed by ``Calc_Am_and_derivs``.

Internal energy and Enthalpy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similarly mixture internal energy for SRK EOS is given by

.. math::
   e_m = \sum_k Y_k e_k^{id} + \left( T  \frac{\partial a_m}{\partial T}  - a_m \right)\frac{1}{b_m} ln \left( 1 + \frac{b_m}{\tau}\right)

and mixture enthalpy :math:`h_m = e + p \tau`

.. math::
   h_m = \sum_k Y_k h_k^{id} + \left ( T \frac{\partial a_m}{\partial T} - a_m \right) \frac{1}{b_m} \ln \left( 1 + \frac{b_m}{\tau}\right) + R T \sum \frac{Y_k}{W_k} \frac{b_m}{\tau -b_m} - \frac{a_m}{\tau + b_m}

and the implementations can be found in the functions ``RTY2E`` and ``RTY2H``.

Speed of Sound
^^^^^^^^^^^^^^

The sound speed for SRK EOS is given by

.. math::

   a^2 = -\frac{c_p}{c_v} \tau^2  \frac{\partial p}{\partial \tau}

which is computed by the function ``RTY2Cs``.

Species enthalpy
^^^^^^^^^^^^^^^^

For computation of kinetics and transport fluxes we will also need the species partial enthalpies and the chemical potential.  The species enthalpies for SRK EOS, computed by the function ``RTY2Hi``, are given by

.. math::

   h_k = \frac{\partial h_m}{\partial Y_k } - \frac {\frac{\partial h}{\partial \tau}} {\frac{\partial p}{\partial \tau}} \frac{\partial p}{\partial Y_k}

where

.. math::
   \frac{\partial h_m}{\partial Y_k } &=  h_k^{id} + (T \frac{\partial^2 a_m}{\partial T \partial Y_k}  - \frac{\partial a_m }{\partial Y_k}) \frac{1}{b_m} \ln\left(1+ \frac{b_m}{\tau}\right) \\&-\left(T \frac{\partial a_m}{\partial T}  - a_m \right) \left[ \frac{1}{b_m^2} \ln\left(1+ \frac{b_m}{\tau}\right) - \frac{1}{b_m(\tau+b_m)} \right ] \frac{\partial b_m}{\partial Y_k} \nonumber \\&+ \frac{a_m}{(\tau+b_m)^2}  \frac{\partial b_m}{\partial Y_k} - \frac{1}{\tau+b_m}  \frac{\partial a_m}{\partial Y_k} + 1 / W_k  \frac{R T b_m}{\tau-b_m}\\&+\sum_i \frac{Y_i}{W_i} R T \left( \frac{1}{\tau -b_m} + \frac{b_m}{(\tau-b_m)^2} \right) \frac{ \partial b_m}{\partial Y_k}

.. math::

   \frac{\partial p}{\partial Y_k} &= R T \frac{1}{W_k} \frac{1}{\tau - b_m} - \frac{\partial a_m}{\partial Y_k} \frac{1}{\tau(\tau + b_m)} \\&+\left(R T \sum \frac{Y_i}{W_i} \frac{1}{(\tau - b_m)^2} + \frac{a_m}{\tau(\tau + b_m)^2} \right ) \frac{\partial b_m}{\partial Y_k}

Chemical potential
^^^^^^^^^^^^^^^^^^

The chemical potentials are the derivative of the free energy with respect to composition.  Here the free energy :math:`f` is given by

.. math::
   f &= \sum_i Y_i (e_i^{id} - T s_i^{id,*}) +  \sum_i \frac{Y_i R T}{W_i} ln (\frac{Y_i R T}{W_i \tau p^{st}})  \nonumber \\ &+ \sum_i \frac{Y_i R T}{W_i} ln (\frac{\tau}{\tau-b_m}) -  a_m \frac{1}{b_m}ln (1+ \frac{b_m}{\tau})  \nonumber \\ &= \sum_i Y_i (e_i^{id} - T s_i^{id,*}) +  \sum_i \frac{Y_i R T}{W_i} ln (\frac{Y_i R T}{W_i (\tau-b_m) p^{st}} )- a_m \frac{1}{b_m} ln (1+ \frac{b_m}{\tau})  \nonumber

Then

.. math::

   \mu_k &= \frac{\partial f}{\partial Y_k} = e_k^{id} - T s_k^{id,*}  + \frac{RT}{W_k} ln (\frac{Y_k R T}{W_k (\tau-b_m) p^{st}}) + \frac{RT}{W_k} +  \frac{RT}{\bar{W}} \frac{1}{\tau-b_m} \frac {\partial b_m}{\partial Y_k} \nonumber \\
   &- \frac{1}{b_m} ln(1 + \frac{b_m}{\tau}) \frac{\partial a_m}{\partial Y_k}+ \frac{a_m}{b_m^2} ln(1 + \frac{b_m}{\tau}) \frac{\partial b_m}{\partial Y_k}- \frac{a_m}{b_m} \frac{1}{\tau+b_m} \frac{\partial b_m}{\partial Y_k}

Other primitive variable derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Godunov (FV) algorithm also needs some derivatives to express source terms in terms of primitive variables. In particular one needs

.. math::

   \left . \frac{\partial p}{\partial \rho} \right|_{e,Y} =-\tau^2 \left( \frac{\partial p}{\partial \tau}- \frac {\frac{\partial e}{\partial \tau}} {\frac{\partial e}{\partial T}} \frac{\partial p}{\partial T} \right )

and

.. math::

   \left . \frac{\partial p}{\partial e} \right|_{\rho,Y} = \frac{1}{c_v} \frac{\partial p}{\partial T}

All of the terms needed to evaluate this quantity are known except for

.. math::

   \frac{\partial e}{\partial \tau} = \frac{1}{\tau ( \tau + b_m)} \left( a_m - T  \frac{\partial a_m}{\partial T}  \right) \;\; .

.. _sec:eos:manifold:

Manifold
========

The Manifold EOS is used for combustion models where the thermochemical state is assumed to lie on a low-dimensional manifold. Rather than transporting the species of a chemical mechanism, a small number of manifold parameters :math:`\xi_k` (for instance a mixture fraction and a progress variable) are transported, and all the other quantities, such as the density, the temperature and the source terms of the manifold parameters, are functions of these parameters, :math:`\phi = \phi(\xi_k)`. These functions are provided either as a table, which is interpolated, or as a neural network. The implementation of the tables and neural networks is not part of the EOS: it can be found in ``PelePhysics/Source/Utility/BlackBoxFunction``. Sample table and network files are provided in ``PelePhysics/Support/CMLM`` (see :ref:`sec_cmlm`).

.. note:: The Manifold EOS is based on a nominal pressure, as is appropriate for low Mach number flow solvers. Several functions of the EOS interface are placeholders, as described below.

State
-----

The Manifold EOS is used with ``Chemistry_Model = Null``. In that case the number of "species" is ``NUM_SPECIES = MANIFOLD_DIM + 1``, where ``MANIFOLD_DIM`` is the number of manifold parameters, set at compile time. The first ``MANIFOLD_DIM`` entries of the species state are the manifold parameters, in the order of the dimensions of the table or of the inputs of the network, and the last entry corresponds to the density. Wherever the EOS interface takes the mass fractions ``Y`` as an argument, the Manifold EOS therefore expects the manifold parameters. The function ``speciesNames`` returns the names of the manifold parameters, as read from the table or network metadata, followed by ``XRHO``.

The following quantities are looked up in the manifold. The names are those of the variables that must be present in the table or among the outputs of the network:

* ``RHO``, ``lnRHO`` (or ``logRHO``) or ``invRHO``: the density, its natural logarithm or its inverse, depending on the value of ``manifold.density_lookup_type``. This is used by ``PYT2R``, ``RTY2P`` and ``RY2dRdY``.
* ``T``: the temperature, returned by ``EY2T``, ``REY2T``, ``HY2T`` and ``RHY2T`` when ``manifold.compute_temperature`` is true. Otherwise these functions return a temperature of 1. The variable must be present in both cases.
* ``SRC_<name>``: the source term of the manifold parameter ``<name>`` (if the name of the parameter starts with ``Y-``, this prefix is removed), returned by ``RTY2WDOT``. If the variable is not found, the source term for that parameter is assumed to be zero. The source term for the density entry is always zero.
* ``WBAR``: the mean molecular weight, returned by ``Y2WBAR``. This is optional, unless the code is compiled with sprays.
* ``CP``: the specific heat at constant pressure, returned by ``TY2Cp`` and ``RTY2Cp``. This is optional, unless the code is compiled with sprays. When it is not present, these functions return a large placeholder value.
* ``CP_fuel``: optional, unless the code is compiled with sprays.

Any other variable of the manifold can be evaluated using the function ``Y2GenericManifoldData``, which takes a list of variable indices. The function ``varNames`` returns the names of all the variables.

The density is obtained directly from the manifold with ``PYT2R``, irrespective of the pressure and temperature arguments. The function ``RTY2P`` returns the nominal pressure rescaled by the ratio of the density passed as an argument to the density from the manifold:

.. math::
   p = p_{nom} \frac{\rho}{\rho(\xi_k)}

so that the nominal pressure is recovered when the two densities are consistent. The function ``RY2dRdY`` returns the derivatives of the density with respect to the manifold parameters.

A manifold parameter whose name ends with ``VAR`` is assumed to be the variance of the parameter with the same name without this suffix (for instance, ``ZMIXVAR`` is the variance of ``ZMIX``), which must then exist. This information is stored in the ``is_variance_of`` array of the EOS parameters.

The quantities that are not defined for a manifold model are not available. The functions for the internal energies, the species enthalpies and specific heats, the specific heat at constant volume, the chemical Jacobian, the species concentrations and the conversion from mass to mole fractions call ``amrex::Error``. A few functions return placeholder values instead: ``T2Hi`` and ``TY2H`` return zero, and ``X2Y`` returns its input unchanged.

Compiling
---------

When using GNUmake, the following parameters are set in the ``GNUmakefile``: ::

  Eos_Model       = Manifold
  Chemistry_Model = Null
  Transport_Model = Manifold
  Manifold_Dim    = 2
  Manifold_Type   = Table

``Manifold_Dim`` is the number of manifold parameters. It must match the number of dimensions of the manifold data (table or network) that is loaded at run time, or the initialization will fail. ``Manifold_Type`` can be set to ``Table`` or ``Network``. If it is not set to either of these values, the selection between table and network is made at run time, which has a significant performance cost. When compiling with SYCL, only ``Table`` is supported. The Manifold EOS can be used with the ``Manifold`` and ``Constant`` transport models (see :ref:`sec:transport`).

Input file options
------------------

The input file options for the Manifold EOS are provided below, with their default values. ::

  manifold.model                = [Must be User Defined]  # Table or NeuralNet
  manifold.nominal_pressure_cgs = [Must be User Defined]  # Nominal pressure, in CGS units
  manifold.compute_temperature  = [Must be User Defined]  # Look up the temperature in the manifold
  manifold.density_lookup_type  = linear                  # linear (RHO), log (lnRHO or logRHO) or inverse (invRHO)
  manifold.has_mani_src         = 1                       # Source terms of the manifold parameters are in the manifold
  manifold.has_species_mw       = 0                       # Read species molecular weights
  manifold.metadata_filename    = [Optional]              # File with additional input file options
  manifold.v                    = 1                       # Verbosity level

Setting ``manifold.has_mani_src = 0`` is not supported at present and leads to an error. If ``manifold.has_species_mw`` is true, the variables of the manifold named ``Y-<species>`` are taken to be species mass fractions, and a molecular weight must be provided for each of them with ``manifold.<species>_mw``. The corresponding mass fractions and molecular weights can then be obtained with the functions ``Y2ChemSpecies`` and ``molecular_weight``, and the names with ``chemSpeciesNames``. When the code is compiled with sprays, only the species involved in the spray phase change are considered. If ``manifold.metadata_filename`` is provided, the file is read as an additional set of input file options before the other options are queried.

The table or the network itself is specified using the same ``manifold`` prefix. For a table: ::

  manifold.model          = Table
  manifold.table.filename = ../../../Support/CMLM/peletable_2dim_10grid
  manifold.table.v        = 2

and for a neural network: ::

  manifold.model                 = NeuralNet
  manifold.network.filename      = ../../../Support/CMLM/fgm_net.pnn
  manifold.network.info_filename = ../../../Support/CMLM/fgm_net_info.txt
  manifold.network.cmlm_net      = 1
  manifold.network.v             = 2

These examples are taken from the ``EosEval`` test case (``PelePhysics/Testing/Exec/EosEval``), which can be compiled with the Manifold EOS by setting ``USE_MANIFOLD = TRUE`` in its ``GNUmakefile``.
