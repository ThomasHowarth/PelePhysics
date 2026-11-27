.. highlight:: rst

.. _sec:thermodynamics:

**************
Thermodynamics
**************

PelePhysics supports the computation of thermodynamic properties for different species. The thermodynamic property models currently supported are the NASA7 and NASA9 polynomial parameterizations. The reader is referred to the many discussions about these models to learn more specifics, including `Cantera's discussion about thermodynamics properties <https://cantera.org/science/species-thermo.html>`_.

The polynomial coefficients are read by CEPTR from the Cantera-format mechanism file and are written into the generated mechanism code (see :ref:`sec:ceptr`), so the thermodynamic model is part of the ``Chemistry_Model`` selected at compile time. For each species, the generated code evaluates the polynomials of the temperature range that contains the temperature passed as an argument: two ranges for NASA7, and one or several ranges for NASA9.

The generated ``mechanism.H`` provides the following routines, which return the nondimensional properties of all the species at a given temperature:

* ``cp_R`` and ``cv_R``: specific heats at constant pressure and volume, :math:`c_p/R` and :math:`c_v/R`
* ``speciesEnthalpy`` and ``speciesInternalEnergy``: enthalpy and internal energy, :math:`h/(RT)` and :math:`e/(RT)`
* ``speciesEntropy``: standard-state entropy, :math:`s/R`
* ``gibbs`` and ``helmholtz``: Gibbs and Helmholtz free energies, :math:`g/(RT)` and :math:`a/(RT)`
* ``dcvpRdT``: temperature derivative of the specific heats, :math:`d(c_p/R)/dT`

Dimensional properties are obtained from routines that follow the CHEMKIN-II naming conventions, where the suffix indicates whether the property is for the species or the mixture, and whether it is in molar or mass units. For instance, ``CKCPMS`` and ``CKHMS`` return the species specific heats at constant pressure and enthalpies in mass units, ``CKCPML`` and ``CKHML`` return the same quantities in molar units, and ``CKCPBS`` and ``CKHBMS`` return the mixture-averaged specific heat and enthalpy in mass units for a given set of mass fractions. The routines ``GET_T_GIVEN_EY`` and ``GET_T_GIVEN_HY`` compute the temperature from the mixture internal energy or enthalpy and the mass fractions. All these routines use CGS units.

In most cases these routines are not called directly, but through the equation of state (see :ref:`sec:eos`), which adds the departure from the ideal gas behavior when the Soave-Redlich-Kwong EOS is used.
