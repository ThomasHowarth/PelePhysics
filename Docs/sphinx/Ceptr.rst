.. highlight:: rst

.. _sec:ceptr:

CEPTR: Chemistry Evaluation for Pele Through Recasting
======================================================

We use CEPTR to generate C++ mechanism files from `Cantera <https://cantera.org>`_ YAML chemistry files. CEPTR is a python package part of the PelePhysics source code.

.. _sec_ceptr_software:

Software requirements
---------------------

The CEPTR package uses `uv <https://docs.astral.sh/uv/getting-started/installation/>`_ to manage the Python dependencies. uv can typically be installed through system package managers (e.g. HomeBrew) or following the instructions in uv's documentation.

To install CEPTR and its dependencies::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ uv sync --locked

This installs the exact versions recorded in ``Support/ceptr/uv.lock``, so that every user generates the same mechanisms.

CEPTR supports Python 3.13 and 3.14. You do not need to install these yourself: if the system Python is older (e.g. 3.12 on Ubuntu 24.04), uv downloads and uses a managed Python 3.13 (set in ``Support/ceptr/.python-version``), leaving the system Python untouched. Only a ``pip`` install (see the note below) needs a 3.13 or 3.14 interpreter already present.

The plotting packages (matplotlib, networkx and SciPy) are an optional ``plot`` extra: they are only needed for
``qssa --visualize`` and for the plotting scripts elsewhere in PelePhysics (``plotPMF.py``, ``gen_hit_ic.py``).
To install them too::

  $ uv sync --locked --extra plot

``uv sync`` removes what it is not asked for, so a later ``uv sync --locked`` uninstalls them again: give
``--extra plot`` to every sync (or to ``uv run``) when you need them.

The committed mechanisms are exactly what CEPTR generates with these versions: the CI regenerates every mechanism in
``Mechanisms/list_mech`` and ``Mechanisms/list_qss_mech`` on Linux and fails if any differs. When the check fails, the
``regenerated-mechanisms`` artifact of the CI run holds the differences as a patch::

  $ gh run download <run-id> -n regenerated-mechanisms
  $ git apply mechanisms.patch

.. note:: CEPTR is a standard Python package, so it can also be installed with ``pip install -e .`` from ``Support/ceptr``. pip then picks any versions allowed by ``pyproject.toml``, which may generate mechanisms that differ from the committed ones.

Usage
-----

Generating YAML chemistry files: converting from CHEMKIN files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. _sec_convertCK:

We rely on Cantera's ``ck2yaml`` utility to convert CHEMKIN files to the Cantera YAML format::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ uv run ck2yaml --input=${PATH_TO_CHEMKIN_DIR}/mechanism.inp --thermo=${PATH_TO_CHEMKIN_DIR}/therm.dat --transport=${PATH_TO_CHEMKIN_DIR}/tran.dat --permissive

Of course, the files ``tran.dat`` and ``therm.dat`` are optional if already included in the ``mechanism.inp`` file.


Generating Pele-compatible mechanisms for a single chemistry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three ways to use CEPTR to generate C++ mechanism files for a given chemistry.

1. Using CEPTR directly:

   Executed from the ``${PELE_PHYSICS_HOME}/Support/ceptr`` directory, the most general usage of CEPTR is::

     $ uv run convert -f ${PELE_PHYSICS_HOME}/Mechanisms/${chemistry}/mechanism.yaml \
       --chemistry {chemistry-type} \
       --gas_name {gas-name} \
       --interface_name {interface-name} \
       --plog_pressure {value in Pa}

   The ``--chemistry`` (or ``-c``) argument allows users to convey if the ``${chemistry}`` of interest is either one of two valid options, namely, ``homogeneous`` or ``heterogeneous``. The default value for ``{chemistry-type}`` is ``homogeneous``.
   The ``--gas_name`` argument allows users to specify the names of the homogeneous phase to use, as several different ones can be defined in the corresponding ``mechanism.yaml`` file (under the ``phases:`` item). The default value for ``{gas-name}`` is ``gas``.
   Finally, the ``--interface_name`` arguments allow users to specify the name of the gas-solid interface, also prescribed in the corresponding ``mechanism.yaml`` file. The default value for ``{interface-name}`` is ``None``.

   Note that if ``--chemistry heterogeneous`` is specified, the user must necessarily specify a corresponding ``{interface-name}``.

   An example of directly using CEPTR for homogeneous mechanisms is::

     $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
     $ uv run convert -f ${PELE_PHYSICS_HOME}/Mechanisms/LiDryer/mechanism.yaml

   .. note:: CEPTR interpretations of heterogeneous mechanisms is currently a work in progress.

   .. note:: CEPTR offers only partial support for mechanisms with pressure-dependent Arrhenius (PLOG) reactions.
             These reactions are fixed to a constant pressure, specified with the ``--plog_pressure`` argument,
             in the CEPTR-generated source code. Mechanisms with PLOG reactions generated in this way are limited
             in validity to that specific pressure and would not be suitable for compressible flows with significant
             pressure variation. For that reason, source code for these mechanisms is saved in pressure-specific directories, for
             example ``${PELE_PHYSICS_HOME}/Mechanisms/POLIMI2020/1_000atm/``.


2. Using a helper script in the directory containing the ``mechanism.yaml`` file::

     $ ./convert.sh

   .. note:: The helper scripts run CEPTR with ``uv run --locked``, installing the dependencies on first use. Without uv,
             they fall back to ``python3`` in the current environment, which then needs CEPTR and its dependencies
             installed (e.g. with ``pip install -e ${PELE_PHYSICS_HOME}/Support/ceptr``).

3. Using a helper script in the ``Mechanisms`` directory::

     $ bash ${PELE_PHYSICS_HOME}/Mechanisms/converter.sh -f ./LiDryer/mechanism.yaml


Generating a reduced, QSS chemistry file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step consists in generating a QSS chemistry YAML file from a skeletal or a detailed YAML file. To do so, one needs: the mechanism YAML file ``skeletal.yaml``, as well as a list of non-QSS species, ``non_qssa_list.yaml``. The following command will generate a QSS YAML file, ``qssa.yaml``::

    $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
    $ uv run qssa -f ${PATH_TO_YAML}/skeletal.yaml -n ${PATH_TO_YAML}/non_qssa_list.yaml

The full list of options is::

  $ uv run qssa -h
  usage: qssa [-h] -f FNAME -n NQSSA [-m {0,1,2}] [-v]

  Mechanism converter

  optional arguments:
    -h, --help            show this help message and exit
    -f FNAME, --fname FNAME
                          Mechanism file
    -n NQSSA, --nqssa NQSSA
                          Non-QSSA species list
    -m {0,1,2}, --method {0,1,2}
                          QSSA method (default: 2)
    -v, --visualize       Visualize quadratic coupling and QSSA dependencies

For a detailed description of these options and a further information on the way QSS mechanism are treated in CEPTR the reader may consult :ref:`the QSS section <sec_qss>`.

To generate a QSS C++ mechanism from the ``.yaml`` file thus created, tailored to your needs, please refer to Tutorials :ref:`Generating NC12H26 QSS mechanism with analytical jacobian <sec_tutqss1>` and :ref:`Generating NC12H26 QSS mechanism without analytical jacobian <sec_tutqss2>`.

Batched generation of Pele-compatible mechanisms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   If you are using batched generation as outlined here, it will automatically use multiprocessing to generate the files in parallel using all CPUs detected on the machine. If you want to change that you can pass the optional argument ``-n NPCU``, wheren ``NCPU`` is an integer indicating the number of processes you want to use.


For non-reduced chemistries, CEPTR can take a file with a list of ``mechanism.yaml`` files to convert::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ uv run convert -l ${PELE_PHYSICS_HOME}/Mechanisms/list_mech

The list files take one mechanism per line; blank lines and lines starting with ``#`` are skipped. A mechanism with PLOG
reactions needs its pressure in Pa, given on its line as ``mechanism.yaml --plog=101325.0``; for the lines without it,
``--plog_pressure`` gives the pressure for the whole list.

For reduced chemistries, CEPTR can take a file with a list of ``qssa.yaml`` and ``qssa_input.toml`` to convert::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ uv run convert -lq ${PELE_PHYSICS_HOME}/Mechanisms/list_qss_mech

For generating ``qssa.yaml`` for reduced chemistries, CEPTR can take a file with a list of ``skeletal.yaml`` and ``non_qssa_list.yaml``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ uv run qssa -lq ${PELE_PHYSICS_HOME}/Mechanisms/list_qss_mech

To generate all mechanisms::

  $ uv run convert -l ${PELE_PHYSICS_HOME}/Mechanisms/list_mech
  $ uv run qssa -lq ${PELE_PHYSICS_HOME}/Mechanisms/list_qss_mech
  $ uv run convert -lq ${PELE_PHYSICS_HOME}/Mechanisms/list_qss_mech

Planck-mean absorption coefficients for radiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``--radiation``, CEPTR also writes the Planck-mean absorption coefficients of the radiating gases in the mechanism
into ``mechanism.H``::

  $ uv run convert -f ${PELE_PHYSICS_HOME}/Mechanisms/dodecane_lu/mechanism.yaml --radiation

``get_planck_mean_kp(T, kp)`` then gives, for each gas species ``k``, its Planck-mean absorption coefficient per unit
partial pressure ``kp[k]``, in 1/(cm bar), and 0 for species without data. The absorption coefficient of the mixture is
:math:`\kappa_P = \sum_k X_k P k_{p,k}(T)`. ``NUM_PLANCK_MEAN_SPECIES`` is the number of species with data. Mechanisms
converted without ``--radiation`` are unchanged.

The data are fits of :math:`\ln k_p` against :math:`x = \ln(T / 1000\,\mathrm{K})`, one coefficient file per species in
``Support/ceptr/ceptr/radiation_data``. A species of the mechanism uses a file when it has the same name or one of its
aliases (in any case), and the same elemental composition. As with the NASA polynomials of the thermodynamics, nothing
is computed from spectroscopic data during a conversion. Each fit has

* the fewest polynomial ranges (up to three) of degree up to 6 between :math:`T_\mathrm{low}` and
  :math:`T_\mathrm{max}`, continuous at the breaks, within 0.5% of the data; :math:`T_\mathrm{low}` is 250 K, or 300 K
  for C2H4 and NH3, whose data are given every 100 K;
* :math:`k_p` held at its :math:`T_\mathrm{low}` value below :math:`T_\mathrm{low}`;
* above :math:`T_\mathrm{max}`, the highest temperature of the line list (at most 5000 K), the power law
  :math:`k_p(T) = k_p(T_\mathrm{max}) (T / T_\mathrm{max})^s`, with :math:`s` the slope of
  :math:`\ln k_p` against :math:`\ln T` at :math:`T_\mathrm{max}` (0 if the fit still rises there).

The coefficient files, and the comments of the generated code, give the source, the fit ranges, the exponent :math:`s`
and the largest fit error of each species.

The coefficient files are written by the ``planck-mean`` utility, which downloads ExoMol data [ExoMol2016]_ and fits
them; it is only needed to add a species or to update the data::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ uv run planck-mean            # the species of the table below
  $ uv run planck-mean CO2 H2O    # some of them
  $ uv run planck-mean SO2        # another ExoMol molecule

The species of the table are those of the PelePhysics mechanisms with ExoMol data, with the line lists chosen in
``DATASETS`` of ``ceptr/planck_mean.py`` (where SO2, H2S and OCS are also set up). Any other ExoMol molecule can be
given by its formula: ``planck-mean`` then takes the cooling function (or super-lines) of the recommended line list of
its main isotopologue, or of another of its line lists, and the references listed on its ExoMol page, which should be
checked. Left out of the table are N2, O2 and H2, whose :math:`k_p`, from weak electric-quadrupole and electronic
transitions, is 7 to over 100 orders of magnitude below that of H2O, and CH2O and H2O2, which ExoMol only has as full
line lists (120-190 GB).

In local thermodynamic equilibrium the Planck-mean absorption coefficient is the total spontaneous emission of the gas
over that of a black body (Kirchhoff's law). ExoMol tabulates the former as the cooling function :math:`W(T)`, the
power emitted per molecule and steradian, computed with ExoCross [ExoCross2018]_ from line lists of up to billions of
lines, so that

.. math::

   k_p(T) = \frac{N \pi W(T)}{\sigma T^4}, \qquad N = \frac{P}{k_B T}, \quad P = 1\ \mathrm{bar}.

For NH3 and C2H4, which have no cooling function, :math:`k_p = N \sum S\, E_b(\nu, T) / (\sigma T^4)` is computed from
the ExoMol super-lines (line intensities :math:`S` summed in narrow wavenumber bins) at each temperature they are given
for; these files are read as they download and are not stored. Where ExoMol has no natural-abundance data, the main
isotopologue stands for the gas (for CO, the two agree to 0.2%).

.. list-table:: Data of the Planck-mean absorption coefficients
   :header-rows: 1

   * - Species
     - ExoMol line list
     - Data
     - Line list up to
     - References
   * - H2O
     - POKAZATEL (natural abundance)
     - cooling function
     - 6000 K
     - [POKAZATEL]_
   * - CO2
     - UCL-4000 (12C-16O2)
     - cooling function
     - 4000 K
     - [UCL4000]_
   * - CO
     - Li2015 (natural abundance)
     - cooling function
     - 9000 K
     - [Li2015]_
   * - CH4
     - MM (12C-1H4)
     - cooling function
     - 2000 K
     - [MM]_
   * - C2H4
     - MaYTY (12C2-1H4)
     - super-lines, 200-700 K
     - 700 K
     - [MaYTY]_
   * - C2H2
     - aCeTY (12C2-1H2)
     - cooling function
     - 2200 K
     - [aCeTY]_
   * - NH3
     - CoYuTe (14N-1H3)
     - super-lines, 200-2000 K
     - 2000 K
     - [CoYuTe]_ [CoYuTeM]_
   * - NO
     - XABC (14N-16O)
     - cooling function
     - 3500 K
     - [XABC]_
   * - N2O
     - TYM (natural abundance)
     - cooling function
     - 2000 K
     - [TYM]_
   * - OH
     - MYTHOS (16O-1H)
     - cooling function
     - 10000 K
     - [MYTHOS]_
   * - HCN
     - Harris (1H-12C-14N)
     - cooling function
     - 4000 K
     - [Harris]_ [HarrisB]_
   * - CH3
     - AYYJ (12C-1H3)
     - cooling function
     - 2000 K
     - [AYYJ]_
   * - CH
     - MoLLIST (12C-1H)
     - cooling function
     - 8000 K
     - [MoLLIST]_ [MoLLISTCH]_
   * - NH
     - kNigHt (14N-1H)
     - cooling function
     - 10000 K
     - [kNigHt]_
   * - CN
     - KTPSYT (12C-14N)
     - cooling function
     - 10000 K
     - [KTPSYT]_
   * - HF
     - Coxon-Hajig (1H-19F)
     - cooling function
     - 5000 K
     - [CoxonHajig]_ [CoxonHajigL]_ [CoxonHajigS]_
   * - CH3F
     - OYKYT (12C-1H3-19F)
     - cooling function
     - 500 K
     - [OYKYT]_
   * - HNO3 (also HONO2)
     - AIJS (1H-14N-16O3)
     - cooling function
     - 500 K
     - [AIJS]_

.. note::

   **Licence and attribution.** The ExoMol data are released under the Creative Commons Attribution-ShareAlike 4.0
   International licence (`CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0/>`_). The coefficient files in
   ``Support/ceptr/ceptr/radiation_data``, and the fits written into ``mechanism.H`` with ``--radiation``, are derived
   from them and are released under CC BY-SA 4.0 too; the rest of PelePhysics keeps its own licence. Work using these
   coefficients should cite the ExoMol database [ExoMol2016]_ and the line lists of the species used, listed above, in
   the coefficient files and in the generated code; those from cooling functions also ExoCross [ExoCross2018]_.

.. [ExoMol2016] J. Tennyson et al., The ExoMol database: molecular line lists for exoplanet and other hot atmospheres,
   J. Mol. Spectrosc. 327, 73-94 (2016), doi:10.1016/j.jms.2016.05.002
.. [ExoCross2018] S. N. Yurchenko, A. F. Al-Refaie, J. Tennyson, ExoCross: a general program for generating spectra
   from molecular line lists, Astron. Astrophys. 614, A131 (2018), doi:10.1051/0004-6361/201732531
.. [POKAZATEL] O. L. Polyansky, A. A. Kyuberis, N. F. Zobov, J. Tennyson, S. N. Yurchenko, L. Lodi, ExoMol molecular
   line lists XXX: a complete high-accuracy line list for water, Mon. Not. R. Astron. Soc. 480, 2597-2608 (2018),
   doi:10.1093/mnras/sty1877
.. [UCL4000] S. N. Yurchenko, T. M. Mellor, R. S. Freedman, J. Tennyson, ExoMol line lists - XXXIX. Ro-vibrational
   molecular line list for CO2, Mon. Not. R. Astron. Soc. 496, 5282-5291 (2020), doi:10.1093/mnras/staa1874
.. [Li2015] G. Li, I. E. Gordon, L. S. Rothman, Y. Tan, S.-M. Hu, S. Kassi, A. Campargue, E. S. Medvedev,
   Rovibrational line lists for nine isotopologues of the CO molecule in the X 1Sigma+ ground electronic state,
   Astrophys. J. Suppl. Ser. 216, 15 (2015), doi:10.1088/0067-0049/216/1/15
.. [MM] S. N. Yurchenko, A. Owens, K. Kefala, J. Tennyson, ExoMol line lists - LVII: High accuracy ro-vibrational line
   list for methane (CH4), Mon. Not. R. Astron. Soc. 528, 3719-3729 (2024), doi:10.1093/mnras/stae148
.. [MaYTY] B. P. Mant, A. Yachmenev, J. Tennyson, S. N. Yurchenko, ExoMol molecular line lists - XXVII: spectra of
   C2H4, Mon. Not. R. Astron. Soc. 478, 3220-3232 (2018), doi:10.1093/mnras/sty1239
.. [aCeTY] K. L. Chubb, J. Tennyson, S. N. Yurchenko, ExoMol molecular line lists - XXXVII: spectra of acetylene,
   Mon. Not. R. Astron. Soc. 493, 1531-1545 (2020), doi:10.1093/mnras/staa229
.. [CoYuTe] P. A. Coles, S. N. Yurchenko, J. Tennyson, ExoMol molecular line lists - XXXV. A rotation-vibration line
   list for hot ammonia, Mon. Not. R. Astron. Soc. 490, 4638-4647 (2019), doi:10.1093/mnras/stz2778
.. [CoYuTeM] A. R. Al Derzi, T. Furtenbacher, J. Tennyson, S. N. Yurchenko, A. G. Csaszar, MARVEL analysis of the
   measured high-resolution spectra of 14NH3, J. Quant. Spectrosc. Radiat. Transf. 161, 117-130 (2015),
   doi:10.1016/j.jqsrt.2015.03.034
.. [XABC] Q. Qu, S. N. Yurchenko, J. Tennyson, ExoMol molecular line lists - XLII. Rovibronic molecular line list for
   the low-lying states of NO, Mon. Not. R. Astron. Soc. 504, 5768-5777 (2021), doi:10.1093/mnras/stab1154
.. [TYM] S. N. Yurchenko, T. M. Mellor, J. Tennyson, ExoMol line lists - LIX. High-temperature line list for N2O,
   Mon. Not. R. Astron. Soc. 534, 1364-1375 (2024), doi:10.1093/mnras/stae2201
.. [MYTHOS] G. B. Mitev, C. A. Bowesman, J. Zhang, S. N. Yurchenko, J. Tennyson, ExoMol line lists - LXI. A trihybrid
   line list for rovibronic transitions of the hydroxyl radical (OH), Mon. Not. R. Astron. Soc. 536, 3401-3420
   (2025), doi:10.1093/mnras/stae2803
.. [Harris] G. J. Harris, J. Tennyson, B. M. Kaminsky, Ya. V. Pavlenko, H. R. A. Jones, Improved HCN/HNC linelist,
   model atmospheres and synthetic spectra for WZ Cas, Mon. Not. R. Astron. Soc. 367, 400-406 (2006),
   doi:10.1111/j.1365-2966.2005.09960.x
.. [HarrisB] R. J. Barber, J. K. Strange, C. Hill, O. L. Polyansky, G. Ch. Mellau, S. N. Yurchenko, J. Tennyson,
   ExoMol line lists - III. An improved hot rotation-vibration line list for HCN and HNC, Mon. Not. R. Astron. Soc.
   437, 1828-1835 (2014), doi:10.1093/mnras/stt2011
.. [AYYJ] A. Y. Adam, A. Yachmenev, S. N. Yurchenko, P. Jensen, Variationally computed IR line list for the methyl
   radical CH3, J. Phys. Chem. A 123, 4755-4763 (2019), doi:10.1021/acs.jpca.9b02919
.. [MoLLIST] P. F. Bernath, MoLLIST: Molecular line lists, intensities and spectra, J. Quant. Spectrosc. Radiat. Transf.
   240, 106687 (2020), doi:10.1016/j.jqsrt.2019.106687
.. [MoLLISTCH] T. Masseron et al., CH in stellar atmospheres: an extensive linelist, Astron. Astrophys. 571, A47 (2014),
   doi:10.1051/0004-6361/201423956
.. [kNigHt] A. N. Perri, L. K. McKemmish, Full spectroscopic model and trihybrid experimental-perturbative-variational
   line list for NH, Mon. Not. R. Astron. Soc. 531, 3023-3033 (2024), doi:10.1093/mnras/stae1340
.. [KTPSYT] S. V. Kozlov, V. A. Terashkevich, E. A. Pazyuk, A. V. Stolyarov, S. N. Yurchenko, J. Tennyson, An upgraded
   line list for radiative transitions between the B 2Sigma+, A 2Pi, and X 2Sigma+ states of the CN radical,
   Astrophys. J. Suppl. Ser. 275, 29 (2024), doi:10.3847/1538-4365/ad8b40
.. [CoxonHajig] J. A. Coxon, P. G. Hajigeorgiou, Improved direct potential fit analyses for the ground electronic states
   of the hydrogen halides: HF/DF/TF, HCl/DCl/TCl, HBr/DBr/TBr and HI/DI/TI, J. Quant. Spectrosc. Radiat. Transf. 151,
   133-154 (2015), doi:10.1016/j.jqsrt.2014.08.028
.. [CoxonHajigL] G. Li, I. E. Gordon, R. J. Le Roy, P. G. Hajigeorgiou, J. A. Coxon, P. F. Bernath, L. S. Rothman,
   Reference spectroscopic data for hydrogen halides. Part I: Construction and validation of the ro-vibrational dipole
   moment functions, J. Quant. Spectrosc. Radiat. Transf. 121, 78-90 (2013), doi:10.1016/j.jqsrt.2013.02.005
.. [CoxonHajigS] W. Somogyi, S. N. Yurchenko, A. Yachmenev, Calculation of electric quadrupole linestrengths for
   diatomic molecules: Application to the H2, CO, HF, and O2 molecules, J. Chem. Phys. 155, 214303 (2021),
   doi:10.1063/5.0063256
.. [OYKYT] A. Owens, A. Yachmenev, J. Kupper, S. N. Yurchenko, W. Thiel, The rotation-vibration spectrum of methyl
   fluoride from first principles, Phys. Chem. Chem. Phys. 21, 3496-3505 (2019), doi:10.1039/C8CP01721B
.. [AIJS] A. I. Pavlyuchko, S. N. Yurchenko, J. Tennyson, ExoMol molecular line lists - XI. The spectrum of nitric
   acid, Mon. Not. R. Astron. Soc. 452, 1702-1706 (2015), doi:10.1093/mnras/stv1376
