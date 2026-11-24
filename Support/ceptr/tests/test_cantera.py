"""Pointwise comparison of CEPTR-generated mechanisms with Cantera.

Each mechanism is converted with CEPTR into a temporary directory, compiled
with the evaluator of Testing/Exec/CeptrCantera (AMReX headers only, nothing
to build first) and evaluated at reproducible random states. Thermodynamics,
rates of progress, production rates and transport are compared with Cantera,
in CGS units, quantity by quantity.

These tests need a C++ compiler, make and the AMReX submodule, and take a
minute or two, so they are not run by default: run them with

    uv run pytest -m cantera -s -rx

-s prints the measured errors, -rx the known discrepancies. CEPTR_TEST_JOBS
sets how many mechanisms are converted and compiled at once (default: half
the number of CPUs, each compiling two files at once).
"""

import concurrent.futures
import dataclasses
import os
import pathlib
import shutil
import subprocess
import sys

import cantera as ct
import numpy as np
import pytest

pytestmark = pytest.mark.cantera

REPO = pathlib.Path(__file__).resolve().parents[3]
EVALUATOR = REPO / "Testing" / "Exec" / "CeptrCantera"

# A mechanism with the SRI falloff (3 and 5 parameters), which none of the
# committed mechanisms use, with Troe and Lindemann falloffs and third bodies
# around it. Rates from GRI-Mech 3.0; the SRI parameters are illustrative.
SRI_YAML = """
units: {length: cm, quantity: mol, activation-energy: cal/mol}
phases:
- name: gas
  thermo: ideal-gas
  elements: [O, H, C, N, Ar]
  species:
  - gri30.yaml/species: [H2, H, O, O2, OH, H2O, HO2, H2O2, CH3, CH4, CO, CO2,
      N2, AR]
  kinetics: gas
  transport: mixture-averaged
  reactions: all
reactions:
- equation: CH3 + H (+M) <=> CH4 (+M)
  type: falloff
  low-P-rate-constant: {A: 2.62e+33, b: -4.76, Ea: 2440.0}
  high-P-rate-constant: {A: 1.39e+16, b: -0.534, Ea: 536.0}
  SRI: {A: 0.783, B: 74.0, C: 2941.0}
  efficiencies: {H2: 2.0, H2O: 6.0, CH4: 2.0, CO: 1.5, CO2: 2.0, AR: 0.7}
- equation: H + O2 (+M) <=> HO2 (+M)
  type: falloff
  low-P-rate-constant: {A: 2.8e+18, b: -0.86, Ea: 0.0}
  high-P-rate-constant: {A: 4.65e+12, b: 0.44, Ea: 0.0}
  SRI: {A: 0.45, B: 797.0, C: 979.0, D: 1.2, E: 0.1}
  efficiencies: {O2: 0.0, H2O: 0.0, CO: 0.75, CO2: 1.5, N2: 0.0, AR: 0.0}
- equation: 2 OH (+M) <=> H2O2 (+M)
  type: falloff
  low-P-rate-constant: {A: 2.3e+18, b: -0.9, Ea: -1700.0}
  high-P-rate-constant: {A: 7.4e+13, b: -0.37, Ea: 0.0}
  Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
  efficiencies: {H2: 2.0, H2O: 6.0, CH4: 2.0, CO: 1.5, CO2: 2.0, AR: 0.7}
- equation: CO + O (+M) <=> CO2 (+M)
  type: falloff
  low-P-rate-constant: {A: 6.02e+14, b: 0.0, Ea: 3000.0}
  high-P-rate-constant: {A: 1.8e+10, b: 0.0, Ea: 2385.0}
  efficiencies: {H2: 2.0, O2: 6.0, H2O: 6.0, CH4: 2.0, CO: 1.5, CO2: 3.5,
    AR: 0.5}
- equation: H + O2 + AR <=> HO2 + AR
  rate-constant: {A: 7.0e+17, b: -0.8, Ea: 0.0}
- equation: H + OH + M <=> H2O + M
  type: three-body
  rate-constant: {A: 2.2e+22, b: -2.0, Ea: 0.0}
  efficiencies: {H2: 0.73, H2O: 3.65, CH4: 2.0, AR: 0.38}
- equation: H2 + O <=> H + OH
  rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
- equation: CH4 + OH => CH3 + H2O
  rate-constant: {A: 1.0e+08, b: 1.6, Ea: 3120.0}
"""


# Falloff reactions with a specific collider, as in C3MechLite
COLLIDER_YAML = """
units: {length: cm, quantity: mol, activation-energy: cal/mol}
phases:
- name: gas
  thermo: ideal-gas
  elements: [O, H, N, Ar]
  species:
  - gri30.yaml/species: [H2, H, O, O2, OH, H2O, HO2, N2, AR]
  kinetics: gas
  transport: mixture-averaged
  reactions: all
reactions:
- equation: H + O2 (+N2) <=> HO2 (+N2)
  type: falloff
  low-P-rate-constant: {A: 1.91e+20, b: -1.5568, Ea: 253.86}
  high-P-rate-constant: {A: 4.66e+12, b: 0.44, Ea: 0.0}
  Troe: {A: 0.5, T3: 1.0, T1: 1.0e+10, T2: 1.0e+30}
- equation: H + O2 (+AR) <=> HO2 (+AR)
  type: falloff
  low-P-rate-constant: {A: 1.216e+19, b: -1.23, Ea: 0.0}
  high-P-rate-constant: {A: 4.66e+12, b: 0.44, Ea: 0.0}
- equation: H2 + O <=> H + OH
  rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
"""

# A PLOG reaction with two Arrhenius expressions per pressure, which add up,
# from C3MechLite
PLOG_SUM_YAML = """
units: {length: cm, quantity: mol, activation-energy: cal/mol}
phases:
- name: gas
  thermo: ideal-gas
  elements: [O, H, C, N]
  species:
  - gri30.yaml/species: [H2, H, O, O2, OH, H2O, HO2, CO, HCO, CH2O, C2H3, N2]
  kinetics: gas
  transport: mixture-averaged
  reactions: all
reactions:
- equation: C2H3 + O2 <=> CH2O + HCO
  type: pressure-dependent-Arrhenius
  rate-constants:
  - {P: 0.316 atm, A: 2.17e+36, b: -7.57, Ea: 1.249e+04}
  - {P: 0.316 atm, A: 5.34e+15, b: -1.29, Ea: 520.6}
  - {P: 1.0 atm, A: 3.03e+35, b: -7.32, Ea: 1.182e+04}
  - {P: 1.0 atm, A: 6.76e+15, b: -1.31, Ea: 645.7}
  - {P: 3.16 atm, A: 1.59e+36, b: -7.47, Ea: 1.246e+04}
  - {P: 3.16 atm, A: 1.05e+16, b: -1.36, Ea: 1066.0}
- equation: H2 + O <=> H + OH
  rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
"""


def specific_collider_falloff(reaction):
    """Whether a falloff reaction has a specific collider, as (+N2)."""
    return (
        reaction.rate.type == "falloff"
        and reaction.third_body.default_efficiency == 0.0
    )


def summed_plog(reaction):
    """Whether a PLOG reaction has several Arrhenius expressions per pressure."""
    if reaction.rate.type != "pressure-dependent-Arrhenius":
        return False
    pressures = [p for p, _ in reaction.rate.rates]
    return len(set(pressures)) < len(pressures)


@dataclasses.dataclass
class Case:
    """A mechanism to compare, and what it covers."""

    covers: str
    # In Mechanisms/, or the mechanism itself
    source: str
    # The PLOG pressure, for CEPTR, which evaluates PLOG at one pressure
    plog: float | None = None
    # Which reactions have a known discrepancy, each checked by a case of its
    # own: their rates of progress are left out of the comparison, and their
    # contribution out of the production rates
    known: tuple = ()
    # Discrepancies CEPTR still has in whole quantities: {quantity: (the
    # error they are known to give, reason)}. The test fails if one goes away.
    known_quantities: dict = dataclasses.field(default_factory=dict)
    # A known discrepancy: (reason, the error CEPTR fails with). Without the
    # error, CEPTR converts the mechanism but the forward rates of its known
    # reactions must still differ.
    xfail: tuple | None = None


# Together they cover the reaction types CEPTR supports: elementary,
# three-body (with efficiencies and with a specific collider), falloff (Troe,
# Lindemann, SRI), PLOG (with and without a third body, at and between the
# tabulated pressures), reversible and irreversible reactions (ck2yaml writes
# explicit reverse parameters as pairs of irreversible reactions),
# duplicates, negative pre-exponential factors, non-unity and non-reactant
# orders and fractional stoichiometric coefficients. CEPTR has no Chebyshev
# nor chemically-activated reactions. QSS mechanisms are left out: QSS is an
# approximation, which Cantera does not make.
CASES = {
    "grimech30": Case(
        "elementary, three-body, specific three-body collider, Troe,"
        " Lindemann, irreversible, duplicates",
        "grimech30/mechanism.yaml",
    ),
    "C3MechLite": Case(
        "PLOG at a tabulated pressure, three-body PLOG, negative A,"
        " duplicates, Troe, Lindemann",
        "C3MechLite/mechanism.yaml",
        plog=ct.one_atm,
        known=(specific_collider_falloff, summed_plog),
    ),
    "POLIMI2020": Case(
        "PLOG between tabulated pressures, Troe, Lindemann",
        "POLIMI2020/mechanism.yaml",
        plog=2.0 * ct.one_atm,
    ),
    "posf10325": Case(
        "fractional stoichiometric coefficients",
        "posf10325/mechanism.yaml",
        known_quantities={
            "net_production_rates": (
                1e-5,
                "productionRate writes the fractional stoichiometric"
                " coefficients with 6 decimals (0.0871338 iC4H8 as 0.087134)",
            ),
        },
    ),
    "chem-CH4-2step": Case(
        "non-unity and non-reactant orders", "chem-CH4-2step/mechanism.yaml"
    ),
    "JL4": Case("orders, reversible global reactions", "JL4/mechanism.yaml"),
    "falloff-collider": Case(
        "falloff with a specific collider (Troe, Lindemann)",
        COLLIDER_YAML,
        known=(specific_collider_falloff,),
        xfail=(
            "CEPTR ignores the specific collider of falloff reactions such as"
            " H + O2 (+N2) <=> HO2 (+N2): their third-body concentration is"
            " the whole mixture, not [N2] (utilities.enhancement_d ignores the"
            " default efficiency of 0)",
            None,
        ),
    ),
    "plog-sum": Case(
        "PLOG with several expressions per pressure",
        PLOG_SUM_YAML,
        plog=ct.one_atm,
        known=(summed_plog,),
        xfail=(
            "CEPTR uses one of the Arrhenius expressions given at a pressure"
            " of a PLOG reaction instead of their sum (utilities.evaluate_plog)",
            None,
        ),
    ),
    "SRI": Case(
        "SRI falloff (3 and 5 parameters), Troe, Lindemann, third bodies",
        SRI_YAML,
        xfail=(
            "CEPTR looks for the falloff type 'Sri', Cantera names it 'SRI':"
            " SRI reactions are rejected as unrecognized",
            "Unrecognized reaction rate type falloff, SRI",
        ),
    ),
}

N_RANDOM = 12
N_EQUILIBRIUM = 8
SEED = 20261124

# Maximum relative error per quantity, set from the errors measured over
# 210 states per mechanism, with a margin. Thermodynamics and kinetics are
# evaluated from the same data in both codes, so their errors are round-off:
# up to 1e-14 for the thermodynamics, 3e-13 for the forward rates, and 1.2e-12
# for the reverse rates, through exp(-g/RT) with g/RT of a few hundred (N2 + M
# <=> 2 N + M at 320 K). The transport properties are cubic fits in log(T) of
# CEPTR's kinetic theory, which differs from Cantera's: up to 2.4e-2 for the
# viscosities and 1.6e-2 for the binary diffusion coefficients (helium near
# 3000 K, where the reduced temperature is beyond the collision integral
# tables), 4.4e-2 for the species conductivities (the atoms, whose electronic
# heat capacity Cantera counts as internal energy, near 300 K).
TOLERANCES = dict.fromkeys(
    (
        "cp_R",
        "h_RT",
        "s_R",
        "g_RT",
        "cp_mass_k",
        "cv_mass_k",
        "h_mass_k",
        "u_mass_k",
        "s_mass_k",
        "density",
        "mean_molecular_weight",
        "cp_mass",
        "cv_mass",
        "enthalpy_mass",
        "int_energy_mass",
        "entropy_mass",
        "forward_rates_of_progress",
    ),
    1e-12,
) | {
    "reverse_rates_of_progress": 5e-12,
    "net_production_rates": 5e-12,
    "species_viscosities": 4e-2,
    "species_conductivities": 8e-2,
    "binary_diff_coeffs_1atm": 3e-2,
    "viscosity": 1e-2,
    "thermal_conductivity": 5e-2,
    "rho_mix_diff_coeffs": 2e-2,
}

# Rates more than this many times smaller than the largest at a state are not
# compared: they can underflow differently and do not matter
RATE_FLOOR = 1e-30


def sample_states(gas, rng, pressure):
    """Return random (T, P, X) states, a third of them at equilibrium.

    The temperatures stay in the range where every species' thermodynamic
    fit is valid. The pressure is fixed for PLOG mechanisms, which CEPTR
    evaluates at one pressure.
    """
    t_min = max(300.0, max(s.thermo.min_temp for s in gas.species()))
    t_max = min(3000.0, min(s.thermo.max_temp for s in gas.species()))
    states = []
    for n in range(N_RANDOM + N_EQUILIBRIUM):
        temp = rng.uniform(t_min, t_max)
        pres = pressure or ct.one_atm * 10.0 ** rng.uniform(-1.0, 2.0)
        # Every species present, over eight orders of magnitude
        x = 10.0 ** rng.uniform(-8.0, 0.0, gas.n_species)
        gas.TPX = temp, pres, x
        if n >= N_RANDOM:
            gas.equilibrate("TP")
        # Not gas.P: equilibrate changes it in the last digits, which moves
        # PLOG off the pressure CEPTR evaluates it at
        states.append((temp, pres, gas.X))
    return states


def cantera_values(gas, temp, pres, x):
    """Return Cantera's values of the evaluator's quantities, in CGS."""
    gas.TPX = temp, pres, x
    rgas = ct.gas_constant / gas.molecular_weights  # J/kg/K
    cgs_energy = 1e4  # J/kg -> erg/g
    # Cantera's standard state is at the pressure of the gas, CEPTR's at the
    # reference pressure (1 atm)
    s_r = gas.standard_entropies_R + np.log(pres / gas.reference_pressure)
    ref = {
        "cp_R": gas.standard_cp_R,
        "h_RT": gas.standard_enthalpies_RT,
        "s_R": s_r,
        "g_RT": gas.standard_enthalpies_RT - s_r,
        "cp_mass_k": gas.standard_cp_R * rgas * cgs_energy,
        "cv_mass_k": (gas.standard_cp_R - 1.0) * rgas * cgs_energy,
        "h_mass_k": gas.standard_enthalpies_RT * rgas * temp * cgs_energy,
        "u_mass_k": (gas.standard_enthalpies_RT - 1.0) * rgas * temp * cgs_energy,
        "s_mass_k": s_r * rgas * cgs_energy,
        "density": gas.density * 1e-3,
        "mean_molecular_weight": gas.mean_molecular_weight,
        "cp_mass": gas.cp_mass * cgs_energy,
        "cv_mass": gas.cv_mass * cgs_energy,
        "enthalpy_mass": gas.enthalpy_mass * cgs_energy,
        "int_energy_mass": gas.int_energy_mass * cgs_energy,
        "entropy_mass": gas.entropy_mass * cgs_energy,
    }
    # Error scales: the enthalpies are compared relative to R T (they cross
    # zero), the net production rates relative to the gross rates
    scale = {
        "h_RT": np.maximum(np.abs(ref["h_RT"]), 1.0),
        "g_RT": np.maximum(np.abs(ref["g_RT"]), 1.0),
        "h_mass_k": np.maximum(np.abs(ref["h_mass_k"]), rgas * temp * cgs_energy),
        "u_mass_k": np.maximum(np.abs(ref["u_mass_k"]), rgas * temp * cgs_energy),
        "enthalpy_mass": max(abs(ref["enthalpy_mass"]), ref["cp_mass"] * temp),
        "int_energy_mass": max(abs(ref["int_energy_mass"]), ref["cv_mass"] * temp),
    }
    if gas.n_reactions > 0:
        kmol_m3 = 1e-3  # kmol/m^3/s -> mol/cm^3/s
        for name in ("forward_rates_of_progress", "reverse_rates_of_progress"):
            ref[name] = getattr(gas, name) * kmol_m3
            scale[name] = np.maximum(
                np.abs(ref[name]), RATE_FLOOR * np.abs(ref[name]).max()
            )
        ref["net_production_rates"] = gas.net_production_rates * kmol_m3
        gross = (gas.creation_rates + gas.destruction_rates) * kmol_m3
        scale["net_production_rates"] = np.maximum(gross, RATE_FLOOR * gross.max())

    # Transport, mixture-averaged. Pure species conductivities are those of
    # the pure gases. Simple transport's mixture rules for the viscosity and
    # conductivity are not Cantera's: they are applied to Cantera's species
    # values.
    mu_k = gas.species_viscosities * 10.0  # Pa s -> g/cm/s
    ref["species_viscosities"] = mu_k
    ref["binary_diff_coeffs_1atm"] = (
        gas.binary_diff_coeffs * pres / ct.one_atm * 1e4  # m^2/s -> cm^2/s
    ).ravel()
    ref["rho_mix_diff_coeffs"] = (
        gas.mix_diff_coeffs
        * gas.density
        * gas.molecular_weights
        / gas.mean_molecular_weight
        * 10.0  # kg/m/s -> g/cm/s
    )
    ref["viscosity"] = np.sum(x * mu_k**6) ** (1.0 / 6.0)
    lam_k = np.empty(gas.n_species)
    for k in range(gas.n_species):
        gas.TPX = temp, pres, {gas.species_name(k): 1.0}
        lam_k[k] = gas.thermal_conductivity * 1e5  # W/m/K -> erg/cm/s/K
    ref["species_conductivities"] = lam_k
    ref["thermal_conductivity"] = np.sum(x * lam_k**0.25) ** 4
    return ref, scale


def parse_output(text):
    """Return the evaluator's results, a dict of arrays per state."""
    states = []
    for line in text.splitlines():
        name, *values = line.split()
        if name == "state":
            states.append({})
        else:
            states[-1][name] = np.array([float(v) for v in values])
    return states


def run(cmd, **kwargs):
    """Run a command, with its output in the error if it fails."""
    proc = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{' '.join(map(str, cmd))} failed:\n{proc.stdout}\n{proc.stderr}"
        )
    return proc.stdout


def evaluate(name, case, workdir):
    """Convert, build and evaluate a mechanism, and compare it with Cantera.

    Returns the largest relative error of each quantity over the states (an
    Error), and that of the forward rates of the known reactions.
    """
    mechdir = workdir / name
    mechdir.mkdir()
    yaml = mechdir / "mechanism.yaml"
    if "\n" in case.source:
        yaml.write_text(case.source)
    else:
        shutil.copy(REPO / "Mechanisms" / case.source, yaml)
    # CEPTR copies it into the folder of a PLOG pressure
    (mechdir / "Make.package").touch()
    cmd = [sys.executable, "-m", "ceptr", "-f", yaml]
    if case.plog:
        cmd += ["-p", str(case.plog)]
        run(cmd)
        (mechdir,) = [d for d in mechdir.iterdir() if d.name.endswith("atm")]
    else:
        run(cmd)
    run(
        [
            "make",
            "-j2",
            "-C",
            EVALUATOR,
            f"MECH_DIR={mechdir}",
            f"PYTHON={sys.executable}",
        ]
    )

    gas = ct.Solution(yaml)
    rng = np.random.default_rng([SEED, sum(name.encode())])
    states = sample_states(gas, rng, case.plog)
    states_file = mechdir / "states.txt"
    with open(states_file, "w") as f:
        f.write(f"{gas.n_species} {' '.join(gas.species_names)}\n")
        for temp, pres, x in states:
            values = [temp, pres * 10.0, *x]  # Pa -> dyn/cm^2
            f.write(" ".join(f"{v:.17g}" for v in values) + "\n")
    results = parse_output(run([mechdir / "evaluate.ex", states_file]))
    assert len(results) == len(states)

    known = np.array(
        [any(k(r) for k in case.known) for r in gas.reactions()], dtype=bool
    )
    nu = gas.product_stoich_coeffs - gas.reactant_stoich_coeffs
    fwd, rev = "forward_rates_of_progress", "reverse_rates_of_progress"
    errors, known_errors = {}, {}

    def record(errors, quantity, err, n, temp, pres):
        # A NaN would never compare as the largest: count it as infinite
        err = np.where(np.isfinite(err), err, np.inf)
        i = int(np.argmax(err))
        if quantity not in errors or err[i] > errors[quantity].error:
            errors[quantity] = Error(
                float(err[i]), n, temp, pres, where(gas, quantity, i, err.size)
            )

    for n, ((temp, pres, x), result) in enumerate(zip(states, results, strict=True)):
        ref, scale = cantera_values(gas, temp, pres, x)
        assert set(ref) == set(result), set(ref) ^ set(result)
        if known.any():
            # The production rates without the known reactions, on both sides
            for values in (ref, result):
                q = values[fwd] - values[rev]
                values["net_production_rates"] = (
                    values["net_production_rates"] - nu[:, known] @ q[known]
                )
        for quantity, expected in ref.items():
            expected = np.atleast_1d(expected)
            denom = np.atleast_1d(scale.get(quantity, np.abs(expected)))
            diff = np.abs(result[quantity] - expected)
            with np.errstate(divide="ignore", invalid="ignore"):
                err = np.where(diff == 0.0, 0.0, diff / denom)
            if quantity in (fwd, rev) and known.any():
                # The forward rates show whether the discrepancy is still
                # there (an irreversible known reaction has no reverse rate)
                if quantity == fwd:
                    record(
                        known_errors, quantity, np.where(known, err, 0.0), n, temp, pres
                    )
                err = np.where(known, 0.0, err)
            record(errors, quantity, err, n, temp, pres)
    return errors, known_errors


@dataclasses.dataclass
class Error:
    """The largest error of a quantity, and where it is."""

    error: float
    state: int
    temp: float
    pres: float
    location: str

    def __str__(self):
        """Describe the error and where it is."""
        return (
            f"{self.error:9.2e}  state {self.state:2d} T={self.temp:4.0f} K"
            f" P={self.pres / ct.one_atm:<5.3g} atm  {self.location}"
        )


def where(gas, quantity, index, size):
    """Name the species, species pair or reaction of an error."""
    if "rates_of_progress" in quantity:
        return gas.reaction(index).equation
    if quantity == "binary_diff_coeffs_1atm":
        i, j = divmod(index, gas.n_species)
        return f"{gas.species_name(i)}-{gas.species_name(j)}"
    return gas.species_name(index) if size == gas.n_species else ""


@pytest.fixture(scope="module")
def measured(tmp_path_factory):
    """Convert, build and evaluate all the mechanisms, a few at a time."""
    if shutil.which("make") is None:
        pytest.fail("make is needed to build the evaluator")
    if not (REPO / "Submodules" / "amrex" / "Src" / "Base" / "AMReX.H").exists():
        pytest.fail("The AMReX submodule is needed: git submodule update --init")
    workdir = tmp_path_factory.mktemp("cantera")
    jobs = int(os.environ.get("CEPTR_TEST_JOBS", max(1, (os.cpu_count() or 1) // 2)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            name: pool.submit(evaluate, name, case, workdir)
            for name, case in CASES.items()
        }
    # A mechanism that fails fails only its own test
    return {
        name: future.exception() or future.result() for name, future in futures.items()
    }


@pytest.mark.parametrize("name", CASES)
def test_cantera(measured, name):
    """The generated mechanism agrees with Cantera, quantity by quantity."""
    case = CASES[name]
    reason, error = case.xfail or (None, None)
    if isinstance(measured[name], Exception):
        # Only the known failure: anything else is a failure of the test
        if error and error in str(measured[name]):
            pytest.xfail(reason)
        raise measured[name]
    assert not error, f"{name} is converted now: remove its xfail"
    errors, known_errors = measured[name]
    print(f"\n{name}: {case.covers}")
    failures, known, fixed = [], [], []
    for quantity, err in errors.items():
        tol = TOLERANCES[quantity]
        bound, why = case.known_quantities.get(quantity, (0.0, None))
        if err.error <= tol:
            status = ""
            if why:
                fixed.append(quantity)
        elif why and err.error <= bound:
            status = "known"
            known.append(why)
        else:
            status = "FAILED"
            failures.append(quantity)
        print(f"  {quantity:34s} {status:8s} (tol {tol:.0e}) {err}")
    # The known reactions of a case of their own must still differ
    for quantity, err in known_errors.items():
        tol = TOLERANCES[quantity]
        status = "left out"
        if reason:
            status = "known"
            if err.error <= tol:
                fixed.append(f"{quantity} of the known reactions")
            else:
                known.append(reason)
        print(f"  {quantity:34s} {status:8s} (tol {tol:.0e}) {err}")
    assert not failures, f"{name} differs from Cantera in {failures}"
    assert not fixed, f"{name} agrees now in {fixed}: update its known discrepancies"
    if known:
        pytest.xfail("; ".join(sorted(set(known))))
