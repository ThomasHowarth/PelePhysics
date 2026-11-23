"""Tests for the Planck-mean absorption coefficients (no downloads)."""

import io
import shutil
import subprocess
import types

import cantera as ct
import numpy as np
import pytest

from ceptr import planck_mean as pm
from ceptr import radiation as crad


def line(temp, nu=2349.1, a_ul=400.0, g_u=3.0, e_l=1500.0, q=50.0):
    """Intensity S (cm/molecule) and cooling function W of one line at temp."""
    e_u = e_l + nu
    s = (
        g_u
        * a_ul
        / (8.0 * np.pi * pm.C * nu**2)
        * np.exp(-pm.C2 * e_l / temp)
        * (-np.expm1(-pm.C2 * nu / temp))
        / q
    )
    w = g_u * a_ul * pm.H * pm.C * nu * np.exp(-pm.C2 * e_u / temp) / (4.0 * np.pi * q)
    return s, w


def test_kirchhoff():
    """A line gives the same kp from its intensity and its cooling function."""
    nu = 2349.1
    for temp in (300.0, 1000.0, 2800.0):
        s, w = line(temp, nu=nu)
        kp_w = pm.kp_cooling_function([temp - 1, temp + 1], [w, w], [temp])
        kp_s = pm.kp_planck_sum(pm.planck_sum([nu], [s], temp), temp)
        np.testing.assert_allclose(kp_s, kp_w, rtol=1e-12)


def synthetic(temps):
    """Return a smooth kp(T) with a maximum, like that of CO2."""
    x = np.log(temps / 1000.0)
    return np.exp(-1.0 - 0.8 * x - 1.5 * x**2 + 0.3 * x**3 + 0.2 * np.sin(4.0 * x))


def test_fit_one_range():
    """A smooth kp fits in one range; outside, held below and power law above."""
    temps = np.arange(250.0, 2000.0 + 1.0, 10.0)
    x = np.log(temps / 1000.0)
    kp = np.exp(-1.0 - 0.8 * x - 1.5 * x**2 + 0.3 * x**3)
    fit = pm.fit_ln_kp(temps, kp)
    assert len(fit["ranges"]) == 1
    assert fit["max_relative_error"] < 1e-9
    np.testing.assert_allclose(pm.evaluate(fit, temps), kp, rtol=1e-9)
    np.testing.assert_allclose(
        pm.evaluate(fit, [100.0, 200.0]), pm.evaluate(fit, 250.0)[0]
    )
    k_max = pm.evaluate(fit, 2000.0)[0]
    np.testing.assert_allclose(
        pm.evaluate(fit, [3000.0, 4000.0]),
        k_max * (np.array([3000.0, 4000.0]) / 2000.0) ** fit["exponent"],
        rtol=1e-12,
    )
    assert fit["exponent"] <= 0.0


def test_fit_two_ranges():
    """A kink needs two ranges, continuous at the break."""
    temps = np.arange(250.0, 3000.0 + 1.0, 10.0)
    kp = synthetic(temps) * np.where(temps < 800.0, 1.0, (temps / 800.0) ** -1.5)
    fit = pm.fit_ln_kp(temps, kp, tol=1e-2)
    assert len(fit["ranges"]) == 2
    t_b = fit["ranges"][0]["T_max"]
    lo = np.exp(
        np.polynomial.polynomial.polyval(
            np.log(t_b / 1000.0), fit["ranges"][0]["coefficients"]
        )
    )
    hi = np.exp(
        np.polynomial.polynomial.polyval(
            np.log(t_b / 1000.0), fit["ranges"][1]["coefficients"]
        )
    )
    np.testing.assert_allclose(lo, hi, rtol=1e-12)
    np.testing.assert_allclose(pm.evaluate(fit, temps), kp, rtol=1e-2)


def test_fit_three_ranges():
    """A kp falling then rising, like that of CH, needs three ranges."""
    temps = np.arange(250.0, 5000.0 + 1.0, 25.0)
    # Power laws of slope -5, 0 and 3, with kinks at 1000 K and 2500 K
    t = np.clip(temps, 1000.0, None) / 1000.0
    kp = np.minimum(temps / 1000.0, 1.0) ** -5 * np.maximum(t / 2.5, 1.0) ** 3
    fit = pm.fit_ln_kp(temps, kp, tol=1e-2)
    assert len(fit["ranges"]) == 3
    np.testing.assert_allclose(pm.evaluate(fit, temps), kp, rtol=1e-2)
    for lo, hi in zip(fit["ranges"][:-1], fit["ranges"][1:], strict=True):
        x_b = np.log(lo["T_max"] / 1000.0)
        np.testing.assert_allclose(
            pm._polyval(x_b, lo["coefficients"]),
            pm._polyval(x_b, hi["coefficients"]),
            rtol=1e-12,
        )


def test_fit_rising_end_is_held():
    """A kp still rising at T_max is held above it, not extrapolated up."""
    temps = np.arange(250.0, 1000.0 + 1.0, 10.0)
    fit = pm.fit_ln_kp(temps, (temps / 250.0) ** 0.5)
    assert fit["exponent"] == 0.0


def test_isotopes():
    """ExoMol isotopologue formulas are read with their counts."""
    assert pm._isotopes("(1H)2(16O)") == [(1, "H", 2), (16, "O", 1)]
    assert pm._isotopes("(12C)(1H)3(19F)") == [(12, "C", 1), (1, "H", 3), (19, "F", 1)]


def test_references():
    """Every dataset is cited, and every line-list reference is used."""
    used = set()
    for ds in pm.DATASETS.values():
        refs = ds.references()
        assert refs[0] == pm.REFERENCES[ds.linelist][0]
        assert pm.EXOMOL_REF in refs
        assert (pm.EXOCROSS_REF in refs) == ds.cooling_function
        used.update(refs)
    every = {r for refs in pm.REFERENCES.values() for r in refs}
    assert used == every | {pm.EXOMOL_REF, pm.EXOCROSS_REF}


def test_coefficient_files():
    """The committed fits are complete and positive over their range."""
    fits = crad.load_fits()
    assert set(fits) == set(pm.DEFAULT_SPECIES)
    for name, data in fits.items():
        assert data["licence"] == pm.LICENCE
        assert data["references"] == pm.DATASETS[name].references()
        assert data["composition"] == pm.DATASETS[name].composition
        fit = data["fit"]
        assert fit["max_relative_error"] <= pm.TOLERANCE
        kp = pm.evaluate(fit, np.linspace(200.0, 3500.0, 50))
        assert np.all(np.isfinite(kp)) and np.all(kp > 0.0)


def gri30():
    """GRI-Mech 3.0, which has most of the radiating species."""
    gas = ct.Solution("gri30.yaml")
    info = types.SimpleNamespace(nonqssa_species_list=gas.species_names)
    return gas, info


def test_match_species():
    """Species are matched by name and composition."""
    gas, info = gri30()
    pairs = crad.match_species(gas, info.nonqssa_species_list, crad.load_fits())
    names = {name for _, name, _ in pairs}
    assert {"H2O", "CO2", "CO", "CH4", "C2H4", "C2H2", "NH3", "NO", "N2O"} <= names
    assert {"OH", "HCN"} <= names
    for idx, name, _ in pairs:
        assert info.nonqssa_species_list[idx] == name
    # A species of the same name but another composition is not matched
    fits = {"CO2": dict(crad.load_fits()["CO2"], composition={"C": 1, "O": 1})}
    assert crad.match_species(gas, info.nonqssa_species_list, fits) == []


def test_match_alias():
    """A species is matched by an alias of the fit, HONO2 for HNO3."""
    composition = {"H": 1.0, "N": 1.0, "O": 3.0}
    mech = types.SimpleNamespace(
        species=lambda name: types.SimpleNamespace(composition=composition)
    )
    fits = {"HNO3": crad.load_fits()["HNO3"]}
    pairs = crad.match_species(mech, ["N2", "HONO2"], fits)
    assert [(idx, name) for idx, name, _ in pairs] == [(1, "HONO2")]


HARNESS = """
#include <cmath>
#include <cstdio>
using std::exp;
using std::log;
namespace amrex {
using Real = %s;
// As in AMReX: one template parameter, so mixed float/double does not compile
template <class T>
const T& max(const T& a, const T& b) { return a < b ? b : a; }
template <class... Ts>
void ignore_unused(const Ts&...) {}
}
#define AMREX_GPU_HOST_DEVICE
#define AMREX_FORCE_INLINE inline
#define NUM_GAS_SPECIES %d
%s
int main()
{
  double T;
  amrex::Real kp[NUM_GAS_SPECIES];
  while (std::scanf("%%lf", &T) == 1) {
    get_planck_mean_kp(amrex::Real(T), kp);
    for (int k = 0; k < NUM_GAS_SPECIES; ++k) std::printf("%%.17e ", double(kp[k]));
    std::printf("\\n");
  }
}
"""


WARNINGS = ["-Wall", "-Wextra", "-Wshadow", "-Werror"]


@pytest.mark.skipif(shutil.which("c++") is None, reason="no C++ compiler")
@pytest.mark.parametrize("real", ["float", "double"])
@pytest.mark.parametrize("with_fits", [True, False])
def test_generated_code_compiles(tmp_path, real, with_fits):
    """The function compiles without warnings in single and double precision."""
    gas, info = gri30()
    code = io.StringIO()
    crad.planck_mean(code, gas, info, crad.load_fits() if with_fits else {})
    src = tmp_path / "kp.cpp"
    src.write_text(HARNESS % (real, gas.n_species, code.getvalue()))
    subprocess.run(["c++", "-std=c++17", *WARNINGS, "-fsyntax-only", src], check=True)


@pytest.mark.skipif(shutil.which("c++") is None, reason="no C++ compiler")
def test_generated_code(tmp_path):
    """The generated function compiles and agrees with evaluate()."""
    gas, info = gri30()
    fits = crad.load_fits()
    code = io.StringIO()
    crad.planck_mean(code, gas, info, fits)
    src = tmp_path / "kp.cpp"
    src.write_text(HARNESS % ("double", gas.n_species, code.getvalue()))
    exe = tmp_path / "kp"
    subprocess.run(["c++", "-std=c++17", *WARNINGS, "-o", exe, src], check=True)

    temps = [150.0, 300.0, 650.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 6000.0]
    out = subprocess.run(
        [exe],
        input=" ".join(map(str, temps)),
        capture_output=True,
        text=True,
        check=True,
    )
    kp = np.array([[float(v) for v in row.split()] for row in out.stdout.splitlines()])
    pairs = crad.match_species(gas, info.nonqssa_species_list, fits)
    for idx, _, data in pairs:
        np.testing.assert_allclose(
            kp[:, idx], pm.evaluate(data["fit"], temps), rtol=1e-9
        )
    others = set(range(gas.n_species)) - {idx for idx, _, _ in pairs}
    assert np.all(kp[:, sorted(others)] == 0.0)
