"""Tests for heterogeneous (gas + surface) mechanisms."""

import pathlib
import re
import shutil
import subprocess

import cantera as ct
import numpy as np
import pytest
import yaml

import ceptr.converter as converter
import ceptr.utilities as cu

GAS_NAME = "gas"
INTERFACE_NAME = "Pt_surf"

# Minimal stand-ins for the AMReX headers that generated mechanisms include,
# so the generated code can be compiled without an AMReX checkout.
AMREX_STUB = """
#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#define AMREX_GPU_HOST_DEVICE
#define AMREX_FORCE_INLINE inline
#define AMREX_GPU_CONSTANT
#define AMREX_DEVICE_COMPILE 0
namespace amrex {
using Real = double;
template <class T> using Vector = std::vector<T>;
template <class T, std::size_t N> using GpuArray = std::array<T, N>;
inline void Abort(const std::string& msg) { throw std::runtime_error(msg); }
template <class... Ts> inline void ignore_unused(const Ts&...) {}
} // namespace amrex
"""

DRIVER = r"""
#include "mechanism.H"
#include <cstdio>
int main() {
  amrex::Real T, P, x[NUM_ALL_SPECIES], sc[NUM_ALL_SPECIES];
  if (std::scanf("%lf %lf", &T, &P) != 2) return 1;
  for (int i = 0; i < NUM_ALL_SPECIES; ++i)
    if (std::scanf("%lf", &x[i]) != 1) return 1;
  for (int i = 0; i < NUM_ALL_SPECIES; ++i)
    if (std::scanf("%lf", &sc[i]) != 1) return 1;
  amrex::Real qf[NUM_SURFACE_REACTIONS], qr[NUM_SURFACE_REACTIONS];
  comp_qfqr_surface(qf, qr, sc, T, 1.0 / T, std::log(T));
  for (int i = 0; i < NUM_SURFACE_REACTIONS; ++i)
    std::printf("QS %d %.17g %.17g\n", i, qf[i], qr[i]);
  amrex::Real sdot[NUM_ALL_SPECIES];
  surfaceProductionRate(sdot, sc, T);
  for (int i = 0; i < NUM_ALL_SPECIES; ++i) std::printf("SD %d %.17g\n", i, sdot[i]);
  amrex::Real kf[NUM_SURFACE_REACTIONS], kr[NUM_SURFACE_REACTIONS];
  SKKFKR(P, T, x, kf, kr);
  for (int i = 0; i < NUM_SURFACE_REACTIONS; ++i)
    std::printf("KK %d %.17g %.17g\n", i, kf[i], kr[i]);
  return 0;
}
"""


def find_cantera_data(fname):
    """Locate a data file shipped with Cantera."""
    for directory in ct.get_data_directories():
        path = pathlib.Path(directory) / fname
        if path.is_file():
            return path
    pytest.skip(f"{fname} is not available in this Cantera installation")


def convert(workdir):
    """Run CEPTR on workdir/mechanism.yaml.

    The output is not clang-formatted: the checks below do not depend on
    the layout.
    """
    interface = ct.Interface(workdir / "mechanism.yaml", INTERFACE_NAME)
    gas = interface.adjacent[GAS_NAME]
    conv = converter.Converter(gas, interface, "heterogeneous")
    conv.writer()
    return conv, workdir


@pytest.fixture(scope="module")
def ptcombust(tmp_path_factory):
    """Convert Cantera's ptcombust surface mechanism in a scratch directory."""
    workdir = tmp_path_factory.mktemp("ptcombust")
    shutil.copy(find_cantera_data("ptcombust.yaml"), workdir / "mechanism.yaml")
    return convert(workdir)


def modify_reaction(reactions, match, **changes):
    """Update the single reaction whose equation is match, in place."""
    (reaction,) = [r for r in reactions if r["equation"] == match]
    reaction.update(changes)
    return reaction


@pytest.fixture(scope="module")
def ptcombust_variant(tmp_path_factory):
    """Return ptcombust with rate features the stock mechanism does not exercise."""
    workdir = tmp_path_factory.mktemp("ptcombust_variant")
    data = yaml.safe_load(find_cantera_data("ptcombust.yaml").read_text())
    reactions = data["reactions"]
    # negative coverage exponent, evaluated at zero H(S) coverage below
    modify_reaction(
        reactions,
        "2 H(S) => H2 + 2 PT(S)",
        **{"coverage-dependencies": {"H(S)": {"a": 0.0, "m": -0.5, "E": -6000.0}}},
    )
    # negative order on a species that is not a reactant
    modify_reaction(
        reactions,
        "2 O(S) => O2 + 2 PT(S)",
        orders={"O(S)": 2.0, "PT(S)": -1.0},
        **{"negative-orders": True, "nonreactant-orders": True},
    )
    modify_reaction(reactions, "H2O + PT(S) => H2O(S)", **{"Motz-Wise": True})
    # reversible adsorption/desorption couples gas and surface Gibbs energies
    modify_reaction(
        reactions, "H2O(S) => H2O + PT(S)", equation="H2O(S) <=> H2O + PT(S)"
    )
    (workdir / "mechanism.yaml").write_text(yaml.safe_dump(data, sort_keys=False))
    return convert(workdir)


def defines(header):
    """Return the integer #defines of a generated header."""
    pattern = re.compile(r"^#define (\w+) (-?\d+)\b", re.M)
    return {k: int(v) for k, v in pattern.findall(header)}


def test_ptcombust_generation(ptcombust):
    """Generated sizes and element layout are consistent."""
    conv, workdir = ptcombust
    gas, interface = conv.mechanism, conv.interface
    hdr = (workdir / "mechanism.H").read_text()
    cpp = (workdir / "mechanism.cpp").read_text()
    d = defines(hdr)

    assert d["NUM_GAS_SPECIES"] == gas.n_species
    assert d["NUM_SURFACE_SPECIES"] == interface.n_species
    assert d["NUM_GAS_REACTIONS"] == gas.n_reactions
    assert d["NUM_SURFACE_REACTIONS"] == interface.n_reactions
    surface_only = [e for e in interface.element_names if e not in gas.element_names]
    assert d["NUM_GAS_ELEMENTS"] == gas.n_elements
    assert d["NUM_SURFACE_ELEMENTS"] == len(surface_only)

    # The gas-only API used by PelePhysics' EOS/transport keeps gas sizes.
    compact = re.sub(r"\s+", "", hdr)
    for name in ("ELEMENTS", "SPECIES", "REACTIONS"):
        assert f"#defineNUM_{name}NUM_GAS_{name}" in compact
        assert f"#defineNUM_ALL_{name}(NUM_GAS_{name}+NUM_SURFACE_{name})" in compact
    gas_names = dict(enumerate(gas.species_names))
    surf_names = dict(enumerate(interface.species_names))
    cksyms = cpp[cpp.index("CKSYMS_STR") : cpp.index("SKSYMS_STR")]
    sksyms = cpp[cpp.index("SKSYMS_STR") :]
    assert {
        int(i): n for i, n in re.findall(r'kname\[(\d+)\] = "([^"]+)";', cksyms)
    } == gas_names
    assert {
        int(i): n for i, n in re.findall(r'kname\[(\d+)\] = "([^"]+)";', sksyms)
    } == surf_names

    # Element names: CKSYME_STR gas only; SKSYME_STR gas elements first, then
    # surface-only elements.
    cksyme = cpp[cpp.index("CKSYME_STR") : cpp.index("CKSYMS_STR")]
    sksyme = cpp[cpp.index("SKSYME_STR") : cpp.index("SKSYMS_STR")]
    pattern = r'ename\[(\d+)\] = "(\w+)";'
    assert {int(i): n for i, n in re.findall(pattern, cksyme)} == dict(
        enumerate(gas.element_names)
    )
    assert {int(i): n for i, n in re.findall(pattern, sksyme)} == dict(
        enumerate(gas.element_names + surface_only)
    )

    # SKNCF rows are the surface species, columns the global element order.
    pt = len(gas.element_names) + surface_only.index("Pt")
    skncf = compact[compact.index("voidSKNCF(") :]
    skncf = skncf[: skncf.index("void", 4)]
    assert "intkd=NUM_ALL_ELEMENTS;" in skncf
    assert f"ncf[{interface.species_index('PT(S)')}*kd+{pt}]=1;" in skncf
    for function in (
        "comp_qfqr_surface",
        "surfaceProductionRate",
        "gibbs_surface",
        "SKKFKR",
        "SKAWT",
    ):
        assert f"void{function}(" in compact


class FakePhase:
    """Minimal stand-in for a Cantera phase's element data."""

    def __init__(self, element_names):
        self.element_names = list(element_names)


def test_element_numbering():
    """Surface-only elements follow the gas elements, in interface order."""
    gas = FakePhase(["O", "H", "C", "N", "Ar"])
    # Surface-only elements deliberately not listed first, in non-sorted order.
    interface = FakePhase(["H", "O", "Pt", "C", "Ni", "Al"])
    names = cu.get_element_names(gas, interface)
    assert names == ["O", "H", "C", "N", "Ar", "Pt", "Ni", "Al"]
    assert cu.get_element_names(gas, None) == gas.element_names


def compile_driver(workdir):
    """Compile the rate driver against the generated mechanism."""
    compiler = shutil.which("c++") or shutil.which("g++")
    if compiler is None:
        pytest.skip("no C++ compiler available")
    stub = workdir / "stub"
    stub.mkdir(exist_ok=True)
    for name in ("AMReX_Gpu.H", "AMReX_REAL.H"):
        (stub / name).write_text(AMREX_STUB)
    (workdir / "driver.cpp").write_text(DRIVER)
    exe = workdir / "driver"
    subprocess.run(
        [
            compiler,
            "-std=c++17",
            "-O0",
            f"-I{stub}",
            f"-I{workdir}",
            str(workdir / "driver.cpp"),
            str(workdir / "mechanism.cpp"),
            "-o",
            str(exe),
        ],
        check=True,
    )
    return exe


def max_rel_err(actual, expected):
    """Return the max relative error, ignoring entries negligible vs the max."""
    scale = 1e-12 * max(np.max(np.abs(expected)), 1e-300)
    return np.max(np.abs(actual - expected) / (np.abs(expected) + scale))


def check_rates_match_cantera(workdir, exe, bare_species=None):
    """Compiled surface rates agree with Cantera at random states.

    If bare_species is given, its coverage is set to zero.
    """
    interface = ct.Interface(workdir / "mechanism.yaml", INTERFACE_NAME)
    gas = interface.adjacent[GAS_NAME]
    names = gas.species_names + interface.species_names
    site_density = interface.site_density * 1000.0  # mol/m^2
    sizes = np.array([sp.size for sp in interface.species()])
    rng = np.random.default_rng(0)

    for _ in range(10):
        temp = rng.uniform(500.0, 1800.0)
        pres = ct.one_atm * rng.uniform(0.2, 5.0)
        mole_frac = rng.uniform(0.01, 1.0, gas.n_species)
        mole_frac /= mole_frac.sum()
        theta = rng.uniform(0.01, 1.0, interface.n_species)
        if bare_species is not None:
            theta[interface.species_index(bare_species)] = 0.0
        theta /= theta.sum()
        gas.TPX = temp, pres, mole_frac
        interface.TP = temp, pres
        interface.coverages = theta

        # Gas concentrations in mol/m^3, surface concentrations in mol/m^2.
        sc = np.concatenate([gas.concentrations * 1000.0, theta * site_density / sizes])
        # SKKFKR takes P in dyn/cm^2
        values = [temp, pres * 10.0, *mole_frac, *theta, *sc]
        out = subprocess.run(
            [str(exe)],
            input=" ".join(f"{v:.17g}" for v in values),
            capture_output=True,
            text=True,
            check=True,
        ).stdout

        rows = [line.split() for line in out.splitlines()]
        qs = {int(r[1]): (float(r[2]), float(r[3])) for r in rows if r[0] == "QS"}
        sd = {int(r[1]): float(r[2]) for r in rows if r[0] == "SD"}
        kk = {int(r[1]): (float(r[2]), float(r[3])) for r in rows if r[0] == "KK"}

        # surface rates are reported in mechanism order
        qf, qr = np.array([qs[i] for i in range(interface.n_reactions)]).T
        kf, kr = np.array([kk[i] for i in range(interface.n_reactions)]).T

        # Cantera works in kmol; the generated code in mol.
        ct_qf = interface.forward_rates_of_progress * 1000.0
        ct_qr = interface.reverse_rates_of_progress * 1000.0
        ct_sdot = dict(
            zip(interface.kinetics_species_names, interface.net_production_rates)
        )
        ct_sdot = np.array([ct_sdot[n] * 1000.0 for n in names])
        sdot = np.array([sd[i] for i in range(len(names))])

        tol = 1e-10
        assert max_rel_err(qf, ct_qf) < tol
        assert max_rel_err(qr, ct_qr) < tol
        assert max_rel_err(sdot, ct_sdot) < tol
        # SKKFKR returns surface rates in mol/cm^2/s.
        assert max_rel_err(kf * 1e4, ct_qf) < tol
        assert max_rel_err(kr * 1e4, ct_qr) < tol


def test_ptcombust_rates_match_cantera(ptcombust):
    """Stock ptcombust rates agree with Cantera."""
    workdir = ptcombust[1]
    check_rates_match_cantera(workdir, compile_driver(workdir))


def test_ptcombust_variant_rates_match_cantera(ptcombust_variant):
    """Rate features absent from the stock ptcombust agree with Cantera.

    Covers negative orders, negative coverage exponents on a bare surface,
    Motz-Wise sticking and reversible desorption.
    """
    _, workdir = ptcombust_variant
    exe = compile_driver(workdir)
    check_rates_match_cantera(workdir, exe)
    check_rates_match_cantera(workdir, exe, bare_species="H(S)")
    # zero concentration with a negative order gives a zero rate
    check_rates_match_cantera(workdir, exe, bare_species="PT(S)")


def test_unsupported_coverage_dependency(tmp_path):
    """Polynomial coverage activation energies are rejected before writing."""
    data = yaml.safe_load(find_cantera_data("ptcombust.yaml").read_text())
    modify_reaction(
        data["reactions"],
        "2 H(S) => H2 + 2 PT(S)",
        **{
            "coverage-dependencies": {
                "H(S)": {"a": 0.0, "m": 0.0, "E": [-6000.0, 1000.0, 0.0, 0.0]}
            }
        },
    )
    fname = tmp_path / "mechanism.yaml"
    fname.write_text(yaml.safe_dump(data, sort_keys=False))
    interface = ct.Interface(fname, INTERFACE_NAME)
    with pytest.raises(NotImplementedError, match="coverage"):
        converter.Converter(interface.adjacent[GAS_NAME], interface, "heterogeneous")
    assert not (tmp_path / "mechanism.H").exists()
