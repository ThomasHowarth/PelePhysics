"""Tests for heterogeneous (gas + surface) mechanisms."""

import pathlib
import re
import shutil

import cantera as ct
import pytest

import ceptr.converter as converter
import ceptr.utilities as cu

GAS_NAME = "gas"
INTERFACE_NAME = "Pt_surf"


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
    assert "voidSKAWT(" in compact


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
