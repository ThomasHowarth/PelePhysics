"""Tests for the committed-mechanism check."""

from ceptr.check_mechanisms import compare


def test_identical():
    """Identical files have no differences."""
    assert compare("m.H", "a = 1.0;\n", "a = 1.0;\n", 1e-6) == []


def test_last_digit_tolerated():
    """A last-digit change of a printed coefficient is tolerated."""
    old = "  COFD[6263] = 1.93525132E-06;\n"
    new = "  COFD[6263] = 1.93525134E-06;\n"
    assert compare("m.H", old, new, 1e-6) == []
    assert compare("m.H", old, new, 1e-9) != []


def test_value_change_detected():
    """A real change of a value is reported."""
    assert compare("m.H", "k = 3.17e-13;\n", "k = 3.18e-13;\n", 1e-6) != []


def test_reordering_detected():
    """Reordered operations are reported even if the numbers match."""
    old = "wdot[12] -= qdot;\nwdot[12] += 0.7 * qdot;\n"
    new = "wdot[12] += 0.7 * qdot;\nwdot[12] -= qdot;\n"
    assert compare("m.H", old, new, 1e-6) != []


def test_index_change_detected():
    """Integers (array indices, species ids) must match exactly in practice."""
    assert compare("m.H", "sc[12] * x;\n", "sc[13] * x;\n", 1e-6) != []


def test_line_count_change_detected():
    """Added or removed lines are reported."""
    assert compare("m.H", "a;\n", "a;\nb;\n", 1e-6) != []


def test_qssa_yaml_header_ignored():
    """Cantera's volatile qssa.yaml header lines are ignored."""
    old = "cantera-version: 3.1.0\ngit-commit: '6e00275'\ndate: Fri Jul 18\nx: 1\n"
    new = "cantera-version: 3.2.0\ngit-commit: 4a8358e\ndate: Wed Nov 11\nx: 1\n"
    assert compare("qssa.yaml", old, new, 1e-6) == []
    assert compare("qssa.yaml", old, new.replace("x: 1", "x: 2"), 1e-6) != []


def test_large_integer_exact():
    """Integer literals are compared exactly, whatever the tolerance."""
    assert compare("m.H", "J[1000000] = x;\n", "J[1000001] = x;\n", 1e-3) != []
