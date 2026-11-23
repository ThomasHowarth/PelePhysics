"""Planck-mean absorption coefficients in mechanism.H (convert --radiation).

The fits come from the coefficient files in ceptr/radiation_data, written by
the planck-mean utility (ceptr/planck_mean.py) from the ExoMol database. They are
derived from ExoMol data and released under CC BY-SA 4.0: the generated code
carries the licence notice and the references.
"""

import math
import pathlib

import yaml

import ceptr.writer as cw

FIT_DIR = pathlib.Path(__file__).parent / "radiation_data"


def load_fits(directory=FIT_DIR):
    """Coefficient files of ``directory``, by species name."""
    fits = {}
    for path in sorted(pathlib.Path(directory).glob("*.yaml")):
        with open(path) as f:
            data = yaml.safe_load(f)
        fits[data["species"]] = data
    return fits


def match_species(mechanism, species_list, fits):
    """Pair gas species with fits: same name or alias (any case), same composition.

    Returns
    -------
    list of (int, str, dict)
        Index in ``species_list``, mechanism name, and coefficient file data
    """
    pairs = []
    by_name = {
        alias.upper(): data
        for name, data in fits.items()
        for alias in [name, *data.get("aliases", [])]
    }
    for idx, name in enumerate(species_list):
        data = by_name.get(name.upper())
        if data is None:
            continue
        composition = {
            el: n for el, n in mechanism.species(name).composition.items() if n
        }
        if composition == {el: float(n) for el, n in data["composition"].items()}:
            pairs.append((idx, name, data))
    return pairs


def _horner(coeffs, var="x"):
    expr = f"{coeffs[-1]:.17g}"
    for c in reversed(coeffs[:-1]):
        expr = f"{c:.17g} + {var} * ({expr})"
    return expr


def _ln_kp_code(fit):
    """C++ statements setting ln_kp from x = ln(T / 1000 K) for one fit."""
    x_low = math.log(fit["T_low"] / 1000.0)
    x_max = math.log(fit["T_max"] / 1000.0)
    ranges = fit["ranges"]
    last = ranges[-1]["coefficients"]
    y_max = sum(c * x_max**k for k, c in enumerate(last))
    lines = [
        f"const amrex::Real xs = amrex::max(x, amrex::Real({x_low:.17g}));",
        "amrex::Real ln_kp = 0.0;",
        f"if (xs > {x_max:.17g}) {{",
        f"ln_kp = {y_max:.17g} + {fit['exponent']:.17g} * (xs - {x_max:.17g});",
    ]
    for k, r in enumerate(ranges):
        if k < len(ranges) - 1:
            x_hi = math.log(r["T_max"] / 1000.0)
            lines.append(f"}} else if (xs <= {x_hi:.17g}) {{")
        else:
            lines.append("} else {")
        lines.append(f"ln_kp = {_horner(r['coefficients'], 'xs')};")
    lines.append("}")
    return lines


def planck_mean(fstream, mechanism, species_info, fits=None):
    """Write get_planck_mean_kp for the mechanism species that have fits."""
    fits = load_fits() if fits is None else fits
    species_list = species_info.nonqssa_species_list
    pairs = match_species(mechanism, species_list, fits)

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Planck-mean absorption coefficients kp(T)"))
    cw.writer(
        fstream,
        cw.comment(
            "kp[k]: species k per unit partial pressure, 1/(cm bar); 0 for"
            " species without data"
        ),
    )
    if pairs:
        cw.writer(
            fstream,
            cw.comment(
                "Fits of ExoMol data (www.exomol.com), released like it under"
                " CC BY-SA 4.0"
            ),
        )
        cw.writer(
            fstream,
            cw.comment("(https://creativecommons.org/licenses/by-sa/4.0/)."),
        )
        cw.writer(fstream, cw.comment("Cite the references of each species."))
    cw.writer(fstream, f"#define NUM_PLANCK_MEAN_SPECIES {len(pairs)}")
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
        " get_planck_mean_kp(const amrex::Real T, amrex::Real kp[])",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "for (int k = 0; k < NUM_GAS_SPECIES; ++k) {")
    cw.writer(fstream, "kp[k] = 0.0;")
    cw.writer(fstream, "}")
    if pairs:
        cw.writer(fstream, "const amrex::Real x = log(T * amrex::Real(1.0e-3));")
    else:
        cw.writer(fstream, "amrex::ignore_unused(T);")
    for idx, name, data in pairs:
        fit = data["fit"]
        src = data["source"]
        cw.writer(fstream)
        cw.writer(
            fstream,
            cw.comment(
                f"{name}: ExoMol {src['linelist']} ({src['isotopologue']}),"
                f" fit {fit['T_low']:.0f}-{fit['T_max']:.0f} K, max error"
                f" {fit['max_relative_error']:.2%}; held below, power law"
                " above"
            ),
        )
        for ref in data["references"]:
            cw.writer(fstream, cw.comment(f"  {ref}"))
        cw.writer(fstream, "{")
        for line in _ln_kp_code(fit):
            cw.writer(fstream, line)
        cw.writer(fstream, f"kp[{idx}] = exp(ln_kp);")
        cw.writer(fstream, "}")
    cw.writer(fstream, "}")
