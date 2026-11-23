"""Fit Planck-mean absorption coefficients from the ExoMol database.

This is a maintenance utility (``planck-mean``): it downloads ExoMol data,
computes the Planck-mean absorption coefficient kp(T) of each species, fits
it, and writes one small coefficient file per species into
``ceptr/radiation_data``.
``convert --radiation`` then writes those fits into mechanism.H (see
ceptr/radiation.py). Nothing is downloaded during a conversion.

The Planck-mean absorption coefficient of a gas at temperature T is

    kp(T) = int kappa_nu(T) E_b,nu(T) dnu / (sigma T^4)

with kappa_nu net of stimulated emission. In local thermodynamic equilibrium,
Kirchhoff's law makes the numerator the total spontaneous emission of the gas,
which ExoMol tabulates as the cooling function W(T), the power emitted per
molecule and per steradian, erg/(s sr molecule):

    W(T) = 1/(4 pi Q(T)) sum_{u,l} A_ul h c nu_ul g_u exp(-c2 E_u / T)

so kp = N pi W(T) / (sigma T^4), N = P / (k_B T). Molecules without a cooling
function have super-lines instead: line intensities S (cm/molecule) at a
temperature, summed in narrow wavenumber bins, and

    kp = N sum_bins S E_b(nu, T) / (sigma T^4)

The two agree (for a line, S E_b(nu, T) = pi W), which the tests check.
Super-line files are read as they download and are not stored.

kp is per unit partial pressure at 1 bar, in 1/(cm bar). The fit is of
ln kp against x = ln(T / 1000 K), see fit_ln_kp: one or two polynomial
ranges between T_low and T_max (the maximum temperature of the line list,
capped at T_CAP), kp held constant below T_low, and continued as a power law
of T above T_max.

Data
----
The data come from the ExoMol database (www.exomol.com), released under the
Creative Commons Attribution-ShareAlike 4.0 International licence
(CC BY-SA 4.0). The fits derived here are CC BY-SA 4.0 as well, and the
references written with each fit must be cited when it is used.
"""

import argparse
import bz2
import html
import json
import pathlib
import re
import urllib.parse
import urllib.request

import numpy as np
import yaml

# CGS constants (CODATA 2018)
H = 6.62607015e-27  # erg s
C = 2.99792458e10  # cm/s
KB = 1.380649e-16  # erg/K
C2 = H * C / KB  # second radiation constant, cm K
SIGMA = 5.670374419e-5  # erg/(s cm^2 K^4)
P_REF = 1.0e6  # reference pressure, dyn/cm^2 (1 bar)

T_LOW = 250.0  # K, kp held constant below
T_CAP = 5000.0  # K, highest T_max of a fit
DT_SAMPLE = 10.0  # K, sampling of cooling functions for the fits
DEGREE = 6  # polynomial degree of a range
TOLERANCE = 5.0e-3  # largest relative error of a fit
MAX_RANGES = 3  # most polynomial ranges of a fit

EXOMOL = "https://www.exomol.com"
EXOMOL_DB = f"{EXOMOL}/db"
LICENCE = "CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)"
FIT_DIR = pathlib.Path(__file__).parent / "radiation_data"

REFERENCES = {
    "ExoMol2016": (
        "J. Tennyson et al., The ExoMol database: molecular line lists for "
        "exoplanet and other hot atmospheres, J. Mol. Spectrosc. 327, 73-94 "
        "(2016), doi:10.1016/j.jms.2016.05.002"
    ),
    "ExoCross2018": (
        "S. N. Yurchenko, A. F. Al-Refaie, J. Tennyson, ExoCross: a general "
        "program for generating spectra from molecular line lists, Astron. "
        "Astrophys. 614, A131 (2018), doi:10.1051/0004-6361/201732531"
    ),
    "POKAZATEL": (
        "O. L. Polyansky, A. A. Kyuberis, N. F. Zobov, J. Tennyson, S. N. "
        "Yurchenko, L. Lodi, ExoMol molecular line lists XXX: a complete "
        "high-accuracy line list for water, Mon. Not. R. Astron. Soc. 480, "
        "2597-2608 (2018), doi:10.1093/mnras/sty1877"
    ),
    "UCL-4000": (
        "S. N. Yurchenko, T. M. Mellor, R. S. Freedman, J. Tennyson, ExoMol "
        "line lists - XXXIX. Ro-vibrational molecular line list for CO2, Mon. "
        "Not. R. Astron. Soc. 496, 5282-5291 (2020), doi:10.1093/mnras/staa1874"
    ),
    "Li2015": (
        "G. Li, I. E. Gordon, L. S. Rothman, Y. Tan, S.-M. Hu, S. Kassi, "
        "A. Campargue, E. S. Medvedev, Rovibrational line lists for nine "
        "isotopologues of the CO molecule in the X 1Sigma+ ground electronic "
        "state, Astrophys. J. Suppl. Ser. 216, 15 (2015), "
        "doi:10.1088/0067-0049/216/1/15"
    ),
    "MM": (
        "S. N. Yurchenko, A. Owens, K. Kefala, J. Tennyson, ExoMol line lists "
        "- LVII: High accuracy ro-vibrational line list for methane (CH4), "
        "Mon. Not. R. Astron. Soc. 528, 3719-3729 (2024), "
        "doi:10.1093/mnras/stae148"
    ),
    "MaYTY": (
        "B. P. Mant, A. Yachmenev, J. Tennyson, S. N. Yurchenko, ExoMol "
        "molecular line lists - XXVII: spectra of C2H4, Mon. Not. R. Astron. "
        "Soc. 478, 3220-3232 (2018), doi:10.1093/mnras/sty1239"
    ),
    "aCeTY": (
        "K. L. Chubb, J. Tennyson, S. N. Yurchenko, ExoMol molecular line lists "
        "- XXXVII: spectra of acetylene, Mon. Not. R. Astron. Soc. 493, "
        "1531-1545 (2020), doi:10.1093/mnras/staa229"
    ),
    "CoYuTe": (
        "P. A. Coles, S. N. Yurchenko, J. Tennyson, ExoMol molecular line lists "
        "- XXXV. A rotation-vibration line list for hot ammonia, Mon. Not. R. "
        "Astron. Soc. 490, 4638-4647 (2019), doi:10.1093/mnras/stz2778"
    ),
    "CoYuTe-MARVEL": (
        "A. R. Al Derzi, T. Furtenbacher, J. Tennyson, S. N. Yurchenko, "
        "A. G. Csaszar, MARVEL analysis of the measured high-resolution spectra "
        "of 14NH3, J. Quant. Spectrosc. Radiat. Transf. 161, 117-130 (2015), "
        "doi:10.1016/j.jqsrt.2015.03.034"
    ),
    "XABC": (
        "Q. Qu, S. N. Yurchenko, J. Tennyson, ExoMol molecular line lists - "
        "XLII. Rovibronic molecular line list for the low-lying states of NO, "
        "Mon. Not. R. Astron. Soc. 504, 5768-5777 (2021), "
        "doi:10.1093/mnras/stab1154"
    ),
    "TYM": (
        "S. N. Yurchenko, T. M. Mellor, J. Tennyson, ExoMol line lists - LIX. "
        "High-temperature line list for N2O, Mon. Not. R. Astron. Soc. 534, "
        "1364-1375 (2024), doi:10.1093/mnras/stae2201"
    ),
    "MYTHOS": (
        "G. B. Mitev, C. A. Bowesman, J. Zhang, S. N. Yurchenko, J. Tennyson, "
        "ExoMol line lists - LXI. A trihybrid line list for rovibronic "
        "transitions of the hydroxyl radical (OH), Mon. Not. R. Astron. Soc. "
        "536, 3401-3420 (2025), doi:10.1093/mnras/stae2803"
    ),
    "Harris": (
        "G. J. Harris, J. Tennyson, B. M. Kaminsky, Ya. V. Pavlenko, "
        "H. R. A. Jones, Improved HCN/HNC linelist, model atmospheres and "
        "synthetic spectra for WZ Cas, Mon. Not. R. Astron. Soc. 367, 400-406 "
        "(2006), doi:10.1111/j.1365-2966.2005.09960.x"
    ),
    "Harris-Barber": (
        "R. J. Barber, J. K. Strange, C. Hill, O. L. Polyansky, G. Ch. Mellau, "
        "S. N. Yurchenko, J. Tennyson, ExoMol line lists - III. An improved hot "
        "rotation-vibration line list for HCN and HNC, Mon. Not. R. Astron. "
        "Soc. 437, 1828-1835 (2014), doi:10.1093/mnras/stt2011"
    ),
    "ExoAmes": (
        "D. S. Underwood, J. Tennyson, S. N. Yurchenko, X. Huang, D. W. "
        "Schwenke, T. J. Lee, S. Clausen, A. Fateev, ExoMol molecular line "
        "lists - XIV. The rotation-vibration spectrum of hot SO2, Mon. Not. R. "
        "Astron. Soc. 459, 3890-3899 (2016), doi:10.1093/mnras/stw849"
    ),
    "ExoAmes-MARVEL": (
        "R. Tobias et al., Critical evaluation of measured rotational-"
        "vibrational transitions of four sulphur isotopologues of S16O2, J. "
        "Quant. Spectrosc. Radiat. Transf. 208, 152-163 (2018), "
        "doi:10.1016/j.jqsrt.2018.01.006"
    ),
    "AYT2": (
        "A. A. A. Azzam, J. Tennyson, S. N. Yurchenko, O. V. Naumenko, ExoMol "
        "molecular line lists - XVI. The rotation-vibration spectrum of hot "
        "H2S, Mon. Not. R. Astron. Soc. 460, 4063-4074 (2016), "
        "doi:10.1093/mnras/stw1133"
    ),
    "AYT2-MARVEL": (
        "K. L. Chubb et al., Marvel analysis of the measured high-resolution "
        "rovibrational spectra of H2 32S, J. Quant. Spectrosc. Radiat. Transf. "
        "218, 178-186 (2018), doi:10.1016/j.jqsrt.2018.07.012"
    ),
    "OYT8": (
        "A. Owens, S. N. Yurchenko, J. Tennyson, ExoMol line lists - LVIII. "
        "High-temperature molecular line list of carbonyl sulfide (OCS), Mon. "
        "Not. R. Astron. Soc. 530, 4004-4015 (2024), doi:10.1093/mnras/stae1110"
    ),
    "AYYJ": (
        "A. Y. Adam, A. Yachmenev, S. N. Yurchenko, P. Jensen, Variationally "
        "computed IR line list for the methyl radical CH3, J. Phys. Chem. A 123, "
        "4755-4763 (2019), doi:10.1021/acs.jpca.9b02919"
    ),
    "MoLLIST": (
        "P. F. Bernath, MoLLIST: Molecular line lists, intensities and spectra, "
        "J. Quant. Spectrosc. Radiat. Transf. 240, 106687 (2020), "
        "doi:10.1016/j.jqsrt.2019.106687"
    ),
    "MoLLIST-CH": (
        "T. Masseron et al., CH in stellar atmospheres: an extensive linelist, "
        "Astron. Astrophys. 571, A47 (2014), doi:10.1051/0004-6361/201423956"
    ),
    "kNigHt": (
        "A. N. Perri, L. K. McKemmish, Full spectroscopic model and trihybrid "
        "experimental-perturbative-variational line list for NH, Mon. Not. R. "
        "Astron. Soc. 531, 3023-3033 (2024), doi:10.1093/mnras/stae1340"
    ),
    "KTPSYT": (
        "S. V. Kozlov, V. A. Terashkevich, E. A. Pazyuk, A. V. Stolyarov, S. N. "
        "Yurchenko, J. Tennyson, An upgraded line list for radiative transitions "
        "between the B 2Sigma+, A 2Pi, and X 2Sigma+ states of the CN radical, "
        "Astrophys. J. Suppl. Ser. 275, 29 (2024), doi:10.3847/1538-4365/ad8b40"
    ),
    "Coxon-Hajig": (
        "J. A. Coxon, P. G. Hajigeorgiou, Improved direct potential fit analyses "
        "for the ground electronic states of the hydrogen halides: HF/DF/TF, "
        "HCl/DCl/TCl, HBr/DBr/TBr and HI/DI/TI, J. Quant. Spectrosc. Radiat. "
        "Transf. 151, 133-154 (2015), doi:10.1016/j.jqsrt.2014.08.028"
    ),
    "Coxon-Hajig-Li": (
        "G. Li, I. E. Gordon, R. J. Le Roy, P. G. Hajigeorgiou, J. A. Coxon, P. "
        "F. Bernath, L. S. Rothman, Reference spectroscopic data for hydrogen "
        "halides. Part I: Construction and validation of the ro-vibrational "
        "dipole moment functions, J. Quant. Spectrosc. Radiat. Transf. 121, "
        "78-90 (2013), doi:10.1016/j.jqsrt.2013.02.005"
    ),
    "Coxon-Hajig-Somogyi": (
        "W. Somogyi, S. N. Yurchenko, A. Yachmenev, Calculation of electric "
        "quadrupole linestrengths for diatomic molecules: Application to the H2, "
        "CO, HF, and O2 molecules, J. Chem. Phys. 155, 214303 (2021), "
        "doi:10.1063/5.0063256"
    ),
    "OYKYT": (
        "A. Owens, A. Yachmenev, J. Kupper, S. N. Yurchenko, W. Thiel, The "
        "rotation-vibration spectrum of methyl fluoride from first principles, "
        "Phys. Chem. Chem. Phys. 21, 3496-3505 (2019), doi:10.1039/C8CP01721B"
    ),
    "AIJS": (
        "A. I. Pavlyuchko, S. N. Yurchenko, J. Tennyson, ExoMol molecular line "
        "lists - XI. The spectrum of nitric acid, Mon. Not. R. Astron. Soc. 452, "
        "1702-1706 (2015), doi:10.1093/mnras/stv1376"
    ),
}


class Dataset:
    """An ExoMol line list and the file(s) kp is computed from.

    Parameters
    ----------
    composition : dict
        Elemental composition, to match mechanism species
    linelist : str
        ExoMol line list name
    path : str
        Directory of the line list under EXOMOL_DB
    files : str or dict
        Cooling function file name, or {temperature: super-lines file name}
    t_max : float
        Maximum temperature of the line list, K
    isotopologue : str
        What the file covers: "natural abundance", or the main isotopologue,
        whose kp then stands for the gas (for CO the two agree to 0.2%)
    refs : list of str
        References of the line list: keys of REFERENCES, or full references
    aliases : list of str
        Other names of the species in mechanisms (CHEMKIN style)
    """

    def __init__(
        self, composition, linelist, path, files, t_max, isotopologue, refs, aliases=()
    ):
        self.composition = composition
        self.linelist = linelist
        self.path = path
        self.files = files
        self.t_max = t_max
        self.isotopologue = isotopologue
        self.refs = refs
        self.aliases = list(aliases)

    def url(self, name):
        """Download URL of one of the dataset files."""
        return f"{EXOMOL_DB}/{self.path}/{name}"

    def sources(self):
        """URLs of the files kp is computed from."""
        if isinstance(self.files, str):
            return [self.url(self.files)]
        return [self.url(f) for f in self.files.values()]

    def references(self):
        """Full references to cite for this dataset."""
        keys = [*self.refs, "ExoMol2016"]
        if isinstance(self.files, str):
            keys.append("ExoCross2018")
        return [REFERENCES.get(k, k) for k in keys]


def _super(iso, linelist, wavenumbers, temps):
    return {t: f"{iso}__{linelist}__{wavenumbers}__{t}K.super.bz2" for t in temps}


DATASETS = {
    "H2O": Dataset(
        {"H": 2, "O": 1},
        "POKAZATEL",
        "H2O/1H2-16O/POKAZATEL",
        "H2O__POKAZATEL__natabund.cf",
        6000.0,
        "natural abundance",
        ["POKAZATEL"],
    ),
    "CO2": Dataset(
        {"C": 1, "O": 2},
        "UCL-4000",
        "CO2/12C-16O2/UCL-4000",
        "12C-16O2__UCL-4000.cf",
        4000.0,
        "12C-16O2",
        ["UCL-4000"],
    ),
    "CO": Dataset(
        {"C": 1, "O": 1},
        "Li2015",
        "CO/12C-16O/Li2015",
        "CO__Li2015__natabund.cf",
        9000.0,
        "natural abundance",
        ["Li2015"],
    ),
    "CH4": Dataset(
        {"C": 1, "H": 4},
        "MM",
        "CH4/12C-1H4/MM",
        "12C-1H4__MM.cf",
        2000.0,
        "12C-1H4",
        ["MM"],
    ),
    "C2H4": Dataset(
        {"C": 2, "H": 4},
        "MaYTY",
        "C2H4/12C2-1H4/MaYTY",
        _super("12C2-1H4", "MaYTY", "00000-07100", range(200, 701, 100)),
        700.0,
        "12C2-1H4",
        ["MaYTY"],
    ),
    "C2H2": Dataset(
        {"C": 2, "H": 2},
        "aCeTY",
        "C2H2/12C2-1H2/aCeTY",
        "12C2-1H2__aCeTY.cf",
        2200.0,
        "12C2-1H2",
        ["aCeTY"],
    ),
    "NH3": Dataset(
        {"N": 1, "H": 3},
        "CoYuTe",
        "NH3/14N-1H3/CoYuTe",
        _super("14N-1H3", "CoYuTe", "0-20000", range(200, 2001, 100)),
        2000.0,
        "14N-1H3",
        ["CoYuTe", "CoYuTe-MARVEL"],
    ),
    "NO": Dataset(
        {"N": 1, "O": 1},
        "XABC",
        "NO/14N-16O/XABC",
        "14N-16O__XABC.cf",
        3500.0,
        "14N-16O",
        ["XABC"],
    ),
    "N2O": Dataset(
        {"N": 2, "O": 1},
        "TYM",
        "N2O/14N2-16O/TYM",
        "N2O__TYM__natabund.cf",
        2000.0,
        "natural abundance",
        ["TYM"],
    ),
    "OH": Dataset(
        {"O": 1, "H": 1},
        "MYTHOS",
        "OH/16O-1H/MYTHOS",
        "16O-1H__MYTHOS.cf",
        10000.0,
        "16O-1H",
        ["MYTHOS"],
    ),
    "HCN": Dataset(
        {"H": 1, "C": 1, "N": 1},
        "Harris",
        "HCN/1H-12C-14N/Harris",
        "1H-12C-14N__Harris.cf",
        4000.0,
        "1H-12C-14N",
        ["Harris", "Harris-Barber"],
    ),
    "CH3": Dataset(
        {"C": 1, "H": 3},
        "AYYJ",
        "CH3/12C-1H3/AYYJ",
        "12C-1H3__AYYJ.cf",
        2000.0,
        "12C-1H3",
        ["AYYJ"],
    ),
    "CH": Dataset(
        {"C": 1, "H": 1},
        "MoLLIST",
        "CH/12C-1H/MoLLIST",
        "12C-1H__MoLLIST.cf",
        8000.0,
        "12C-1H",
        ["MoLLIST", "MoLLIST-CH"],
    ),
    "NH": Dataset(
        {"N": 1, "H": 1},
        "kNigHt",
        "NH/14N-1H/kNigHt",
        "14N-1H__kNigHt.cf",
        10000.0,
        "14N-1H",
        ["kNigHt"],
    ),
    "CN": Dataset(
        {"C": 1, "N": 1},
        "KTPSYT",
        "CN/12C-14N/KTPSYT",
        "12C-14N__KTPSYT.cf",
        10000.0,
        "12C-14N",
        ["KTPSYT"],
    ),
    "HF": Dataset(
        {"H": 1, "F": 1},
        "Coxon-Hajig",
        "HF/1H-19F/Coxon-Hajig",
        "1H-19F__Coxon-Hajig.cf",
        5000.0,
        "1H-19F",
        ["Coxon-Hajig", "Coxon-Hajig-Li", "Coxon-Hajig-Somogyi"],
    ),
    "CH3F": Dataset(
        {"C": 1, "H": 3, "F": 1},
        "OYKYT",
        "CH3F/12C-1H3-19F/OYKYT",
        "12C-1H3-19F__OYKYT.cf",
        500.0,
        "12C-1H3-19F",
        ["OYKYT"],
    ),
    "HNO3": Dataset(
        {"H": 1, "N": 1, "O": 3},
        "AIJS",
        "HNO3/1H-14N-16O3/AIJS",
        "1H-14N-16O3__AIJS.cf",
        500.0,
        "1H-14N-16O3",
        ["AIJS"],
        aliases=["HONO2"],
    ),
    "SO2": Dataset(
        {"S": 1, "O": 2},
        "ExoAmes",
        "SO2/32S-16O2/ExoAmes",
        "32S-16O2__ExoAmes.cf",
        2000.0,
        "32S-16O2",
        ["ExoAmes", "ExoAmes-MARVEL"],
    ),
    "H2S": Dataset(
        {"H": 2, "S": 1},
        "AYT2",
        "H2S/1H2-32S/AYT2",
        "1H2-32S__AYT2.cf",
        1500.0,
        "1H2-32S",
        ["AYT2", "AYT2-MARVEL"],
    ),
    "OCS": Dataset(
        {"O": 1, "C": 1, "S": 1},
        "OYT8",
        "OCS/16O-12C-32S/OYT8",
        "16O-12C-32S__OYT8.cf",
        2000.0,
        "16O-12C-32S",
        ["OYT8"],
    ),
}


# Most abundant isotope of the elements of ExoMol molecules, to pick the main
# isotopologue
MAIN_ISOTOPE = {
    "H": 1, "He": 4, "Li": 7, "Be": 9, "B": 11, "C": 12, "N": 14, "O": 16,
    "F": 19, "Na": 23, "Mg": 24, "Al": 27, "Si": 28, "P": 31, "S": 32,
    "Cl": 35, "K": 39, "Ca": 40, "Sc": 45, "Ti": 48, "V": 51, "Cr": 52,
    "Fe": 56, "Ni": 58, "As": 75, "Br": 79, "Y": 89, "Zr": 90, "La": 139,
}  # fmt: skip

# DOIs of the ExoMol and ExoCross papers, which dataset pages list besides
# those of the line list
GENERIC_DOIS = {
    "10.1016/j.jms.2016.05.002",
    "10.1051/0004-6361/201732531",
    "10.1051/0004-6361/202038350",
    "10.1093/rasti/rzae016",
}


def _get(url):
    with _download(url) as response:
        return response.read()


def _isotopes(iso_formula):
    """[(mass number, element, count)] of an ExoMol formula like (1H)2(16O)."""
    return [
        (int(a), el, int(n or 1))
        for a, el, n in re.findall(r"\((\d+)([A-Z][a-z]?)\)(\d*)", iso_formula)
    ]


def citation(doi):
    """Format the reference of a DOI from Crossref."""
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    work = json.loads(_get(url))["message"]
    authors = work.get("author", [])
    names = [
        " ".join(
            [*(f"{p[0]}." for p in a.get("given", "").replace("-", " ").split()),
             a.get("family", "")]
        )
        for a in authors
    ]  # fmt: skip
    who = ", ".join(names) if len(names) <= 8 else f"{names[0]} et al."
    title = html.unescape(re.sub(r"<[^>]+>", "", " ".join(work.get("title", [""]))))
    title = re.sub(r"\s+", " ", title).replace(" ,", ",").strip()
    journal = html.unescape(
        (work.get("short-container-title") or work.get("container-title") or [""])[0]
    )
    year = (work.get("issued", {}).get("date-parts") or [[""]])[0][0]
    pages = work.get("page") or work.get("article-number", "")
    return (
        f"{who}, {title}, {journal} {work.get('volume', '')}, {pages} ({year}),"
        f" doi:{doi}"
    ).replace(" ,", ",")


def discover(formula):
    """Find ExoMol data for a molecule absent from DATASETS.

    The main isotopologue's datasets are tried, the recommended one first:
    a natural-abundance cooling function, then the isotopologue's, then
    super-lines. The references are those its ExoMol page lists, formatted
    from Crossref: check them, as pages do not always list every paper.
    """
    query = urllib.parse.urlencode({"molecule": formula, "datatype": "linelist"})
    isos = json.loads(_get(f"{EXOMOL}/api/?{query}"))
    if not isos:
        raise ValueError(f"{formula}: not in ExoMol")
    main = [
        k for k in isos if all(MAIN_ISOTOPE.get(el) == a for a, el, _ in _isotopes(k))
    ]
    iso_formula = (main or list(isos))[0]
    entry = isos[iso_formula]
    composition = {}
    for _, el, n in _isotopes(iso_formula):
        composition[el] = composition.get(el, 0) + n
    datasets = [
        (name, ds) for name, ds in entry["linelist"].items()
        if isinstance(ds, dict) and ds.get("files")
    ]  # fmt: skip
    datasets.sort(key=lambda d: not d[1].get("recommended"))
    for name, ds in datasets:
        path = re.search(r"/db/(.+)/[^/]+$", ds["files"][0]["url"])[1]
        slug = path.split("/")[1]
        page = _get(f"{EXOMOL}/data/molecules/{path}/").decode()
        links = sorted(
            set(re.findall(r'href="[^"]*/([^/"]+\.(?:cf|super\.bz2))"', page))
        )
        cfs = [f for f in links if f.endswith(".cf")]
        natabund = [f for f in cfs if "natabund" in f]
        supers = {
            int(m[1]): f for f in links if (m := re.search(r"__(\d+)K\.super\.bz2$", f))
        }
        if not (cfs or supers):
            continue
        info = json.loads(_get(f"{EXOMOL_DB}/{path}/{slug}__{name}.def.json"))
        t_max = info["dataset"].get("max_temperature")
        if not t_max:
            continue
        if natabund or cfs:
            files, isotopologue = (
                (natabund[0], "natural abundance") if natabund else (cfs[0], slug)
            )
        else:
            files = {t: f for t, f in sorted(supers.items()) if t <= t_max}
            isotopologue = slug
        dois = sorted(
            set(re.findall(r"doi(?:\.org)?/(10\.[^\"<> ]+)", page)) - GENERIC_DOIS
        )
        refs = [citation(d) for d in dois] or [
            f"See {EXOMOL}/data/molecules/{path}/ for the references"
        ]
        return Dataset(composition, name, path, files, t_max, isotopologue, refs)
    raise ValueError(f"{formula}: no ExoMol cooling function or super-lines")


def dataset(species):
    """Dataset of ``species``: from DATASETS, or found in ExoMol."""
    return DATASETS[species] if species in DATASETS else discover(species)


# Species of the PelePhysics mechanisms with ExoMol data, written to
# ceptr/radiation_data by default. Left out: N2, O2 and H2 (homonuclear, kp
# 7 to 100+ orders of magnitude below that of H2O), H2CO (CH2O) and H2O2
# (only full line lists, 120-190 GB).
DEFAULT_SPECIES = [
    "H2O", "CO2", "CO", "CH4", "C2H4", "C2H2", "NH3", "NO", "N2O", "OH",
    "HCN", "CH3", "CH", "NH", "CN", "HF", "CH3F", "HNO3",
]  # fmt: skip


def kp_cooling_function(t_cf, w_cf, temps):
    """Return kp, 1/(cm bar), from a cooling function W(t_cf), erg/(s sr molec)."""
    temps = np.asarray(temps, dtype=np.float64)
    w = np.interp(temps, t_cf, w_cf)
    return P_REF / (KB * temps) * np.pi * w / (SIGMA * temps**4)


def planck_sum(nu, s, temp):
    """Sum of S_i E_b(nu_i, T), erg/(s molecule), for S in cm/molecule."""
    nu = np.asarray(nu, dtype=np.float64)
    keep = nu > 0.0
    nu = nu[keep]
    eb = 2.0 * np.pi * H * C**2 * nu**3 / np.expm1(C2 * nu / temp)
    return float(np.sum(np.asarray(s, dtype=np.float64)[keep] * eb))


def kp_planck_sum(total, temp):
    """Return kp, 1/(cm bar), from a planck_sum at temperature temp."""
    return P_REF / (KB * temp) * total / (SIGMA * temp**4)


def _download(url):
    request = urllib.request.Request(url, headers={"User-Agent": "ceptr"})
    return urllib.request.urlopen(request, timeout=600)


def read_columns(text: bytes):
    """Read two whitespace-separated columns of numbers."""
    data = np.array(text.split(), dtype=np.float64)
    return data[0::2], data[1::2]


def super_line_sum(url, temp, chunk_bytes=1 << 24):
    """Planck sum of a bzip2 super-lines file, streamed from ``url``."""
    total = 0.0
    rest = b""
    decomp = bz2.BZ2Decompressor()
    with _download(url) as response:
        while chunk := response.read(chunk_bytes):
            data = rest + decomp.decompress(chunk)
            cut = data.rfind(b"\n") + 1
            rest = data[cut:]
            if cut:
                total += planck_sum(*read_columns(data[:cut]), temp)
    if rest.strip():
        total += planck_sum(*read_columns(rest), temp)
    return total


def sample(ds):
    """Temperatures (K) and kp (1/(cm bar)) of a Dataset from ExoMol."""
    t_max = min(ds.t_max, T_CAP)
    if isinstance(ds.files, str):
        with _download(ds.url(ds.files)) as response:
            t_cf, w_cf = read_columns(response.read())
        temps = np.arange(T_LOW, t_max + 0.5 * DT_SAMPLE, DT_SAMPLE)
        return temps, kp_cooling_function(t_cf, w_cf, temps)
    temps = np.array(sorted(ds.files), dtype=np.float64)
    kp = [kp_planck_sum(super_line_sum(ds.url(ds.files[t]), t), t) for t in temps]
    return temps, np.array(kp)


def _polyfit(x, y, degree):
    # Least squares on the scaled domain [-1, 1], then plain coefficients in x
    return np.polynomial.Polynomial.fit(x, y, degree).convert().coef


def _polyval(x, coeffs):
    return np.polynomial.polynomial.polyval(x, coeffs)


def _range_error(x, y, coeffs):
    return float(np.max(np.abs(np.expm1(_polyval(x, coeffs) - y))))


def _fit_ranges(x, y, deg, tol, n_ranges, y0=None):
    """Smallest largest error of ``n_ranges`` continuous ranges fitting y(x).

    The first range starts at y0 if given. Ranges whose error exceeds tol are
    not split further.

    Returns
    -------
    (float, list of (int, int, ndarray))
        Error, and the (first, last) sample index and coefficients of each range
    """
    if n_ranges == 1:
        coeffs = _polyfit(x, y, deg)
        if y0 is not None:
            coeffs[0] += y0 - _polyval(x[0], coeffs)
        return _range_error(x, y, coeffs), [(0, len(x) - 1, coeffs)]
    best = (np.inf, None)
    for i in range(deg + 2, len(x) - (deg + 2) * (n_ranges - 1)):
        err_lo, lo = _fit_ranges(x[: i + 1], y[: i + 1], deg, tol, 1, y0)
        if err_lo > tol:
            continue
        y_i = _polyval(x[i], lo[0][2])
        err_hi, hi = _fit_ranges(x[i:], y[i:], deg, tol, n_ranges - 1, y_i)
        err = max(err_lo, err_hi)
        if err < best[0]:
            best = (err, lo + [(a + i, b + i, c) for a, b, c in hi])
    return best


def fit_ln_kp(
    temps,
    kp,
    t_low=T_LOW,
    t_max=None,
    degree=DEGREE,
    tol=TOLERANCE,
    max_ranges=MAX_RANGES,
):
    """Fit ln kp as polynomials in x = ln(T / 1000 K).

    The fewest ranges (up to ``max_ranges``) of degree ``degree`` (at most the
    number of samples minus 2) within ``tol``: the breaks are at sample
    temperatures, chosen for the smallest error, and each range is shifted to
    be continuous with the one below.

    Returns
    -------
    dict
        T_low, T_max, ranges ([{T_min, T_max, coefficients}]), the exponent of
        the power law above T_max, and the largest relative error at the samples
    """
    temps = np.asarray(temps, dtype=np.float64)
    kp = np.asarray(kp, dtype=np.float64)
    t_max = temps[-1] if t_max is None else t_max
    use = (temps >= t_low - 1e-9) & (temps <= t_max + 1e-9)
    temps, x, y = temps[use], np.log(temps[use] / 1000.0), np.log(kp[use])
    t_low, t_max = float(temps[0]), float(temps[-1])

    deg = min(degree, len(x) - 2)
    for n_ranges in range(1, max_ranges + 1):
        err, pieces = _fit_ranges(x, y, deg, tol, n_ranges)
        if err <= tol:
            break
    else:
        raise ValueError(f"no fit of ln kp within {tol:.1%} in {max_ranges} ranges")
    ranges = [(temps[a], temps[b], c) for a, b, c in pieces]

    last = ranges[-1][2]
    deriv = np.polynomial.polynomial.polyder(last)
    exponent = min(0.0, float(_polyval(np.log(t_max / 1000.0), deriv)))
    return {
        "T_low": t_low,
        "T_max": t_max,
        "ranges": [
            {
                "T_min": float(a),
                "T_max": float(b),
                "coefficients": [float(v) for v in c],
            }
            for a, b, c in ranges
        ],
        "exponent": exponent,
        "max_relative_error": float(f"{err:.2e}"),
    }


def evaluate(fit, temps):
    """Evaluate a fit: kp, 1/(cm bar), at temperatures ``temps``, K."""
    temps = np.atleast_1d(np.asarray(temps, dtype=np.float64))
    x = np.log(np.maximum(temps, fit["T_low"]) / 1000.0)
    x_max = np.log(fit["T_max"] / 1000.0)
    ranges = fit["ranges"]
    breaks = np.log(np.array([r["T_max"] for r in ranges[:-1]]) / 1000.0)
    which = np.searchsorted(breaks, x, side="left")
    ln_kp = np.empty_like(x)
    for k, r in enumerate(ranges):
        sel = which == k
        ln_kp[sel] = _polyval(x[sel], r["coefficients"])
    above = x > x_max
    if above.any():
        y_max = _polyval(x_max, ranges[-1]["coefficients"])
        ln_kp[above] = y_max + fit["exponent"] * (x[above] - x_max)
    return np.exp(ln_kp)


HEADER = """\
# Planck-mean absorption coefficient kp(T) of {species}, per unit partial
# pressure, 1/(cm bar): a fit of ln kp to ExoMol data, written by planck-mean
# (Support/ceptr/ceptr/planck_mean.py). convert --radiation writes it into
# mechanism.H.
#
# Licence: derived from the ExoMol database (www.exomol.com) and, like it,
# released under {licence}.
# Cite the references below when using it.
"""


def write_fit(species, ds, temps, kp, directory=FIT_DIR):
    """Fit kp of ``species`` (Dataset ``ds``) and write its coefficient file."""
    fit = fit_ln_kp(temps, kp, t_max=min(ds.t_max, T_CAP))
    data = {
        "species": species,
        "aliases": ds.aliases,
        "composition": ds.composition,
        "units": "1/(cm bar)",
        "fit": {
            "variable": "x = ln(T / 1000 K)",
            "form": (
                "ln kp = polynomial in x on each range (increasing powers); "
                "kp held at its T_low value below T_low; "
                "kp = kp(T_max) (T / T_max)^exponent above T_max"
            ),
            **fit,
        },
        "source": {
            "database": "ExoMol",
            "linelist": ds.linelist,
            "isotopologue": ds.isotopologue,
            "linelist_max_temperature": ds.t_max,
            "files": ds.sources(),
        },
        "licence": LICENCE,
        "references": ds.references(),
    }
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{species}.yaml"
    with open(path, "w") as f:
        f.write(HEADER.format(species=species, licence=LICENCE))
        yaml.safe_dump(data, f, sort_keys=False, width=88)
    return path, fit


def main():
    """Fit Planck-mean absorption coefficients from ExoMol data."""
    parser = argparse.ArgumentParser(
        description=main.__doc__,
        epilog=(
            "Species are ExoMol formulas (H2O, CO2, ...). Those of DATASETS use the"
            " line lists chosen there; for others, the ExoMol data are looked up"
            " (check the references written for them)."
        ),
    )
    parser.add_argument(
        "species",
        nargs="*",
        default=DEFAULT_SPECIES,
        help=f"species (default: {' '.join(DEFAULT_SPECIES)})",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        default=FIT_DIR,
        help="directory of the coefficient files (default: ceptr/radiation_data)",
    )
    args = parser.parse_args()
    for sp in args.species:
        ds = dataset(sp)
        print(f"{sp}: ExoMol {ds.linelist} ({ds.isotopologue})", flush=True)
        temps, kp = sample(ds)
        path, fit = write_fit(sp, ds, temps, kp, args.output)
        print(
            f"  {len(fit['ranges'])} range(s), {fit['T_low']:.0f}-{fit['T_max']:.0f} K,"
            f" max error {fit['max_relative_error']:.2%},"
            f" exponent above {fit['exponent']:.2f}: {path}"
        )
        if sp not in DATASETS:
            print("  references found on the ExoMol page, check them:")
            for ref in ds.references():
                print(f"    {ref}")


if __name__ == "__main__":
    main()
