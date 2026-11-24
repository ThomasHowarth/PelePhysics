// Evaluates a CEPTR-generated mechanism at given states, for comparison with
// Cantera (Support/ceptr/tests/test_cantera.py).
//
// Usage: evaluate <states file>
//
// The states file starts with the species names, in the order the caller
// uses (Cantera's), then has one state per line: T [K], P [dyn/cm^2] and the
// mole fractions in that order. For every state the results are written to
// standard output, in CGS units, one quantity per line ("name v0 v1 ..."),
// with species in the caller's order and reactions in the mechanism file's
// order.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "mechanism.H"
#include "Simple.H"

namespace {

using pele::physics::EosType;
using pele::physics::transport::binary_pair_index;
using SimpleParm = pele::physics::transport::
  TransParm<EosType, pele::physics::transport::SimpleTransport>;

void
write(const std::string& name, const std::vector<amrex::Real>& v)
{
  std::printf("%s", name.c_str());
  for (const auto x : v) {
    std::printf(" %.17g", x);
  }
  std::printf("\n");
}

// The transport parameters that InitParm<SimpleTransport>::host_initialize
// sets, without the runtime options (defaults: mixture rules, no Soret)
void
init_simple_transport(SimpleParm& tparm)
{
  egtransetWT(tparm.wt.data());
  egtransetEPS(tparm.eps.data());
  egtransetSIG(tparm.sig.data());
  egtransetDIP(tparm.dip.data());
  egtransetPOL(tparm.pol.data());
  egtransetZROT(tparm.zrot.data());
  egtransetNLIN(tparm.nlin.data());
  egtransetCOFETA(tparm.fitmu.data());
  egtransetCOFLAM(tparm.fitlam.data());
  egtransetCOFD(tparm.fitdbin.data());
  pele::physics::transport::simple_transport_pair_epsilon(&tparm);
  for (int i = 0; i < NUM_SPECIES; ++i) {
    tparm.iwt[i] = 1.0 / tparm.wt[i];
  }
}

// Reads a value with strtod: operator>> rejects subnormal numbers, which
// equilibrium mole fractions can be
bool
read_value(std::istream& in, amrex::Real& v)
{
  std::string token;
  if (!(in >> token)) {
    return false;
  }
  char* end = nullptr;
  v = std::strtod(token.c_str(), &end);
  return end == token.c_str() + token.size();
}

// The species thermodynamics, nondimensional and per mass
using SpeciesThermo = void (*)(amrex::Real T, amrex::Real* v);
const std::pair<const char*, SpeciesThermo> species_thermo[] = {
  {"cp_R", [](amrex::Real T, amrex::Real* v) { cp_R(v, T); }},
  {"h_RT", [](amrex::Real T, amrex::Real* v) { speciesEnthalpy(v, T); }},
  {"s_R", [](amrex::Real T, amrex::Real* v) { speciesEntropy(v, T); }},
  {"g_RT", [](amrex::Real T, amrex::Real* v) { gibbs(v, T); }},
  {"cp_mass_k", [](amrex::Real T, amrex::Real* v) { CKCPMS(T, v); }},
  {"cv_mass_k", [](amrex::Real T, amrex::Real* v) { CKCVMS(T, v); }},
  {"h_mass_k", [](amrex::Real T, amrex::Real* v) { CKHMS(T, v); }},
  {"u_mass_k", [](amrex::Real T, amrex::Real* v) { CKUMS(T, v); }},
  {"s_mass_k", [](amrex::Real T, amrex::Real* v) { CKSMS(T, v); }},
};

} // namespace

int
main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <states file>\n";
    return 1;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "Cannot open " << argv[1] << "\n";
    return 1;
  }

  // Map the caller's species order to the mechanism's
  amrex::Vector<std::string> names;
  CKSYMS_STR(names);
  std::map<std::string, int> index;
  for (int k = 0; k < NUM_SPECIES; ++k) {
    index[names[k]] = k;
  }
  int n_species = 0;
  in >> n_species;
  if (n_species != NUM_SPECIES) {
    std::cerr << "Expected " << NUM_SPECIES << " species, got " << n_species
              << "\n";
    return 1;
  }
  std::vector<int> sp(NUM_SPECIES);
  for (int k = 0; k < NUM_SPECIES; ++k) {
    std::string name;
    in >> name;
    const auto it = index.find(name);
    if (it == index.end()) {
      std::cerr << "Unknown species " << name << "\n";
      return 1;
    }
    sp[k] = it->second;
  }
  // GET_RMAP and CKKFKR fill NUM_REACTIONS entries (gas and surface)
#if NUM_REACTIONS > 0
  std::vector<int> rmap(NUM_REACTIONS);
  GET_RMAP(rmap.data());
#endif

  SimpleParm tparm;
  init_simple_transport(tparm);

  // Species quantities in the caller's order
  auto species = [&](const amrex::Real* v) {
    std::vector<amrex::Real> out(NUM_SPECIES);
    for (int k = 0; k < NUM_SPECIES; ++k) {
      out[k] = v[sp[k]];
    }
    return out;
  };

  // A state that cannot be read (truncated line, not a number) is an error,
  // not the end of the states
  auto read_error = [&](const int n) {
    std::cerr << "Cannot read state " << n << " of " << argv[1] << "\n";
    return 1;
  };
  amrex::Real T = 0.0, P = 0.0;
  int n_state = 0;
  for (; read_value(in, T); ++n_state) {
    amrex::Real x[NUM_SPECIES], y[NUM_SPECIES], c[NUM_SPECIES];
    if (!read_value(in, P)) {
      return read_error(n_state);
    }
    for (int k = 0; k < NUM_SPECIES; ++k) {
      if (!read_value(in, x[sp[k]])) {
        return read_error(n_state);
      }
    }
    CKXTY(x, y);
    CKXTCP(P, T, x, c);

    std::printf("state\n");

    // Thermodynamics, per species and of the mixture
    amrex::Real v[NUM_SPECIES];
    for (const auto& [name, eval] : species_thermo) {
      eval(T, v);
      write(name, species(v));
    }
    amrex::Real rho = 0.0, s = 0.0;
    CKRHOX(P, T, x, rho);
    write("density", {rho});
    CKMMWX(x, s);
    write("mean_molecular_weight", {s});
    CKCPBS(T, y, s);
    write("cp_mass", {s});
    CKCVBS(T, y, s);
    write("cv_mass", {s});
    CKHBMS(T, y, s);
    write("enthalpy_mass", {s});
    CKUBMS(T, y, s);
    write("int_energy_mass", {s});
    CKSBMS(P, T, y, s);
    write("entropy_mass", {s});

    // Kinetics: rates of progress in the mechanism file's reaction order and
    // net production rates
#if NUM_REACTIONS > 0
    {
      std::vector<amrex::Real> qf(NUM_REACTIONS), qr(NUM_REACTIONS);
      CKKFKR(P, T, x, qf.data(), qr.data());
      std::vector<amrex::Real> fwd(NUM_REACTIONS), rev(NUM_REACTIONS);
      for (int j = 0; j < NUM_REACTIONS; ++j) {
        fwd[rmap[j]] = qf[j];
        rev[rmap[j]] = qr[j];
      }
      write("forward_rates_of_progress", fwd);
      write("reverse_rates_of_progress", rev);
      amrex::Real wdot[NUM_SPECIES];
      CKWC(T, c, wdot);
      write("net_production_rates", species(wdot));
    }
#endif

    // Transport: the generated fits (pure species viscosities and
    // conductivities, binary diffusion coefficients at 1 atm), exp of cubics
    // in log(T), and the mixture values of Simple transport
    const amrex::Real logT = std::log(T);
    auto fit = [logT](const amrex::Real* cof) {
      return std::exp(
        cof[0] + logT * (cof[1] + logT * (cof[2] + logT * cof[3])));
    };
    auto species_fit = [&](const amrex::Real* cof) {
      for (int k = 0; k < NUM_SPECIES; ++k) {
        v[k] = fit(&cof[4 * k]);
      }
      return species(v);
    };
    write("species_viscosities", species_fit(tparm.fitmu.data()));
    write("species_conductivities", species_fit(tparm.fitlam.data()));
    std::vector<amrex::Real> dbin(NUM_SPECIES * NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; ++i) {
      for (int j = 0; j < NUM_SPECIES; ++j) {
        dbin[i * NUM_SPECIES + j] =
          fit(&tparm.fitdbin[4 * binary_pair_index(sp[i], sp[j])]);
      }
    }
    write("binary_diff_coeffs_1atm", dbin);

    amrex::Real mu = 0.0, xi = 0.0, lam = 0.0;
    amrex::Real rhoD[NUM_SPECIES] = {0.0}, chi[NUM_SPECIES] = {0.0};
    pele::physics::transport::SimpleTransport::transport(
      false, true, true, true, false, T, rho, y, rhoD, chi, mu, xi, lam,
      &tparm);
    write("viscosity", {mu});
    write("thermal_conductivity", {lam});
    write("rho_mix_diff_coeffs", species(rhoD));
  }
  // The loop ends at the end of the file, or on a temperature it cannot read
  if (!in.eof()) {
    return read_error(n_state);
  }
  return 0;
}
