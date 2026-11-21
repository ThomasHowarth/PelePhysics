
#ifndef SOOTMODEL_H
#define SOOTMODEL_H

// AMReX include statements
#include <AMReX_Gpu.H>

// Soot include statements
#include "Constants_Soot.H"
#include "SootData.H"
#include "SootReactions.H"

class SootModel
{
public:
  //
  // Default constructor
  //
  SootModel();

  //
  // No copy constructors
  //
  SootModel(SootModel const&) = delete;

  //
  // No assignment operator
  //
  SootModel& operator=(SootModel const&) = delete;

  //
  // Destructor
  //
  ~SootModel() = default;

  //
  // Deallocate memory
  //
  void cleanup() // NOLINT(readability-make-member-function-const)
  {
    delete m_sootData;
    delete m_sootReact;
    amrex::The_Arena()->free(d_sootData);
    amrex::The_Arena()->free(d_sootReact);
  }

  //
  // Define instance of class
  //
  void define();

  //
  // Read soot parameters from input file
  //
  void readSootParams();

  //
  // Define member data
  //
  void defineMemberData(const amrex::Real dimerVol);

  //
  // Compute HMOM source term
  //
  void computeSootSourceTerm(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate,
    amrex::Array4<const amrex::Real> const& coeff_mu,
    amrex::Array4<amrex::Real> const& soot_state,
    const amrex::Real time,
    const amrex::Real dt,
    const bool pres_term = true) const;

  //
  // Optional inputs to computeSootSourceTerm, for codes that couple the soot
  // source to the rest of the step (e.g. inside SDC iterations)
  //
  struct SootForcing
  {
    // Constant forcing of the gas (density, species) and of the moments
    // during the integration, in the components and units of soot_state.
    // It is subtracted from the result, so only the soot source is added.
    amrex::Array4<const amrex::Real> forcing;
    // State at the end of the step: the temperature is held at the midpoint
    // of Qstate and Qstate_end instead of at Qstate
    amrex::Array4<const amrex::Real> Qstate_end;
    // Viscosity at the end of the step, averaged likewise
    amrex::Array4<const amrex::Real> coeff_mu_end;
  };

  //
  // Compute HMOM source term, integrating the soot ODEs with a forcing
  //
  void computeSootSourceTerm(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate,
    amrex::Array4<const amrex::Real> const& coeff_mu,
    amrex::Array4<amrex::Real> const& soot_state,
    const amrex::Real time,
    const amrex::Real dt,
    const bool pres_term,
    SootForcing const& forcing) const;

  //
  // Add the instantaneous HMOM source rate at Qstate (no integration), in the
  // same components and units as computeSootSourceTerm
  //
  void computeSootSourceRate(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate,
    amrex::Array4<const amrex::Real> const& coeff_mu,
    amrex::Array4<amrex::Real> const& soot_state,
    const bool pres_term = true) const;

  // What computeSootSourceImpl does in each cell
  enum class SrcMode {
    Subcycle,      // integrate with the explicit subcycles, add the source
    Instantaneous, // add the instantaneous source rate
    CvodeSetup,    // fill the CVODE initial state and cell constants only
    CvodeFinalize  // add the source from the CVODE end state
  };

  // Number of CVODE unknowns per cell: moments (and weight), gas species,
  // density; and of the constants held per cell during the integration
  static constexpr int cvNeq = NUM_SOOT_MOMENTS + 1 + NUM_SOOT_GS + 1;
  static constexpr int cvNconst = 7 + NUM_SOOT_MOMENTS + 1 + NUM_SOOT_GS + 1;

  // Shared kernel of computeSootSourceTerm and computeSootSourceRate (public
  // because it holds a device lambda). cv_y / cv_c are the CVODE state and
  // the cell constants, cell-major over vbox, for the Cvode modes.
  void computeSootSourceImpl(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate,
    amrex::Array4<const amrex::Real> const& coeff_mu,
    amrex::Array4<amrex::Real> const& soot_state,
    const amrex::Real dt,
    const bool pres_term,
    SootForcing const& forcing,
    const SrcMode mode,
    amrex::Real* cv_y = nullptr,
    amrex::Real* cv_c = nullptr) const;

  // soot.integrator = cvode: integrate the soot ODEs of all cells of vbox
  // with CVODE (error controlled) instead of the explicit subcycles
  void computeSootSourceCvode(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate,
    amrex::Array4<const amrex::Real> const& coeff_mu,
    amrex::Array4<amrex::Real> const& soot_state,
    const amrex::Real dt,
    const bool pres_term,
    SootForcing const& forcing) const;

  //
  // Estimate the soot time step
  //
  amrex::Real estSootDt(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate) const;

  //
  // Access gas phase species name
  //
  std::string gasSpeciesName(const int sootGasIndx)
  {
    return m_gasSpecNames[sootGasIndx];
  }

  /***********************************************************************
    Reaction functions (SootModel_react.cpp)
  ***********************************************************************/
  //
  // Initialize reaction and species data
  //
  void initializeReactData();

  /***********************************************************************
    Inline functions
  ***********************************************************************/
  //
  // Set the indices for primitive and conservative variables
  //
  void setIndices(const SootComps& sootIndx) { m_sootIndx = sootIndx; }

  //
  // Return string of the names of each variable
  //
  std::string sootVariableName(const int indx) const
  {
    AMREX_ASSERT(indx >= 0 && indx < NUM_SOOT_MOMENTS + 1);
    return m_sootVarName[indx];
  }

  //
  // Return pointer to soot data
  //
  SootData* getSootData() const { return m_sootData; }
  SootData* getSootData_d() const { return d_sootData; }

  /***********************************************************************
    Member data
  ***********************************************************************/

  /// Units for moments are L^3(x + 2/3y)/L^3
  /// Units for the weight of the delta function are 1/L^3
  /// Where L is cm for PeleC and m for PeleLM

  /// Verbosity for soot model
  int m_sootVerbosity = 0;
  /// Soot component indices
  SootComps m_sootIndx;
  bool m_setIndx = false;
  /// Pointer to soot data and soot specific functions
  SootData* m_sootData = nullptr;
  /// Pointer for soot reaction data
  SootReaction* m_sootReact = nullptr;
  /// Pointers for data on device
  SootData* d_sootData = nullptr;
  SootReaction* d_sootReact = nullptr;
  /// Bool for checking if parameters have been read from input file
  bool m_readSootParams = false;
  /// Bool for checking if member data has all been defined
  bool m_memberDataDefined = false;
  /// Mass is not conserved during OH oxidation
  /// Mass can be conserved again by adding remaining lost mass to H2
  bool m_conserveMass = false;

  /// Species index of PAH
  int m_PAHindx = -1;
  /// PAH inception path to be used, must be A2, A3, or A4
  std::string m_inceptPAH;
  /// Name of PAH used for inception
  std::string m_PAHname;
  amrex::Real m_gammaStick = 0.;

  // Name of soot variables to output (moments and weight of delta function)
  amrex::Vector<std::string> m_sootVarName;

  // Maximum allowable concentration for adjusting time step
  amrex::Real m_Xcutoff = 1.E-12;
  // Temperature cutoff, any temperatures below this will not be solved
  amrex::Real m_Tcutoff = 273.;
  // Maximum possible subcycles
  int m_maxSubcycles = 20;
  // Number of subcycles to use during source calculations
  int m_numSubcycles = 1;
  // Soot ODE integrator: false = explicit subcycles, true = CVODE
  bool m_useCvode = false;
  // CVODE relative tolerance, and absolute tolerance relative to the
  // largest magnitude of each unknown over the box
  amrex::Real m_cvodeRtol = 1.0e-6;
  amrex::Real m_cvodeAtol = 1.0e-10;
  // CVODE solver: true = BDF with matrix-free GMRES, false = Adams with
  // fixed-point iterations
  bool m_cvodeBDF = true;
  int m_cvodeMaxSteps = 100000;

  /***********************************************************************
    Reaction member data
  ***********************************************************************/

  /// Bool to ensure reaction data has been read and filled
  bool m_reactDataFilled = false;
  /// Number and names of species for surface reactions
  amrex::Vector<std::string> m_gasSpecNames;

  /************************************************************************
    Factor member data (to reduce computational costs)
  ************************************************************************/

  /// Frequently used constants
  /// Factors for beta
  amrex::Real m_betaDimerFact = 0.;
  amrex::Real m_betaNuclFact = 0.;
};

#endif
