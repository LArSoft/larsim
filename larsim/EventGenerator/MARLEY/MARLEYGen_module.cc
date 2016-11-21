//////////////////////////////////////////////////////////////////////////////
/// \file MarleyGen_module.cc
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author sjgardiner@ucdavis.edu
//////////////////////////////////////////////////////////////////////////////

// standard library includes
#include <array>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TTree.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/InterpolationGrid.hh"
#include "marley/Logger.hh"
#include "marley/NeutrinoSource.hh"
#include "marley/NuclearReaction.hh"
#include "marley/StructureDatabase.hh"


using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;

// anonynmous namespace for helper functions, etc.
namespace {

  // We need to convert from MARLEY's energy units (MeV) to LArSoft's
  // (GeV) using this conversion factor
  constexpr double MeV_to_GeV = 1e-3;

  void source_check_positive(double x, const char* description,
    const char* source_type)
  {
    if (x <= 0.) throw cet::exception("MarleyGen")
      << "Non-positive " << description << " value defined for a "
      << source_type << " neutrino source";
  }

  void source_check_nonnegative(double x, const char* description,
    const char* source_type)
  {
    if (x < 0.) throw cet::exception("MarleyGen")
      << "Negative " << description << " value defined for a "
      << source_type << " neutrino source";
  }

  std::string neutrino_pdg_to_string(int pdg) {
    if (pdg == marley_utils::ELECTRON_NEUTRINO)
      return std::string("ve");
    else if (pdg == marley_utils::ELECTRON_ANTINEUTRINO)
      return std::string("vebar");
    else if (pdg == marley_utils::MUON_NEUTRINO)
      return std::string("vu");
    else if (pdg == marley_utils::MUON_ANTINEUTRINO)
      return std::string("vubar");
    else if (pdg == marley_utils::TAU_NEUTRINO)
      return std::string("vt");
    else if (pdg == marley_utils::TAU_ANTINEUTRINO)
      return std::string("vtbar");
    else return std::string("?");
  }

}

namespace evgen {
  class MarleyGen;
}

class evgen::MarleyGen : public art::EDProducer {

public:

  explicit MarleyGen(const fhicl::ParameterSet& p);
  virtual ~MarleyGen();

  virtual void produce(art::Event& e) override;
  virtual void beginRun(art::Run& run) override;
  virtual void reconfigure(const fhicl::ParameterSet& p) override;

  void add_marley_particles(simb::MCTruth& truth,
    const std::vector<marley::Particle*>& particles, bool track);
  void prepare_neutrino_source(const fhicl::ParameterSet& p);
  std::string find_file(const std::string& fileName,
    const std::string& fileType);
  void reconfigure_marley(const fhicl::ParameterSet& p);
  void reconfigure_vertex(const fhicl::ParameterSet& p);
  void prepare_reactions(const std::vector<std::string>& files);
  void prepare_structure(const std::vector<std::string>& files);
  void check_pdf_pairs(const std::vector<double>& Es,
    const std::vector<double>& PDs);
  InterpMethod get_interpolation_method(const std::string& rule);
  int neutrino_pdg(const std::string& nu);

  enum class vertex_type_t { kSampled, kFixed };

protected:

  vertex_type_t fVertexType;

  //std::string fMarleyLoggingLevel;

  TLorentzVector fVertexPosition;

  std::unique_ptr<marley::Generator> fMarleyGenerator;

  // unique_ptr to the current event created by MARLEY
  std::unique_ptr<marley::Event> fEvent;

  // the MARLEY event TTree
  TTree* fEventTree;

  // string stream used to capture logger output from MARLEY
  // and redirect it to the LArSoft logger
  std::stringstream fMarleyLogStream;

  // Function that selects a primary vertex location for each event. The
  // result is saved to fVertexPosition.
  void update_vertex_pos();
  // Initializes sampling machinery for vertex position sampling
  void initialize_vertex_pos_sampling(const fhicl::ParameterSet& p);
  // Loads ROOT dictionaries for the MARLEY Event and Particle classes.
  // This allows the module to write the generated events to a TTree.
  void load_marley_dictionaries();
  // Run, subrun, and event numbers from the art::Event being processed
  uint_fast32_t fRunNumber;
  uint_fast32_t fSubRunNumber;
  uint_fast32_t fEventNumber;
  // Discrete distribution object used to sample TPCs based on their active
  // masses
  std::unique_ptr<std::discrete_distribution<size_t> > fTPCDist;
  // RNG object used to sample TPCs
  std::mt19937_64 fTPCEngine;
};

//------------------------------------------------------------------------------
evgen::MarleyGen::MarleyGen(const fhicl::ParameterSet& p)
  : fVertexType(vertex_type_t::kSampled), fEvent(new marley::Event),
  fRunNumber(0), fSubRunNumber(0), fEventNumber(0), fTPCDist(nullptr)
{
  // TODO: set MARLEY logging level from FHiCL parameters
  //fMarleyLoggingLevel = p.get<std::string>("MarleyLoggingLevel", "INFO");
  marley::Logger::Instance().add_stream(fMarleyLogStream,
    marley::Logger::LogLevel::INFO);

  // Delay creating the MARLEY generator until now so that we'll be able to
  // see all of its logger messages.
  fMarleyGenerator = std::make_unique<marley::Generator>();

  // Configure the module (including MARLEY itself) using the FHiCL parameters
  this->reconfigure(p);

  const auto& seed_service = art::ServiceHandle<rndm::NuRandomService>();

  // Register the MARLEY generator with the seed service. For simplicity, we
  // use a lambda as the seeder function (see NuRandomService.h for details).
  // This allows the SeedService to automatically re-seed MARLEY whenever
  // necessary. The user can set an explicit seed for MARLEY in the FHiCL
  // configuration file using the "seed" parameter within the
  // "marley_parameters" table. If you need to get the seed for MARLEY from the
  // SeedService, note that we're using "MARLEY" as its generator instance
  // name.
  rndm::NuRandomService::seed_t marley_seed = seed_service->registerEngine(
    [gen = this->fMarleyGenerator.get()](
      rndm::NuRandomService::EngineId const& /* unused */,
      rndm::NuRandomService::seed_t lar_seed) -> void
    {
      // Since we're capturing a pointer to the MARLEY generator value,
      // double-check that it's still good. This will prevent segfaults.
      if (gen) {
        auto seed = static_cast<uint_fast64_t>(lar_seed);
        gen->reseed(seed);
      }
    },
    "MARLEY", p.get<fhicl::ParameterSet>("marley_parameters"), { "seed" }
  );

  // Unless I'm mistaken, the call to registerEngine should seed the generator
  // with the seed from the FHiCL configuration file if one is included, but it
  // doesn't appear to do so (as of 16 Aug 2016, larsoft v06_03_00). As a
  // workaround, I manually reseed the generator (if needed) here using the
  // result of the call to registerEngine, which will be the seed from the
  // FHiCL file if one was given.
  // TODO: figure out what's going on here, and remove this workaround as
  // needed
  uint_fast64_t marley_cast_seed = static_cast<uint_fast64_t>(marley_seed);
  if (marley_cast_seed != fMarleyGenerator->get_seed()) {
    fMarleyGenerator->reseed(marley_cast_seed);
  }

  // Get the FHiCL parameters controlling the vertex location(s)
  fhicl::ParameterSet vertex_params(p.get<fhicl::ParameterSet>("vertex",
    fhicl::ParameterSet()));

  // Also register the TPC sampling engine with the seed service. Use a similar
  // lambda for simplicity. If you need the seed later, get it from the seed
  // service using "MarleyGenTPCEngine" as the instance name.
  rndm::NuRandomService::seed_t tpc_seed = seed_service->registerEngine(
    [this](rndm::NuRandomService::EngineId const& /* unused */,
      rndm::NuRandomService::seed_t lar_seed) -> void
    {
      // Since we're capturing the this pointer by value, double-check that
      // it's still good (and that the pointer to the TPC engine is still
      // alive) before using it. This will prevent segfaults.
      if (this) {
        auto seed = static_cast<uint_fast64_t>(lar_seed);
        // Use the obtained seed to prepare the random number engine.  This is
        // an attempt to do a decent job, but optimally accomplishing this can
        // be tricky (see, for example,
        // http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
        std::seed_seq seed_sequence{seed};
        this->fTPCEngine.seed(seed_sequence);
      }
    },
    "MarleyGenTPCEngine", vertex_params, { "seed" }
  );

  // TODO: resolve the other workaround mentioned above, then fix this as well
  uint_fast64_t tpc_cast_seed = static_cast<uint_fast64_t>(tpc_seed);
  std::seed_seq tpc_seed_sequence{tpc_cast_seed};
  fTPCEngine.seed(tpc_seed_sequence);

  // Log initialization information from the MARLEY generator
  LOG_INFO("MARLEY") << fMarleyLogStream.str();
  fMarleyLogStream.str("");
  fMarleyLogStream.clear();

  // Do any needed setup of the MARLEY class dictionaries
  load_marley_dictionaries();

  // Create a ROOT TTree using the TFileService that will store the MARLEY
  // event objects
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("MARLEY Event Tree",
    "A tree of marley::Event objects");
  fEventTree->Branch("events", "marley::Event", fEvent.get());

  // Add branches that give the art::Event run, subrun, and event numbers for
  // easy match-ups between the MARLEY and art TTrees. All three are recorded
  // as 32-bit unsigned integers.
  fEventTree->Branch("run_number", &fRunNumber, "run_number/i");
  fEventTree->Branch("subrun_number", &fSubRunNumber, "subrun_number/i");
  fEventTree->Branch("event_number", &fEventNumber, "event_number/i");

  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();
}

//------------------------------------------------------------------------------
evgen::MarleyGen::~MarleyGen()
{
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;
  std::unique_ptr<sumdata::RunData>
    runcol(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runcol));
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::add_marley_particles(simb::MCTruth& truth,
  const std::vector<marley::Particle*>& particles, bool track)
{
  // Loop over the vector of MARLEY particles and add simb::MCParticle
  // versions of each of them to the MCTruth object.
  for (const marley::Particle* p : particles) {
    // Treat all of these particles as primaries, which have negative
    // track IDs by convention
    int trackID = -1*(truth.NParticles() + 1);

    int pdg = p->pdg_code();
    double mass = p->mass() * MeV_to_GeV;
    double px = p->px() * MeV_to_GeV;
    double py = p->py() * MeV_to_GeV;
    double pz = p->pz() * MeV_to_GeV;
    double E = p->total_energy() * MeV_to_GeV;
    TLorentzVector mom(px, py, pz, E);

    int status = 0; // don't track the particles in LArG4 by default
    if (track) status = 1;

    simb::MCParticle part(trackID /* trackID to use in Geant4 */, pdg,
      "MARLEY", -1 /* primary particle */, mass, status);

    part.AddTrajectoryPoint(fVertexPosition, mom);
    truth.Add(part);
  }
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::produce(art::Event& e)
{

  // Get the run, subrun, and event numbers from the current art::Event
  fRunNumber = e.run();
  fSubRunNumber = e.subRun();
  fEventNumber = e.event();

  std::unique_ptr< std::vector<simb::MCTruth> >
    truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;

  truth.SetOrigin(simb::kSuperNovaNeutrino);

  // Get the primary vertex location for this event
  update_vertex_pos();

  *fEvent = fMarleyGenerator->create_event();

  fEventTree->Fill();

  // Add the initial and final state particles to the MCTruth object.
  add_marley_particles(truth, fEvent->get_initial_particles(), false);
  add_marley_particles(truth, fEvent->get_final_particles(), true);

  // calculate a few parameters for the call to SetNeutrino
  const marley::Particle& nu = fEvent->projectile();
  const marley::Particle& lep = fEvent->ejectile();
  double qx = nu.px() - lep.px();
  double qy = nu.py() - lep.py();
  double qz = nu.pz() - lep.pz();
  double Enu = nu.total_energy();
  double Elep = lep.total_energy();
  double Q2 = qx*qx + qy*qy + qz*qz - std::pow(Enu - Elep, 2);

  // For definitions of Bjorken x, etc., a good reference is Mark Thomson's
  // set of slides on deep inelastic scattering (http://tinyurl.com/hcn5n6l)
  double bjorken_x = Q2 / (2 * fEvent->target().mass() * (Enu - Elep));
  double inelasticity_y = 1. - Elep / Enu;

  // Include the initial excitation energy of the final-state nucleus when
  // calculating W (the final-state invariant mass of the hadronic system)
  // since the other parameters (X, Y) also take into account the 2-2
  // scattering reaction only.
  const marley::Particle& res = fEvent->residue();
  double hadronic_mass_W = res.mass() + fEvent->Ex();

  // TODO: do a more careful job of setting the parameters here
  truth.SetNeutrino(
    simb::kCC, // change when MARLEY can handle NC
    simb::kUnknownInteraction, // not sure what the mode should be
    simb::kUnknownInteraction, // not sure what the interaction type should be
    marley_utils::get_nucleus_pid(18, 40), // Ar-40 PDG code
    marley_utils::NEUTRON, // nucleon PDG
    0, // MARLEY handles low enough energies that we shouldn't need HitQuark
    hadronic_mass_W * MeV_to_GeV,
    bjorken_x, // dimensionless
    inelasticity_y, // dimensionless
    Q2 * std::pow(MeV_to_GeV, 2)
  );

  truthcol->push_back(truth);

  e.put(std::move(truthcol));

  // Process the MARLEY logging messages (if any) captured by our
  // stringstream and forward them to the messagefacility logger
  std::string line;
  while(std::getline(fMarleyLogStream, line)) {
    LOG_INFO("MARLEY") << line;
  }

  // Reset the MARLEY log stream
  // This can be done more elegantly by move-assigning a default-constructed
  // std::stringstream object like this:
  // fMarleyLogStream = std::stringstream();
  // However, g++ versions < 5 don't implement the move assignment
  // operator for std::stringstream
  // (see discussion at http://stackoverflow.com/q/27152263/4081973).
  // The approach below is a reasonable workaround.
  fMarleyLogStream.str("");
  fMarleyLogStream.clear();
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::reconfigure(const fhicl::ParameterSet& p)
{
  //fMarleyLoggingLevel = p.get<std::string>("MarleyLoggingLevel", "INFO");

  fhicl::ParameterSet vertex_params(p.get<fhicl::ParameterSet>("vertex"));
  reconfigure_vertex(vertex_params);

  fhicl::ParameterSet marley_params(p.get<fhicl::ParameterSet>(
    "marley_parameters"));

  reconfigure_marley(marley_params);
}

//------------------------------------------------------------------------------
std::string evgen::MarleyGen::find_file(const std::string& fileName,
  const std::string& fileType)
{
  cet::search_path searchPath("FW_SEARCH_PATH");

  std::string fullName;
  searchPath.find_file(fileName, fullName);

  if (fullName.empty())
    throw cet::exception("MarleyGen")
      << "Cannot find MARLEY " << fileType << " data file '"
      << fileName << "\'\n";

  return fullName;
}

//------------------------------------------------------------------------------
int evgen::MarleyGen::neutrino_pdg(const std::string& nu) {

  int pdg = 0;

  bool bad = false;

  // Matches integers
  static const std::regex rx_int = std::regex("[-+]?[0-9]+");
  if (std::regex_match(nu, rx_int)) {
    pdg = std::stoi(nu);
    if (!marley::NeutrinoSource::pdg_is_allowed(pdg)) bad = true;
  }
  else if (!marley_utils::string_to_neutrino_pdg(nu, pdg)) {
    bad = true;
  }

  if (bad) throw cet::exception("MarleyGen") << "Invalid neutrino type"
    << " specification '" << nu << "' given for the MARLEY"
    << " neutrino source.";

  return pdg;
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::reconfigure_marley(const fhicl::ParameterSet& p)
{
  // Remove any previous reaction data and nuclear structure data owned by the
  // MARLEY generator
  fMarleyGenerator->clear_reactions();
  fMarleyGenerator->get_structure_db().clear();

  // Get the incident neutrino direction 3-vector from the FHiCL parameters
  fhicl::ParameterSet dir_params = p.get<fhicl::ParameterSet>("direction",
    fhicl::ParameterSet());
  double Dx = dir_params.get<double>("x", 0.);
  double Dy = dir_params.get<double>("y", 0.);
  double Dz = dir_params.get<double>("z", 1.);
  std::array<double, 3> dir_vec{ Dx, Dy, Dz };

  // Set the incident neutrino direction using the FHiCL settings
  fMarleyGenerator->set_neutrino_direction(dir_vec);

  std::vector<std::string> reactions = p.get<std::vector<std::string> >(
    "reactions");

  std::vector<std::string> structure = p.get<std::vector<std::string> >(
    "structure", std::vector<std::string>());

  fhicl::ParameterSet source_params = p.get<fhicl::ParameterSet>("source");

  prepare_reactions(reactions);
  prepare_structure(structure);
  prepare_neutrino_source(source_params);
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::prepare_reactions(const std::vector<std::string>& files)
{
  if (files.size() == 0) throw cet::exception("MarleyGen") << "At"
    << " least one reaction must be defined in the FHiCL configuration"
    << " file.\n";

  for (auto filename : files) {
    std::string full_name = find_file(filename, "reaction");

    MARLEY_LOG_INFO() << "Loading reaction data from file " << full_name;
    auto nr = std::make_unique<marley::NuclearReaction>(full_name,
      fMarleyGenerator->get_structure_db());

    MARLEY_LOG_INFO() << "Added reaction " << nr->get_description();

    fMarleyGenerator->add_reaction(std::move(nr));
  }
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::prepare_structure(const std::vector<std::string>& files)
{
  // Load data for all nuclides in each structure data file. Assume that
  // all files use the default TALYS format.
  for (auto filename : files) {
    std::string full_name = find_file(filename, "structure");

    auto& sdb = fMarleyGenerator->get_structure_db();

    std::set<int> nucleus_PDGs = sdb.find_all_nuclides(full_name);

    for (int pdg : nucleus_PDGs) {
      int Z = marley_utils::get_particle_Z(pdg);
      int A = marley_utils::get_particle_A(pdg);

      std::string trimmed_nucid = marley_utils::trim_copy(
        marley_utils::nuc_id(Z, A));

      MARLEY_LOG_INFO() << "Loading nuclear structure data for "
        << trimmed_nucid << " from file " << full_name;

      sdb.emplace_decay_scheme(pdg, full_name);
    }
  }
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::check_pdf_pairs(const std::vector<double>& Es,
  const std::vector<double>& PDs)
{
  if (Es.size() < 2) throw cet::exception("MarleyGen")
    << "The energy grid defined for the neutrino source has less than"
    << " two grid points";

  if (PDs.size() != Es.size()) throw cet::exception("MarleyGen")
    << "The energy and probability density vectors used to define the"
    << " neutrino source have different lengths.";

  if (Es.front() < 0.) throw cet::exception("MarleyGen")
    << "Negative minimum energy defined for the neutrino source";

  double sum_of_PDs = 0.;

  for (size_t i = 1; i < Es.size(); ++i) {
    if (Es.at(i) <= Es.at(i - 1)) throw cet::exception("MarleyGen")
      << "The grid point energies defined for the neutrino source"
      << " are not strictly increasing";
    double pd = PDs.at(i);
    if (pd < 0.) throw cet::exception("MarleyGen")
      << "Negative probability density grid point defined for the"
      << " neutrino source";
    sum_of_PDs += pd;
  }

  if (sum_of_PDs <= 0.) throw cet::exception("MarleyGen")
    << "All probability density grid point values are zero for"
    << " the neutrino source";
}


//------------------------------------------------------------------------------
InterpMethod evgen::MarleyGen::get_interpolation_method(
  const std::string& rule)
{
  // Try using the ENDF-style numerical codes first
  static const std::regex rx_nonneg_int("[0-9]+");

  if (std::regex_match(rule, rx_nonneg_int)) {
    int endf_interp_code = std::stoi(rule);
    if (endf_interp_code == 1) return InterpMethod::Constant;
    else if (endf_interp_code == 2) return InterpMethod::LinearLinear;
    else if (endf_interp_code == 3) return InterpMethod::LinearLog;
    else if (endf_interp_code == 4) return InterpMethod::LogLinear;
    else if (endf_interp_code == 5) return InterpMethod::LogLog;
  }

  // Interpolation rules may also be given as strings
  else if (rule == "const" || rule == "constant")
    return InterpMethod::Constant;
  else if (rule == "lin" || rule == "linlin")
    return InterpMethod::LinearLinear;
  else if (rule == "log" || rule == "loglog")
    return InterpMethod::LogLog;
  // linear in energy, logarithmic in probability density
  else if (rule == "linlog")
    return InterpMethod::LinearLog;
  // logarithmic in energy, linear in probability density
  else if (rule == "loglin")
    return InterpMethod::LogLinear;
  else throw cet::exception("MarleyGen") << "Invalid interpolation rule '"
    << rule << "' given in the neutrino source specification";

  // We shouldn't ever end up here, but return something just in case
  return InterpMethod::Constant;
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::prepare_neutrino_source(const fhicl::ParameterSet& p)
{
  std::string type = p.get<std::string>("type");
  std::string nu = p.get<std::string>("neutrino");

  // Particle Data Group code for the neutrino type produced by this source
  int pdg = neutrino_pdg(nu);

  std::unique_ptr<marley::NeutrinoSource> source;

  if (type == "mono" || type == "monoenergetic") {
    double energy = p.get<double>("energy");
    source_check_positive(energy, "energy", "monoenergetic");
    source = std::make_unique<marley::MonoNeutrinoSource>(pdg, energy);
    MARLEY_LOG_INFO() << "Created monoenergetic "
      << neutrino_pdg_to_string(pdg) << " source with"
      << " neutrino energy = " << energy << " MeV";
  }
  else if (type == "dar" || type == "decay-at-rest") {
    source = std::make_unique<marley::DecayAtRestNeutrinoSource>(pdg);
     MARLEY_LOG_INFO() << "Created muon decay-at-rest "
       << neutrino_pdg_to_string(pdg) << " source";
  }
  else if (type == "fd" || type == "fermi-dirac" || type == "fermi_dirac") {
    double Emin = p.get<double>("Emin");
    double Emax = p.get<double>("Emax");
    double temp = p.get<double>("temperature");
    double eta = p.get<double>("eta", 0.);

    source_check_nonnegative(Emin, "Emin", "Fermi-Dirac");
    source_check_positive(temp, "temperature", "Fermi-Dirac");

    if (Emax <= Emin) throw cet::exception("MarleyGen")
      << "Emax <= Emin for a Fermi-Dirac neutrino source";

    source = std::make_unique<marley::FermiDiracNeutrinoSource>(pdg, Emin,
      Emax, temp, eta);
    MARLEY_LOG_INFO() << "Created Fermi-Dirac "
      << neutrino_pdg_to_string(pdg) << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << Emin << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << Emax << " MeV";
    MARLEY_LOG_INFO() << "  temperature = " << temp << " MeV";
    MARLEY_LOG_INFO() << "  eta = " << eta;
  }
  else if (type == "bf" || type == "beta" || type == "beta-fit") {
    double Emin = p.get<double>("Emin");
    double Emax = p.get<double>("Emax");
    double Emean = p.get<double>("Emean");
    double beta = p.get<double>("beta", 4.5);

    source_check_nonnegative(Emin, "Emin", "beta-fit");
    source_check_positive(Emean, "Emean", "beta-fit");

    if (Emax <= Emin) throw cet::exception("MarleyGen")
      << "Emax <= Emin for a beta-fit neutrino source";

    source = std::make_unique<marley::BetaFitNeutrinoSource>(pdg, Emin,
      Emax, Emean, beta);
    MARLEY_LOG_INFO() << "Created beta-fit "
      << neutrino_pdg_to_string(pdg) << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << Emin << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << Emax << " MeV";
    MARLEY_LOG_INFO() << "  average energy = " << Emean << " MeV";
    MARLEY_LOG_INFO() << "  beta = " << beta;
  }
  else if (type == "hist" || type == "histogram") {
    std::vector<double> Es = p.get<std::vector<double> >("E_bin_lefts");
    std::vector<double> weights = p.get<std::vector<double> >("weights");
    double Emax = p.get<double>("Emax");

    // Add Emax to the grid
    Es.push_back(Emax);

    // Set the probability density at E = Emax to be zero (this ensures
    // that no energies outside of the histogram will be sampled)
    weights.push_back(0.);

    // Convert from bin weights to probability densities by dividing by the
    // width of each bin
    int jmax = Es.size() - 1;
    for (int j = 0; j < jmax; ++j) {
      // TODO: The check for monotonically increasing energies that we perform
      // in check_pdf_pairs() should ensure that we don't divide by zero here,
      // but this might be worth looking at more carefully in the future.
      weights.at(j) /= (Es.at(j + 1) - Es.at(j));
    }

    // Check the grid points for the probability density function
    check_pdf_pairs(Es, weights);

    // Create the source
    source = std::make_unique<marley::GridNeutrinoSource>(Es, weights, pdg,
      InterpMethod::Constant);
    MARLEY_LOG_INFO() << "Created histogram "
      << neutrino_pdg_to_string(pdg) << " source";
  }
  else if (type == "grid") {
    std::vector<double> energies = p.get<std::vector<double> >("energies");
    std::vector<double> PDs = p.get<std::vector<double> >("prob_densities");
    std::string rule = p.get<std::string>("rule", "linlin");

    InterpMethod method = get_interpolation_method(rule);

    check_pdf_pairs(energies, PDs);

    source = std::make_unique<marley::GridNeutrinoSource>(energies, PDs, pdg,
      method);
    MARLEY_LOG_INFO() << "Created grid "
      << neutrino_pdg_to_string(pdg) << " source";
  }
  else if (type == "th1" || type == "tgraph") {
    std::string tfile = p.get<std::string>("tfile");
    std::string namecycle = p.get<std::string>("namecycle");

    std::string full_name = find_file(tfile, "source");

    auto file = std::make_unique<TFile>(full_name.c_str(), "read");
    if (!file) throw cet::exception("MarleyGen")
      << "Failed to open the ROOT file '" << tfile << "'";

    if (type == "th1") {
      std::unique_ptr<TH1> th1(dynamic_cast<TH1*>(
        file->Get(namecycle.c_str())));

      // Prevent ROOT from auto-deleting the TH1 when the TFile is closed
      if (th1) th1->SetDirectory(nullptr);
      else throw cet::exception("MarleyGen")
        << "Failed to load the TH1 '" << namecycle << "'" << " from the ROOT"
        << " file '" << full_name << "'";

      // The histogram was successfully recovered from the
      // TFile, so read in the (energy low edge, bin weight)
      // ordered pairs. Keep the overflow bin so that we can
      // use its left edge as the maximum energy value. Skip
      // the unneeded underflow bin.
      const TAxis* x_axis = th1->GetXaxis();
      if (!x_axis) throw cet::exception("MarleyGen") << "Error finding x-axis"
        << " of ROOT histogram (TH1 object) with namecycle '" << namecycle
        << "' from the ROOT file '" << full_name << "'";

      // include the overflow bin but not the underflow bin
      size_t nbins = x_axis->GetNbins() + 1;
      std::vector<double> Es(nbins);
      std::vector<double> PDs(nbins);

      for (size_t b = 1; b <= nbins; ++b) { // underflow bin is bin 0
        Es.at(b - 1) = x_axis->GetBinLowEdge(b);
        // the content of the overflow bin is reset to zero below
        PDs.at(b - 1) = th1->GetBinContent(b);
      }

      // Convert the bin weights to probability densities (used
      // by our GridNeutrinoSource object) by dividing each bin
      // weight by the corresponding bin width.
      for (size_t c = 0; c < nbins - 1; ++c) {
        PDs.at(c) /= Es.at(c + 1) - Es.at(c);
      }

      // Assign zero probability density to the overflow bin's
      // left edge. This ensures that neutrino energies will be
      // sampled on the half-open interval [Elow, Ehigh), where
      // Elow is the left edge of the first bin and Ehigh is the
      // left edge of the overflow bin.
      PDs.back() = 0.;

      // Perform some standard checks to make sure that the probability
      // density function represented by these grid points is sane.
      check_pdf_pairs(Es, PDs);

      // Now that we've processed grid points, create the grid neutrino
      // source
      source = std::make_unique<marley::GridNeutrinoSource>(
        Es, PDs, pdg, InterpMethod::Constant);
      MARLEY_LOG_INFO() << "Created a th1 "
        << neutrino_pdg_to_string(pdg) << " source";
    }
    else if (type == "tgraph") {
      std::unique_ptr<TGraph> tg(dynamic_cast<TGraph*>(
        file->Get(namecycle.c_str())));

      if (!tg) throw cet::exception("MarleyGen")
        << "Failed to load the TGraph '" << namecycle
        << "' from the ROOT file '" << full_name << "'";

      size_t num_points = tg->GetN();
      std::vector<double> Es(num_points);
      std::vector<double> PDs(num_points);

      // Load the energies and PDF values into our vectors
      for(size_t p = 0; p < num_points; ++p) {
        tg->GetPoint(p, Es.at(p), PDs.at(p));
      }

      // Perform some standard checks to make sure that the probability
      // density function represented by these grid points is sane.
      check_pdf_pairs(Es, PDs);

      // Create a neutrino source based on the grid
      source = std::make_unique<marley::GridNeutrinoSource>(Es,
        PDs, pdg, InterpMethod::LinearLinear);
      MARLEY_LOG_INFO() << "Created a tgraph "
        << neutrino_pdg_to_string(pdg) << " source";
    }
    else throw cet::exception("MarleyGen") // shouldn't ever get here
      << "Unknown error encountered while creating ROOT neutrino source"
      << " for MARLEY.";
  }
  else throw cet::exception("MarleyGen") // bad neutrino source type
    << "Unrecognized MARLEY neutrino source type '" << type << "'";

  // Load the generator with the new source object
  fMarleyGenerator->set_source(std::move(source));
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::update_vertex_pos()
{
  // sample a new position if needed
  if (fVertexType == vertex_type_t::kSampled) {
    // Sample a TPC index using the active masses as weights
    size_t tpc_index = fTPCDist->operator()(fTPCEngine);

    // Get the dimensions of this TPC's active volume
    art::ServiceHandle<geo::Geometry> geom;
    const auto& tpc = geom->TPC(tpc_index);
    double minX = tpc.MinX();
    double maxX = tpc.MaxX();
    double minY = tpc.MinY();
    double maxY = tpc.MaxY();
    double minZ = tpc.MinZ();
    double maxZ = tpc.MaxZ();
    std::uniform_real_distribution<double>::param_type x_range(minX, maxX);
    std::uniform_real_distribution<double>::param_type y_range(minY, maxY);
    std::uniform_real_distribution<double>::param_type z_range(minZ, maxZ);

    // Sample a location uniformly over this volume
    static std::uniform_real_distribution<double> uniform_dist;
    double x = uniform_dist(fTPCEngine, x_range);
    double y = uniform_dist(fTPCEngine, y_range);
    double z = uniform_dist(fTPCEngine, z_range);
    LOG_INFO("MarleyGen") << "Sampled primary vertex in TPC #"
      << tpc_index << ", x = " << x << ", y = " << y << ", z = " << z;

    // Update the vertex position 4-vector
    fVertexPosition.SetXYZT(x, y, z, 0.); // TODO: add time sampling
  }

  // if we're using a fixed vertex position, we don't need to do any sampling
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::reconfigure_vertex(const fhicl::ParameterSet& p)
{
  std::string type = p.get<std::string>("type", "sampled");

  if (type == "sampled") {

    fVertexType = vertex_type_t::kSampled;

    // Get the active masses of each of the TPCs in the current geometry. Use
    // them as weights for sampling a TPC to use for the primary vertex.
    art::ServiceHandle<geo::Geometry> geom;
    std::vector<double> tpc_masses;
    size_t num_tpcs = geom->NTPC();
    for (size_t iTPC = 0; iTPC < num_tpcs; ++iTPC) {
      // For each TPC, use its active mass (returned in kg) as its sampling
      // weight
      tpc_masses.push_back(geom->TPC(iTPC).ActiveMass());
    }

    // Load the discrete distribution used to sample TPCs with the up-to-date
    // set of weights
    fTPCDist.reset(new std::discrete_distribution<size_t>(tpc_masses.begin(),
      tpc_masses.end()));
  }
  else if (type == "fixed") {

    fVertexType = vertex_type_t::kFixed;

    double Vx = p.get<double>("x");
    double Vy = p.get<double>("y");
    double Vz = p.get<double>("z");

    fVertexPosition.SetXYZT(Vx, Vy, Vz, 0.); // TODO: add time setting
  }
  else throw cet::exception("MarleyGen") << "Invalid vertex type '"
    << type << "' requested. Allowed values are 'sampled' and 'fixed'";
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::load_marley_dictionaries()
{

  static bool already_loaded_marley_dict = false;

  if (already_loaded_marley_dict) return;

  // Current (24 July 2016) versions of ROOT 6 require runtime
  // loading of headers for custom classes in order to use
  // dictionaries correctly. If we're running ROOT 6+, do the
  // loading here, and give the user guidance if there are any
  // problems.
  //
  // This is the same technique used in the MARLEY source code
  // for the executable (src/marley.cc). If you change how this
  // code works, please sync changes with the executable as well.
  if (gROOT->GetVersionInt() >= 60000) {
    LOG_INFO("MarleyGen") << "ROOT 6 or greater detected. Loading class"
      << " information\nfrom headers \"marley/Particle.hh\""
      << " and \"marley/Event.hh\"";
    TInterpreter::EErrorCode* ec = new TInterpreter::EErrorCode();
    gInterpreter->ProcessLine("#include \"marley/Particle.hh\"", ec);
    if (*ec != 0) throw cet::exception("MarleyGen") << "Error loading"
      << " MARLEY header Particle.hh. For MARLEY headers stored in"
      << " /path/to/include/marley/, please add /path/to/include"
      << " to your ROOT_INCLUDE_PATH environment variable and"
      << " try again.";
    gInterpreter->ProcessLine("#include \"marley/Event.hh\"");
    if (*ec != 0) throw cet::exception("MarleyGen") << "Error loading"
      << " MARLEY header Event.hh. For MARLEY headers stored in"
      << " /path/to/include/marley/, please add /path/to/include"
      << " to your ROOT_INCLUDE_PATH environment variable and"
      << " try again.";
  }

  // No further action is required for ROOT 5 because the compiled
  // dictionaries (which are linked to this module) contain all of
  // the needed information

  already_loaded_marley_dict = true;
}

DEFINE_ART_MODULE(evgen::MarleyGen)
