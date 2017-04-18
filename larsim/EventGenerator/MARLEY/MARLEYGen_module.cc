//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYGen_module.cc
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author sjgardiner@ucdavis.edu
//////////////////////////////////////////////////////////////////////////////

// standard library includes
#include <algorithm>
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
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/detail/ParameterArgumentTypes.h"
#include "fhiclcpp/types/detail/ParameterBase.h"
#include "fhiclcpp/types/detail/TableBase.h"
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
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"


// anonynmous namespace for helper functions, etc.
namespace {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  // We need to convert from MARLEY's energy units (MeV) to LArSoft's
  // (GeV) using this conversion factor
  constexpr double MeV_to_GeV = 1e-3;

  template <typename T> marley::JSON fhicl_atom_to_json(
    const fhicl::Atom<T>* atom)
  {
    return marley::JSON(atom->operator()());
  }

  template <typename T> marley::JSON fhicl_optional_atom_to_json(
    const fhicl::OptionalAtom<T>* opt_atom)
  {
    T value;
    if (opt_atom->operator()(value)) return marley::JSON(value);
    // Return a null JSON object if the fhicl::OptionalAtom wasn't used
    else return marley::JSON();
  }

  template <typename S> marley::JSON fhicl_sequence_to_json(const S* sequence)
  {
    marley::JSON array = marley::JSON::make(
      marley::JSON::DataType::Array);

    for (size_t i = 0; i < sequence->size(); ++i) {
      array.append(sequence->operator()(i));
    }
    return array;
  }

  marley::JSON fhicl_parameter_to_json(
    const fhicl::detail::ParameterBase* par)
  {
    // Don't bother with the rest of the analysis if we've been handed a
    // nullptr or if the parameter is disabled in the current configuration
    if (par && par->should_use()) {

      auto type = par->parameter_type();
      if (fhicl::is_atom(type)) {

        // This is an ugly hack to load FHiCL atoms into JSON objects. We have
        // to do it type-by-type using the current (12/2016) implementation of
        // the fhicl::Atom and fhicl::OptionalAtom template classes.
        if (!par->is_optional()) {
          if (auto atm = dynamic_cast<const fhicl::Atom<double>* >(par))
            return fhicl_atom_to_json<>(atm);
          else if (auto atm = dynamic_cast<const
            fhicl::Atom<std::string>* >(par))
          {
            return fhicl_atom_to_json<>(atm);
          }
          else throw cet::exception("MarleyGen") << "Failed to deduce"
            << " the type of the FHiCL atom " << par->key();
        }
        else { // optional atoms
          if (auto opt_atom = dynamic_cast<const
            fhicl::OptionalAtom<double>* >(par))
          {
            return fhicl_optional_atom_to_json<>(opt_atom);
          }
          else if (auto opt_atom = dynamic_cast<const fhicl::OptionalAtom<
            std::string>* >(par))
          {
            return fhicl_optional_atom_to_json<>(opt_atom);
          }
          else throw cet::exception("MarleyGen") << "Failed to deduce"
            << " the type of the FHiCL optional atom " << par->key();
        }
      }
      else if (fhicl::is_sequence(type)) {

        // This is another ugly hack to allow us to fill the JSON array. We
        // have to do it type-by-type using the current (12/2016)
        // implementation of the fhicl::Sequence template class.
        if (auto seq = dynamic_cast<const fhicl::Sequence<double>* >(par))
          return fhicl_sequence_to_json<>(seq);
        else if (auto seq = dynamic_cast<const
          fhicl::Sequence<double, 3>* >(par))
        {
          return fhicl_sequence_to_json<>(seq);
        }
        else if (auto seq = dynamic_cast<const
          fhicl::Sequence<std::string>* >(par))
        {
          return fhicl_sequence_to_json<>(seq);
        }
        else throw cet::exception("MarleyGen") << "Failed to deduce"
          << " the type of the FHiCL sequence " << par->key();

        return marley::JSON();
      }
      else if (fhicl::is_table(type)) {

        // Downcast the parameter pointer into a TableBase pointer so that we
        // can iterate over its members
        auto table = dynamic_cast<const fhicl::detail::TableBase*>(par);
        if (!table) throw cet::exception("MarleyGen") << "Failed dynamic_cast"
          << " of fhicl::detail::ParameterBase* to fhicl::detail::TableBase*"
          << " when fhicl_is_table() was true";

        // Create an empty JSON object
        marley::JSON object = marley::JSON::make(
          marley::JSON::DataType::Object);

        // Load it with the members of the table recursively
        for (const auto& m : table->members()) {
          if (m && m->should_use()) {
            marley::JSON temp = fhicl_parameter_to_json(m.get());
	    // Skip members that are null (typically unused optional FHiCL
	    // parameters)
            if (!temp.is_null()) object[ m->name() ] = temp;
          }
        }
        return object;
      }
    }

    // Return a null JSON object if something strange happened.
    return marley::JSON();
  }

  /// Collection of configuration parameters used to
  /// determine the vertex location for each event
  struct Vertex_Config {
    fhicl::Atom<std::string> type_ {
      Name("type"),
      Comment("Technique used to choose MARLEY vertex locations"),
      "sampled" // default value
    };

    fhicl::OptionalAtom<std::string> seed_ {
      Name("seed"),
      Comment("Seed used for sampling MARLEY vertex locations"),
      [this]() -> bool { return type_() == "sampled"; }
    };

    fhicl::Sequence<double, 3> position_ {
      Name("position"),
      Comment("Coordinates of the fixed vertex position"),
      [this]() -> bool { return type_() == "fixed"; }
    };

  }; // struct Vertex_Config

  /// Collection of configuration parameters that will be
  /// forwarded to MARLEY and used to define the neutrino
  /// source
  struct Source_Config {
    fhicl::Atom<std::string> type_ {
      Name("type"),
      Comment("Type of neutrino source for MARLEY to use")
    };

    fhicl::Atom<std::string> neutrino_ {
      Name("neutrino"),
      Comment("Kind of neutrino (flavor, matter/antimatter) that the"
        " neutrino source produces")
    };

    fhicl::Atom<double> Emin_ {
      Name("Emin"),
      Comment("Minimum energy (MeV) of the neutrinos produced by this"
        " source"),
      [this]() -> bool {
        auto type = type_();
        return (type == "fermi-dirac") || (type == "beta-fit" );
      }
    };

    fhicl::Atom<double> Emax_ {
      Name("Emax"),
      Comment("Maximum energy (MeV) of the neutrinos produced by this"
        " source"),
      [this]() -> bool {
        auto type = type_();
        return (type == "fermi-dirac") || (type == "beta-fit" )
          || (type == "histogram");
      }
    };

    fhicl::Atom<double> temperature_ {
      Name("temperature"),
      Comment("Effective temperature for the Fermi-Dirac distribution"),
      [this]() -> bool {
        auto type = type_();
        return (type == "fermi-dirac");
      }
    };

    fhicl::OptionalAtom<double> eta_ {
      Name("eta"),
      Comment("Pinching parameter for the Fermi-Dirac distribution"),
      [this]() -> bool {
        auto type = type_();
        return (type == "fermi-dirac");
      }
    };

    fhicl::Atom<double> energy_ {
      Name("energy"),
      Comment("Energy (MeV) of the neutrinos produced by a monoenergetic"
        " source"),
      [this]() -> bool {
        auto type = type_();
        return (type == "monoenergetic");
      }
    };

    fhicl::Atom<double> Emean_ {
      Name("Emean"),
      Comment("Mean energy (MeV) of the neutrinos produced by a beta-fit"
        " source"),
      [this]() -> bool {
        auto type = type_();
        return (type == "beta-fit");
      }
    };

    fhicl::OptionalAtom<double> beta_ {
      Name("beta"),
      Comment("Pinching parameter for a beta-fit source"),
      [this]() -> bool {
        auto type = type_();
        return (type == "beta-fit");
      }
    };

    fhicl::Sequence<double> E_bin_lefts_ {
      Name("E_bin_lefts"),
      Comment("Left edges for each energy bin in the histogram"),
      [this]() -> bool {
        auto type = type_();
        return type == "histogram";
      }
    };

    fhicl::Sequence<double> weights_ {
      Name("weights"),
      Comment("Weights for each energy bin in the histogram"),
      [this]() -> bool {
        auto type = type_();
        return type == "histogram";
      }
    };

    fhicl::Sequence<double> energies_ {
      Name("energies"),
      Comment("Energies (MeV) for each grid point"),
      [this]() -> bool {
        auto type = type_();
        return type == "grid";
      }
    };

    fhicl::Sequence<double> prob_densities_ {
      Name("prob_densities"),
      Comment("Probability densities for each grid point"),
      [this]() -> bool {
        auto type = type_();
        return type == "grid";
      }
    };

    fhicl::Atom<std::string> rule_ {
      Name("rule"),
      Comment("Interpolation rule for computing probability densities"
        " between grid points. Allowed values include \"linlin\","
        " \"linlog\", \"loglin\", \"loglog\", and \"constant\""),
      [this]() -> bool {
        auto type = type_();
        return type == "grid";
      }
    };

    fhicl::Atom<std::string> tfile_ {
      Name("tfile"),
      Comment("Name of the ROOT file that contains a TH1 or TGraph neutrino"
        " source spectrum"),
      [this]() -> bool {
        auto type = type_();
        return (type == "th1") || (type == "tgraph");
      }
    };

    fhicl::Atom<std::string> namecycle_ {
      Name("namecycle"),
      Comment("Namecycle of the ROOT TH1 or TGraph to use"),
      [this]() -> bool {
        auto type = type_();
        return (type == "th1") || (type == "tgraph");
      }
    };

  }; // struct Source_Config

  /// Collection of configuration parameters that will be
  /// forwarded to MARLEY
  struct Marley_Config {

    fhicl::OptionalAtom<std::string> seed_ {
      Name("seed"),
      Comment("Seed used for the MARLEY generator")
    };

    fhicl::Sequence<double, 3> direction_ {
      Name("direction"),
      Comment("3-vector that points in the direction of motion of the"
        " incident neutrinos"),
      std::array<double, 3> { 0., 0., 1. } // default value
    };

    fhicl::Sequence<std::string> reactions_ {
      Name("reactions"),
      Comment("List of matrix element data files to use to define reactions"
        " in MARLEY")
    };

    fhicl::Sequence<std::string> structure_ {
      Name("structure"),
      Comment("List of TALYS format nuclear structure data files to use"
        " in MARLEY")
    };

    fhicl::Table<Source_Config> source_ {
      Name("source"),
      Comment("Neutrino source configuration")
    };

  }; // struct Marley_Config

  /// Collection of configuration parameters for the module
  struct Config {

    fhicl::Table<Vertex_Config> vertex_ {
      Name("vertex"),
      Comment("Configuration for selecting the vertex location(s)")
    };

    fhicl::Table<Marley_Config> marley_parameters_ {
      Name("marley_parameters"),
      Comment("Configuration for the MARLEY generator")
    };

    fhicl::Atom<std::string> module_type_ {
      Name("module_type"),
      Comment(),
      "MARLEYGen" // default value
    };

  }; // struct Config

}

namespace evgen {
  class MarleyGen;
}

class evgen::MarleyGen : public art::EDProducer {

public:

  // Type to enable FHiCL parameter validation by art
  using Parameters = art::EDProducer::Table<Config>;

  // Configuration-checking constructor
  explicit MarleyGen(const Parameters& p);
  virtual ~MarleyGen();

  virtual void produce(art::Event& e) override;
  virtual void beginRun(art::Run& run) override;

  virtual void reconfigure(const Parameters& p);
  void reconfigure_marley(const fhicl::Table<Marley_Config>& conf);
  void reconfigure_vertex(const fhicl::Table<Vertex_Config>& conf);

  void add_marley_particles(simb::MCTruth& truth,
    const std::vector<marley::Particle*>& particles, bool track);

  std::string find_file(const std::string& fileName,
    const std::string& fileType);

  void load_full_paths_into_json(marley::JSON& json,
    const std::string& array_name);

  enum class vertex_type_t { kSampled, kFixed };

protected:

  vertex_type_t fVertexType;

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
evgen::MarleyGen::MarleyGen(const Parameters& p)
  : fVertexType(vertex_type_t::kSampled), fEvent(new marley::Event),
  fRunNumber(0), fSubRunNumber(0), fEventNumber(0), fTPCDist(nullptr)
{
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
    "MARLEY", p().marley_parameters_.get_PSet(), { "seed" }
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
    "MarleyGenTPCEngine", p().vertex_.get_PSet(), { "seed" }
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
void evgen::MarleyGen::reconfigure(const Parameters& p)
{
  reconfigure_vertex(p().vertex_);
  reconfigure_marley(p().marley_parameters_);
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
      << fileName << '\'';

  return fullName;
}


//------------------------------------------------------------------------------
void evgen::MarleyGen::load_full_paths_into_json(
  marley::JSON& json, const std::string& array_name)
{
  if (json.has_key(array_name)) {
    // Replace each file name (which may appear in the FHiCL configuration
    // without a full path) with the full path found using cetlib
    for (auto& element : json.at(array_name).array_range()) {
      element = find_file(element.to_string(), array_name);
    }
  }
  else throw cet::exception("MarleyGen") << "Missing \"" << array_name
    << "\" key in the MARLEY parameters.";
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::reconfigure_marley(
  const fhicl::Table<Marley_Config>& conf)
{
  // Convert the FHiCL parameters into a JSON object that MARLEY can understand
  marley::JSON json = fhicl_parameter_to_json(&conf);

  // Update the reaction and structure data file names to the full paths
  // using cetlib to search for them
  load_full_paths_into_json(json, "reactions");
  load_full_paths_into_json(json, "structure");

  // Create a new MARLEY configuration based on the JSON parameters
  LOG_INFO("MarleyGen") << "MARLEY will now use the JSON configuration\n"
    << json.dump_string() << '\n';
  marley::JSONConfig config(json);

  // Create a new marley::Generator object basd on the current configuration
  fMarleyGenerator = std::make_unique<marley::Generator>(
    config.create_generator());
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
void evgen::MarleyGen::reconfigure_vertex(
  const fhicl::Table<Vertex_Config>& conf)
{
  auto type = conf().type_();
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

    auto vertex_pos = conf().position_();
    double Vx = vertex_pos.at(0);
    double Vy = vertex_pos.at(0);
    double Vz = vertex_pos.at(0);

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
