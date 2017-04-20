//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYGen_module.cc
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author sjgardiner@ucdavis.edu
//////////////////////////////////////////////////////////////////////////////

// standard library includes
#include <memory>
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
#include "fhiclcpp/types/Table.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/EventGenerator/MARLEY/MARLEYGenerator.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// anonymous namespace for helper functions, etc.
namespace {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

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

  /// Collection of configuration parameters for the module
  struct Config {

    fhicl::Table<Vertex_Config> vertex_ {
      Name("vertex"),
      Comment("Configuration for selecting the vertex location(s)")
    };

    fhicl::Table<evgen::MARLEYGenerator::Config> marley_parameters_ {
      Name("marley_parameters"),
      Comment("Configuration for the MARLEY generator")
    };

    fhicl::Atom<std::string> module_type_ {
      Name("module_type"),
      Comment(""),
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

  // Configuration-checking constructors
  explicit MarleyGen(const Parameters& p);

  virtual ~MarleyGen();

  virtual void produce(art::Event& e) override;
  virtual void beginRun(art::Run& run) override;

  virtual void reconfigure(const Parameters& p);

  void reconfigure_marley(
    const fhicl::Table<evgen::MARLEYGenerator::Config>& conf);

  void reconfigure_vertex(const fhicl::Table<Vertex_Config>& conf);

  enum class vertex_type_t { kSampled, kFixed };

protected:

  // Function that selects a primary vertex location for each event. The
  // result is saved to fVertexPosition.
  void update_vertex_pos();

  vertex_type_t fVertexType;

  TLorentzVector fVertexPosition;

  std::unique_ptr<evgen::MARLEYGenerator> fMarleyGenerator;

  // unique_ptr to the current event created by MARLEY
  std::unique_ptr<marley::Event> fEvent;

  // the MARLEY event TTree
  TTree* fEventTree;

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

  // Register the TPC sampling engine with the seed service. Iff you need the
  // seed later, get it from the seed service using "MarleyGenTPCEngine" as the
  // instance name.
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

  // TODO: resolve the other workaround mentioned in the MARLEYGenerator
  // class, then fix this as well
  uint_fast64_t tpc_cast_seed = static_cast<uint_fast64_t>(tpc_seed);
  std::seed_seq tpc_seed_sequence{tpc_cast_seed};
  fTPCEngine.seed(tpc_seed_sequence);

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
void evgen::MarleyGen::produce(art::Event& e)
{

  // Get the run, subrun, and event numbers from the current art::Event
  fRunNumber = e.run();
  fSubRunNumber = e.subRun();
  fEventNumber = e.event();

  std::unique_ptr< std::vector<simb::MCTruth> >
    truthcol(new std::vector<simb::MCTruth>);

  // Get the primary vertex location for this event
  update_vertex_pos();

  // Create the MCTruth object, and retrieve the marley::Event object
  // that was generated as it was created
  simb::MCTruth truth = fMarleyGenerator->create_MCTruth(fVertexPosition,
    fEvent.get());

  // Write the marley::Event object to the event tree
  fEventTree->Fill();

  truthcol->push_back(truth);

  e.put(std::move(truthcol));
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::reconfigure(const Parameters& p)
{
  reconfigure_vertex(p().vertex_);
  reconfigure_marley(p().marley_parameters_);
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::reconfigure_marley(
  const fhicl::Table<evgen::MARLEYGenerator::Config>& conf)
{
  const auto& seed_service = art::ServiceHandle<rndm::NuRandomService>();

  // Create a new marley::Generator object basd on the current configuration
  fMarleyGenerator = std::make_unique<MARLEYGenerator>(conf, *seed_service,
    "MARLEY");
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

DEFINE_ART_MODULE(evgen::MarleyGen)
