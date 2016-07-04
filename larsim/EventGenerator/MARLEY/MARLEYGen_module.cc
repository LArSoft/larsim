//////////////////////////////////////////////////////////////////////////////
/// \file MarleyGen_module.cc
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author sjgardiner@ucdavis.edu
//////////////////////////////////////////////////////////////////////////////

// standard library includes
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "larsim/RandomUtils/LArSeedService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Geant4 includes
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4IonTable.hh"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/Logger.hh"

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

protected:

  std::string fConfigFileName;
  std::string fMarleyLoggingLevel;

  std::unique_ptr<marley::Generator> fMarleyGenerator;

  std::stringstream fMarleyLogStream;

  // Function that samples a primary vertex position uniformly over the active
  // TPC volumes. The variables x, y, and z are loaded with the appropriate
  // values.
  void sample_vertex_pos(double& x, double& y, double& z);
  // Initializes sampling machinery for vertex position sampling
  void initialize_vertex_pos_sampling(const fhicl::ParameterSet& p);
  // RNG seed (currently only used for primary vertex position sampling)
  uint_fast64_t fSeed;
  // Discrete distribution object used to sample TPCs based on their active
  // masses
  std::unique_ptr<std::discrete_distribution<size_t> > fTPCDist;
  // RNG object used to sample TPCs
  std::mt19937_64 fTPCEngine;
};

//------------------------------------------------------------------------------
evgen::MarleyGen::MarleyGen(const fhicl::ParameterSet& p)
  : fSeed(0), fTPCDist(nullptr)
{
  // TODO: set MARLEY logging level from FHiCL parameters
  //fMarleyLoggingLevel = p.get<std::string>("MarleyLoggingLevel", "INFO");
  marley::Logger::Instance().add_stream(fMarleyLogStream,
    marley::Logger::LogLevel::INFO);

  this->reconfigure(p);

  // Log initialization information from the MARLEY generator
  LOG_INFO("MARLEY") << fMarleyLogStream.str();
  fMarleyLogStream.str("");
  fMarleyLogStream.clear();

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

  std::unique_ptr< std::vector<simb::MCTruth> >
    truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;

  static double vertex_x = 0.;
  static double vertex_y = 0.;
  static double vertex_z = 0.;

  // We need to convert from MARLEY's energy units (MeV) to LArSoft's
  // (GeV) using this conversion factor
  static constexpr double MeV_to_GeV = 1e-3;

  // Sample a primary vertex location for this event
  sample_vertex_pos(vertex_x, vertex_y, vertex_z);

  // TODO: include sampling of vertex times
  TLorentzVector pos(vertex_x, vertex_y, vertex_z, 0.);

  marley::Event ev = fMarleyGenerator->create_event();

  size_t fp_count = 0;
  for (const marley::Particle* fp : ev.get_final_particles()) {

    int pdg = fp->pdg_code();
    double mass = fp->mass() * MeV_to_GeV;
    double px = fp->px() * MeV_to_GeV;
    double py = fp->py() * MeV_to_GeV;
    double pz = fp->pz() * MeV_to_GeV;
    double E = fp->total_energy() * MeV_to_GeV;
    TLorentzVector mom(px, py, pz, E);

    simb::MCParticle part(fp_count, pdg, "MARLEY", -1 /* primary particle */,
      mass, 1 /* track this particle in LArG4 */);
    part.AddTrajectoryPoint(pos, mom);
    truth.Add(part);

    ++fp_count;

    // If Geant4 has not yet produced a particle definition for a nucleus
    // generated by MARLEY, then prepare one. The LArG4 module will refuse to
    // track it unless we make a definition in advance.
    if (pdg > 1000000000 && // nucleus
      !G4ParticleTable::GetParticleTable()->FindParticle(pdg))
    {
      int Z = (pdg % 10000000) / 10000;
      int A = (pdg % 10000) / 10;
      G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()
        ->GetIonTable()->GetIon(Z, A, 0.);
      if (ion) LOG_INFO("MarleyGen") << "Added ion Z = " << Z << ", A = " << A;
      else LOG_INFO("MarleyGen") << "Failed to create ion Z = " << Z
        << ", A = " << A;
    }
  }

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
  fConfigFileName = p.get<std::string>("ConfigFileName");
  fMarleyLoggingLevel = p.get<std::string>("MarleyLoggingLevel", "INFO");

  fMarleyGenerator = std::make_unique<marley::Generator>(fConfigFileName);

  initialize_vertex_pos_sampling(p);
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::sample_vertex_pos(double& x, double& y, double& z)
{
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
  x = uniform_dist(fTPCEngine, x_range);
  y = uniform_dist(fTPCEngine, y_range);
  z = uniform_dist(fTPCEngine, z_range);
  LOG_INFO("MarleyGen") << "Sampled primary vertex in TPC #"
    << tpc_index << ", x = " << x << ", y = " << y << ", z = " << z;
}

//------------------------------------------------------------------------------
void evgen::MarleyGen::initialize_vertex_pos_sampling(const fhicl::ParameterSet& p)
{
  // Get the RNG seed from the configuration file (or, if it's not set, use a
  // seed from the LArSeedService)
  fSeed = p.get<uint_fast64_t>("Seed",
    art::ServiceHandle< sim::LArSeedService >()->getSeed());

  // Use the obtained seed to prepare the random number engine.  This is an
  // attempt to do a decent job, but optimally accomplishing this can be tricky
  // (see, for example,
  // http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
  std::seed_seq seed_sequence{fSeed};
  fTPCEngine.seed(seed_sequence);

  // Get the active masses of each of the TPCs in the current geometry. Use
  // them as weights for sampling a TPC to use for the primary vertex.
  art::ServiceHandle<geo::Geometry> geom;
  std::vector<double> tpc_masses;
  size_t num_tpcs = geom->NTPC();
  for (size_t iTPC = 0; iTPC < num_tpcs; ++iTPC) {
    // For each TPC, use its active mass (returned in kg) as its sampling weight
    tpc_masses.push_back(geom->TPC(iTPC).ActiveMass());
  }

  // Load the discrete distribution used to sample TPCs with the up-to-date set
  // of weights
  fTPCDist.reset(new std::discrete_distribution<size_t>(tpc_masses.begin(),
    tpc_masses.end()));
}

DEFINE_ART_MODULE(evgen::MarleyGen)
