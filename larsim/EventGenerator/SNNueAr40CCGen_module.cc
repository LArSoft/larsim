////////////////////////////////////////////////////////////////////////////////
/// \file  SNNueAr40CCGen_module.cc
/// \brief Supernova Nue-Ar CC interaction generator
/// 
/// Module that produces MC particles from supernova electron neutrinos
/// interacting via charge current with Ar40
/// 
/// \author  gleb.sinev@duke.edu
////////////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"

// LArSoft includes
#include "larcoreobj/SummaryData/RunData.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventGenerator/NueAr40CCGenerator.h"

#include "CLHEP/Random/RandomEngine.h"

// C++ includes
#include <map>

namespace evgen {
  
  /// Module producing electrons and gammas from supernova neutrino 
  /// interactions with liquid argon (based on SingleGen)
  class SNNueAr40CCGen : public art::EDProducer {
    public:
      explicit SNNueAr40CCGen(fhicl::ParameterSet const& pset);

  private:
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;

    CLHEP::HepRandomEngine& fEngine;
    evgen::NueAr40CCGenerator fGenerator;
  };

  //____________________________________________________________________________
  SNNueAr40CCGen::SNNueAr40CCGen(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    // Create a default random engine: obtain the random seed
    // freom NuRandomService, unless overriden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, pset, "Seed"))
    , fGenerator{pset.get< fhicl::ParameterSet >("GeneratorAlg")}
  {
    produces< std::vector< simb::MCTruth > >();
    produces< sumdata::RunData, art::InRun >();
  }

  //____________________________________________________________________________
  void SNNueAr40CCGen::beginRun(art::Run& run)
  {
    // Store information about the geometry we are using in run information
    art::ServiceHandle< geo::Geometry > geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
  }

  //____________________________________________________________________________
  void SNNueAr40CCGen::produce(art::Event& event)
  {
    event.put(std::make_unique<std::vector<simb::MCTruth>>(fGenerator.Generate(fEngine)));
  }

}

DEFINE_ART_MODULE(evgen::SNNueAr40CCGen)
