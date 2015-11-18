////////////////////////////////////////////////////////////////////////////////
/// \file  SNNueAr40CCGen_module.cc
/// \brief Supernova Nue-Ar CC interaction generator
/// 
/// Module that produces MC particles from supernova electron neutrinos
/// interacting via charge current with Ar40
/// 
/// \version $Id: SingleGen_module.cc, v 1.0 2015/11/16 17:28:00 gvsinev Exp $
/// \author  gleb.sinev@duke.edu
////////////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_SNNUEARCCGEN
#define EVGEN_SNNUEARCCGEN

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/ParameterSet.h"

// art extensions
#include "artextensions/SeedService/SeedService.hh"

// nutools includes
#include "SimulationBase/MCTruth.h"

// LArSoft includes
#include "SummaryData/RunData.h"
#include "Geometry/Geometry.h"

#include "CLHEP/Random/RandomEngine.h"

// C++ includes
#include <map>

#include "NueAr40CCGenerator.h"

namespace evgen {
  
  /// Module producing electrons and gammas from supernova neutrino 
  /// interactions with liquid argon (based on SingleGen)
  class SNNueAr40CCGen : public art::EDProducer {

    public:

      explicit SNNueAr40CCGen(fhicl::ParameterSet const& pset);

      void beginRun(art::Run& run);

      void produce(art::Event& event);

    private:

    evgen::NueAr40CCGenerator fGenerator;

  };

}

namespace evgen {

  //____________________________________________________________________________
  SNNueAr40CCGen::SNNueAr40CCGen(fhicl::ParameterSet const& pset)
    : fGenerator(evgen::NueAr40CCGenerator
                             (pset.get< fhicl::ParameterSet >("GeneratorAlg")))
  {

    produces< std::vector< simb::MCTruth > >();
    produces< sumdata::RunData, art::InRun >();

    // Create a default random engine: obtain the random seed
    // freom SeedService, unless overriden in configuration with key "Seed"
    art::ServiceHandle< artext::SeedService >()
      ->createEngine(*this, pset, "Seed");

  }

  //____________________________________________________________________________
  void SNNueAr40CCGen::beginRun(art::Run& run)
  {

    // Store information about the geometry we are using in run information
    art::ServiceHandle< geo::Geometry > geo;
    std::unique_ptr< sumdata::RunData > 
      runCol(new sumdata::RunData(geo->DetectorName()));

    run.put(std::move(runCol));

    return;

  }

  //____________________________________________________________________________
  void SNNueAr40CCGen::produce(art::Event& event)
  {

    std::unique_ptr< std::vector< simb::MCTruth > > 
                                truthCol(new std::vector< simb::MCTruth >);

    // Get an engine from the random number generator
    art::ServiceHandle< art::RandomNumberGenerator > randomNumberGenerator;
    CLHEP::HepRandomEngine &engine = randomNumberGenerator->getEngine();

    truthCol->emplace_back(fGenerator.Generate(engine));

    event.put(std::move(truthCol));

    return;

  }

}

namespace evgen {

  DEFINE_ART_MODULE(SNNueAr40CCGen)

}

#endif
////////////////////////////////////////////////////////////////////////////////
