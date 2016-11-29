/**
 * @file   SeedTestRegisterOnProcessEvent_module.cc
 * @brief  Tests engine registration and query of seeds in analyze() context.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 17, 2016
 */


// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/Exception.h"

// supporting libraries
#include "messagefacility/MessageLogger/MessageLogger.h"


// forward declarations
namespace art { class Event; }


namespace testing {
  
  /**
   * @brief Test module for NuRandomService
   * 
   * The test tries to get a seed of an unknown engine during event processing.
   * 
   * It is expected to fail, since the engine should be registered on
   * construction.
   * 
   * Note that the test does not actually get any random number.
   * 
   * Configuration parameters: none
   */
  class SeedTestRegisterOnProcessEvent: public art::EDAnalyzer {
    
      public:

    explicit SeedTestRegisterOnProcessEvent(fhicl::ParameterSet const& pset);

    virtual void analyze(art::Event const&) override;

  }; // class SeedTestRegisterOnProcessEvent
  
  
  SeedTestRegisterOnProcessEvent::SeedTestRegisterOnProcessEvent
    (fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset)
    {}

  void SeedTestRegisterOnProcessEvent::analyze(art::Event const&) {
    
    art::ServiceHandle<rndm::NuRandomService> Seeds;
    auto const seed = Seeds->getSeed(); // should throw
    
    mf::LogError("SeedTestRegisterOnProcessEvent")
      << "Got seed=" << seed
      << " for unregistered (unnamed) engine in analyze(); should have thrown!";
    
  } // SeedTestRegisterOnProcessEvent::analyze()

} // end namespace testing

DEFINE_ART_MODULE(testing::SeedTestRegisterOnProcessEvent)
