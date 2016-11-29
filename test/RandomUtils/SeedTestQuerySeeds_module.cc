/**
 * @file   SeedTestQuerySeeds_module.cc
 * @brief  Tests query of seeds in various contexts.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 18, 2016
 */


// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


// forward declarations
namespace art {
  class Run;
  class SubRun;
  class Event;
} // namespace art


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
  class SeedTestQuerySeeds: public art::EDAnalyzer {
    
      public:

    explicit SeedTestQuerySeeds(fhicl::ParameterSet const& pset);

    virtual void analyze(art::Event const&) override;

    virtual void beginRun   (art::Run const&) override;
    virtual void beginSubRun(art::SubRun const&) override;
    
    virtual void endJob() override;
    
  }; // class SeedTestQuerySeeds
  
  
  //----------------------------------------------------------------------------
  SeedTestQuerySeeds::SeedTestQuerySeeds
    (fhicl::ParameterSet const& pSet)
    : art::EDAnalyzer(pSet)
  {
    // This call has two effects:
    // 1. declares the existence of a random engine, and associates it with
    //    a default (empty) instance name
    //    [we don't have such engine, but that does not matter for the test]
    // 2. gets the seed associated with that engine
    
    auto seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
    mf::LogInfo("SeedTestQuerySeeds")
      << "Declared a default engine, got seed=" << seed;
  }
  
  
  //----------------------------------------------------------------------------
  void SeedTestQuerySeeds::analyze(art::Event const&) {
    
    auto seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
    mf::LogInfo("SeedTestQuerySeeds")
      << "In event loop: got seed=" << seed << " for the default engine";
    
  } // SeedTestQuerySeeds::analyze()
  
  
  void SeedTestQuerySeeds::beginRun(art::Run const&) {
    
    auto seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
    mf::LogInfo("SeedTestQuerySeeds")
      << "At begin run: got seed=" << seed << " for the default engine";
    
  } // SeedTestQuerySeeds::beginRun()
  
  
  void SeedTestQuerySeeds::beginSubRun(art::SubRun const&) {
    
    auto seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
    mf::LogInfo("SeedTestQuerySeeds")
      << "At begin subrun: got seed=" << seed << " for the default engine";
    
  } // SeedTestQuerySeeds::beginSubRun()
  
  
  void SeedTestQuerySeeds::endJob() {
    
    auto seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
    mf::LogInfo("SeedTestQuerySeeds")
      << "At end job: got seed=" << seed << " for the default engine";
    
  } // SeedTestQuerySeeds::endJob()
  
  
  //----------------------------------------------------------------------------
  
} // namespace testing

DEFINE_ART_MODULE(testing::SeedTestQuerySeeds)
