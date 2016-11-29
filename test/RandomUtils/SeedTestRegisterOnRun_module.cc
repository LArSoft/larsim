/**
 * @file   SeedTestRegisterOnRun_module.cc
 * @brief  Tests engine registration and query of seeds in beginRun() context.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 17, 2016
 */


// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// Framework includes.
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// forward declarations
namespace art { class SubRun; }


namespace testing {
  
  /**
   * @brief Test module for NuRandomService
   * 
   * The test tries to get a seed of an unknown engine during begin subrun.
   * 
   * Note that the test does not actually get any random number.
   * 
   * Configuration parameters
   * -------------------------
   * 
   * * *register* (bool, default: true): register the engine
   * 
   */
  class SeedTestRegisterOnRun: public art::EDAnalyzer {
    
      public:

    explicit SeedTestRegisterOnRun(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event&) override {}

    virtual void beginRun(art::Run const&) override;
    
      private:
    
    bool registerEngine; ///< whether to register the engine on construction
    
  }; // class SeedTestRegisterOnRun
  
  
  SeedTestRegisterOnRun::SeedTestRegisterOnRun
    (fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset)
    , registerEngine(pset.get<bool>("register", true))
  {
    if (registerEngine) {
      art::ServiceHandle<rndm::NuRandomService> Seeds;
      Seeds->declareEngine(); // declare an engine with empty instance name
    }
  } // SeedTestRegisterOnRun::SeedTestRegisterOnRun()
  
  
  void SeedTestRegisterOnRun::beginRun(art::Run const&) {
    
    art::ServiceHandle<rndm::NuRandomService> Seeds;
    auto const seed = Seeds->getSeed(); // should throw
    
    mf::LogInfo("SeedTestRegisterOnRun")
      << "Got seed=" << seed
      << " for " << (registerEngine? "registered": "unregistered")
      << " (unnamed) engine in beginRun()";
    
  } // SeedTestRegisterOnRun::beginSubRun()
  
} // namespace testing

DEFINE_ART_MODULE(testing::SeedTestRegisterOnRun)
