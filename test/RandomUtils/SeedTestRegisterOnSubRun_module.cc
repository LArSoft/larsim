/**
 * @file   SeedTestRegisterOnSubRun_module.cc
 * @brief  Tests engine registration and query of seeds in beginSubRun() context.
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
   * It is expected to fail, since the engine should be registered on
   * construction.
   * 
   * Note that the test does not actually get any random number.
   * 
   * Configuration parameters: none
   */
  class SeedTestRegisterOnSubRun: public art::EDAnalyzer {
    
      public:

    explicit SeedTestRegisterOnSubRun(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event&) override {}

//    void beginRun   ( const art::Run&    run) override;
    virtual void beginSubRun(art::SubRun const&) override;
    
  }; // class SeedTestRegisterOnSubRun
  
  
  SeedTestRegisterOnSubRun::SeedTestRegisterOnSubRun
    (fhicl::ParameterSet const& pSet)
    : art::EDAnalyzer(pSet)
    {}
  
  
  void SeedTestRegisterOnSubRun::beginSubRun(art::SubRun const&) {
    
    art::ServiceHandle<rndm::NuRandomService> Seeds;
    auto const seed = Seeds->getSeed(); // should throw
    
    mf::LogError("SeedTestRegisterOnSubRun")
      << "Got seed=" << seed
      << " for unregistered (unnamed) engine in beginSubRun(); should have thrown!";
    
  } // SeedTestRegisterOnSubRun::beginSubRun()
  
} // namespace testing

DEFINE_ART_MODULE(testing::SeedTestRegisterOnSubRun)
