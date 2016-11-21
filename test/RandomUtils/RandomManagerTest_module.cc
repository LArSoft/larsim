/**
 * @file   RandomManagerTest_module.cc
 * @brief  Test of the random engine managing interface of NuRandomService
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 19th, 2015
 */


// art extensions
#define NUTOOLS_RANDOMUTILS_NuRandomService_USECLHEP 1 // to have NuRandomService.h define CLHEPengineSeeder
#include "nutools/RandomUtils/NuRandomService.h"

#include "test/RandomUtils/SeedTestUtils.h"

// C++ includes.
#include <string>
#include <array>
#include <vector>
#include <sstream>
#include <algorithm> // std::generate()
#include <iomanip> // std::setw()

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine
#include "CLHEP/Random/Ranlux64Engine.h" // CLHEP::Ranlux64Engine
#include "CLHEP/Random/RandFlat.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// Framework includes.
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/CurrentModule.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace testing {
  
  /**
   * @brief Test module for random engine managing interface of NuRandomService
   * 
   * The test writes on screen the random seeds it gets.
   * 
   * Configuration parameters:
   * - *instanceNames* (string list, optional): use one random number
   *   generator for each instance name here specified; if not specified,
   *   an anonymous engine is used
   * - *externalInstance* (string, optional): if specified, an engine not
   *   managed by RandomNumberGenerator is also used, with this instance name
   * - *standardInstance* (string, optional): if specified, an engine
   *   is created by RandomNumberGenerator but not registered is NuRandomService
   *   is also used, with this instance name
   * - *Seed*, *Seed_XXX* (strings, optional): set the seed of instance `XXX`
   *   to a set value ("Seed" sets the seed of the anonymous instance)
   * 
   */
  class RandomManagerTest: public art::EDAnalyzer {
      public:
    typedef art::RandomNumberGenerator::seed_t seed_t;
    
    explicit RandomManagerTest(fhicl::ParameterSet const& pset);
    
    void analyze(const art::Event& event) override;
    
      private:
    std::vector<std::string> instanceNames;
    std::string externalInstanceName;
    std::string standardInstanceName;
    
    std::unique_ptr<CLHEP::HepRandomEngine> extEngine;
    CLHEP::HepRandomEngine* stdEngine;
    
    std::string moduleName() const;
    
  }; // class RandomManagerTest
  
  
  
  //****************************************************************************
  //--- RandomManagerTest implementation
  //---
  RandomManagerTest::RandomManagerTest(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    instanceNames(),
    stdEngine(nullptr)
  {
    
    art::ServiceHandle<rndm::NuRandomService> EngineManager;
    
    // check if we want an "external" engine
    if (pset.get_if_present("externalInstance", externalInstanceName)) {
      mf::LogInfo("RandomManagerTest") << "Creating an unmanaged engine '"
        << externalInstanceName << "' in module '" << moduleName() << "'";
      extEngine.reset(new CLHEP::Ranlux64Engine);
      
      EngineManager->registerEngine(
      //  rndm::NuRandomService::CLHEPengineSeeder(extEngine.get()),
        [this](rndm::NuRandomService::EngineId const&, seed_t seed)
          { this->extEngine->setSeed(seed, 0); },
        externalInstanceName, pset, "Seed_" + externalInstanceName
        );
    } // if we have the external engine
    
    // check if we want an unmanaged standard engine
    if (pset.get_if_present("standardInstance", standardInstanceName)) {
      mf::LogInfo("RandomManagerTest") << "Creating a standard engine '"
        << standardInstanceName << "' in module '" << moduleName()
        << "' with RandomNumberGenerator";
      seed_t seed
        = pset.get<unsigned int>("Seed_" + standardInstanceName, 0);
      createEngine(seed, "HepJamesRandom", standardInstanceName);
      stdEngine = &(art::ServiceHandle<art::RandomNumberGenerator>()
        ->getEngine(standardInstanceName));
    } // if we have the external engine
    
    // initialize the standard engines with RandomNumberGenerator
    if (pset.get_if_present("instanceNames", instanceNames)) {
      for (std::string const& instanceName: instanceNames) {
        mf::LogInfo("RandomManagerTest") << "Creating a default engine '"
          << instanceName << "' in module '" << moduleName() << "'";
        EngineManager->createEngine
          (*this, "HepJamesRandom", instanceName, pset, "Seed_" + instanceName);
      } // for
    }
    
    // create a default engine, if needed
    if (instanceNames.empty() && !extEngine && !stdEngine) {
      instanceNames.emplace_back();
      mf::LogInfo("RandomManagerTest")
        << "Creating a nameless default engine in module '"
        << moduleName() << "'";
      EngineManager->createEngine(*this, pset, "Seed");
    }
    
    { // anonymous block
      mf::LogInfo log("RandomManagerTest");
      log << "RandomManagerTest[" << moduleName() << "]: instances:";
      for (std::string const& instanceName: instanceNames)
        log << " " << instanceName;
    } // anonymous block
    
  } // RandomManagerTest::RandomManagerTest()
  
  
  //----------------------------------------------------------------------------
  void RandomManagerTest::analyze(const art::Event& event){
    { // anonymous block
      mf::LogVerbatim log("RandomManagerTest");
      log << "RandomManagerTest[" << moduleName() << "]::analyze "
          << event.id();
    } // anonymous block
    
    // prepare a list of all the engines (including the external one)
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    
    std::vector<std::string> allInstances(instanceNames);
    std::vector<CLHEP::HepRandomEngine*> allEngines;
    for (std::string const& instanceName: instanceNames)
      allEngines.push_back(&rng->getEngine(instanceName));
    
    if (extEngine) {
      allInstances.push_back(externalInstanceName);
      allEngines.push_back(extEngine.get());
    }
    
    if (stdEngine) {
      allInstances.push_back(standardInstanceName);
      allEngines.push_back(stdEngine);
    }
    
    for (size_t iEngine = 0; iEngine < allEngines.size(); ++iEngine) {
      std::string const& instanceName = allInstances[iEngine];
      CLHEP::HepRandomEngine& engine = *(allEngines[iEngine]);
      
      seed_t actualSeed = testing::NuRandomService::readSeed(engine);
      mf::LogVerbatim("RandomManagerTest")
        << std::setw(12) << (instanceName.empty()? "<default>": instanceName)
        << ": " << testing::NuRandomService::CreateCharacter(engine)
        << "   (seed: " << actualSeed << ")";
    } // for
    
  } // RandomManagerTest::analyze()
  
  
  inline std::string RandomManagerTest::moduleName() const
    { return art::ServiceHandle<art::CurrentModule>()->label(); }
  
  
} // end namespace testing

DEFINE_ART_MODULE(testing::RandomManagerTest)
