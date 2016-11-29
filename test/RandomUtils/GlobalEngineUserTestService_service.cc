/**
 * @file    GlobalEngineUserService_service.cc
 * @brief   Test service registering its own engine (implementation file)
 * @author  Gianluca Petrillo
 * @date    March 22, 2016
 * @see     GlobalEngineUserTestService.h
 * 
 */


// my header
#include "test/RandomUtils/GlobalEngineUserTestService.h"

// framework libraries
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// supporting libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TClass.h"
#include "TRandom3.h"


//------------------------------------------------------------------------------
//--- testing::GlobalEngineUserTestService
//---
const std::string testing::GlobalEngineUserTestService::GlobalInstanceName
  { "GlobalEngineUserServiceEngine" };

testing::GlobalEngineUserTestService::GlobalEngineUserTestService
  (fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  : instanceNames       (pset.get<std::vector<std::string>>("instances", {}))
  , perEventSeeds       (pset.get<bool>                    ("perEventSeeds", false))
{
  
  //
  // input parameters check
  //
  if (instanceNames.empty()) {
    instanceNames.push_back("GlobalEngineUserTestService");
    mf::LogInfo("GlobalEngineUserTestService")
      << "Using a default engine instance name: '" << instanceNames.back()
      << "'";
  }
  
  //
  // create all our engines (uninteresting temporary seeds are set)
  //
  for (std::string instanceName: instanceNames) {
    engines.emplace_back(std::make_unique<TRandom3>());
    engines.back()->SetTitle(instanceName.c_str());
  } // for

  //
  // register all the engines
  //
  rndm::NuRandomService& Seeds = *(art::ServiceHandle<rndm::NuRandomService>());
  // NuRandomService::createEngine() can't be called here
  //   because it needs a EngineCreator: good
  
  for (auto& engine: engines) {
    
    // NuRandomService::registerEngine() should instead succeed;
    // rndm::NuRandomService::TRandomSeeder is optionally declared (inline)
    // in NuRandomService.h
    auto seed = Seeds.registerEngine
      (rndm::NuRandomService::TRandomSeeder(engine.get()), engine->GetTitle());
    
    mf::LogInfo("GlobalEngineUserTestService")
      << "Registered my random engine "
      << engine->IsA()->GetName() << "[" << ((void*) engine.get()) << "]"
      << " with seed " << seed;
    
  } // for

  //
  // let's check that the seeds were actually set
  //
  if (perEventSeeds) {
    mf::LogWarning("GlobalEngineUserTestService")
      << "Check of seeds on construction skipped because policy is per event.";
  }
  else CheckAllSeeds();
  
  //
  // register callbacks
  //
  reg.sPreModuleConstruction.watch (this, &GlobalEngineUserTestService::preModuleConstruction  );
  reg.sPostModuleConstruction.watch(this, &GlobalEngineUserTestService::postModuleConstruction );
  reg.sPreModuleBeginRun.watch     (this, &GlobalEngineUserTestService::preModuleBeginRun      );
  reg.sPostModuleBeginRun.watch    (this, &GlobalEngineUserTestService::postModuleBeginRun     );
  reg.sPreProcessEvent.watch       (this, &GlobalEngineUserTestService::preProcessEvent        );
  reg.sPreModule.watch             (this, &GlobalEngineUserTestService::preModule              );
  reg.sPostModule.watch            (this, &GlobalEngineUserTestService::postModule             );
  reg.sPostProcessEvent.watch      (this, &GlobalEngineUserTestService::postProcessEvent       );
  reg.sPreModuleEndJob.watch       (this, &GlobalEngineUserTestService::preModuleEndJob        );
  reg.sPostModuleEndJob.watch      (this, &GlobalEngineUserTestService::postModuleEndJob       );
  
} // testing::GlobalEngineUserTestService::GlobalEngineUserTestService
   

//------------------------------------------------------------------------------

void testing::GlobalEngineUserTestService::CheckSeed(TRandom const& engine) {
  
  auto& Seeds = *(art::ServiceHandle<rndm::NuRandomService>());
  
  auto expectedSeed = Seeds.getGlobalCurrentSeed(engine.GetTitle());
  auto actualSeed = engine.GetSeed();
  
  if (actualSeed != expectedSeed) {
    mf::LogError("GlobalEngineUserTestService")
      << "Engine " << engine.IsA()->GetName() << "[" << ((void*) &engine)
      << "](\"" << engine.GetTitle() << "\") has seed " << actualSeed
      << " (expected: " << expectedSeed << ")";
  }
  else {
    mf::LogVerbatim("GlobalEngineUserTestService")
      << "Engine " << engine.IsA()->GetName() << "[" << ((void*) &engine)
      << "](\"" << engine.GetTitle() << "\") has seed " << actualSeed
      << " (as expected)";
  }
  if (actualSeed != expectedSeed) {
    throw art::Exception(art::errors::LogicError)
      << "Unexpected seed " << actualSeed << "!\n";
  }
} // testing::GlobalEngineUserTestService::CheckSeed()


void testing::GlobalEngineUserTestService::CheckAllSeeds() const {
  mf::LogVerbatim("GlobalEngineUserTestService")
    << "Checking all " << engines.size() << " seeds";
  for (auto& engine: engines) CheckSeed(*engine);
} // testing::GlobalEngineUserTestService::CheckAllSeeds()


//------------------------------------------------------------------------------

void testing::GlobalEngineUserTestService::preModuleConstruction
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::preModuleConstruction()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::preModuleConstruction()

     
void testing::GlobalEngineUserTestService::postModuleConstruction
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::postModuleConstruction()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::postModuleConstruction()

     
void testing::GlobalEngineUserTestService::preModuleBeginRun
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::preModuleBeginRun()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::preModuleBeginRun()

     
void testing::GlobalEngineUserTestService::postModuleBeginRun
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::postModuleBeginRun()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::postModuleBeginRun()

     
void testing::GlobalEngineUserTestService::preProcessEvent(art::Event const&) {
   LOG_DEBUG("GlobalEngineUserTestService")
     << "GlobalEngineUserTestService::preProcessEvent()";
   CheckAllSeeds();
} // testing::GlobalEngineUserTestService::preProcessEvent()

     
void testing::GlobalEngineUserTestService::preModule
  (art::ModuleDescription const&)
{
   LOG_DEBUG("GlobalEngineUserTestService")
     << "GlobalEngineUserTestService::preModule()";
   CheckAllSeeds();
} // testing::GlobalEngineUserTestService::preModule()

     
void testing::GlobalEngineUserTestService::postModule
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::postModule()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::postModule()

     
void testing::GlobalEngineUserTestService::postProcessEvent(art::Event const&) {
   LOG_DEBUG("GlobalEngineUserTestService")
     << "GlobalEngineUserTestService::postProcessEvent()";
   CheckAllSeeds();
} // testing::GlobalEngineUserTestService::postProcessEvent()

     
void testing::GlobalEngineUserTestService::preModuleEndJob
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::preModuleEndJob()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::preModuleBeginRun()

     
void testing::GlobalEngineUserTestService::postModuleEndJob
  (art::ModuleDescription const&)
{
   if (!perEventSeeds) {
     LOG_DEBUG("GlobalEngineUserTestService")
       << "GlobalEngineUserTestService::postModuleEndJob()";
     CheckAllSeeds();
   }
} // testing::GlobalEngineUserTestService::preModuleBeginRun()

     
//------------------------------------------------------------------------------

DEFINE_ART_SERVICE(testing::GlobalEngineUserTestService)

//------------------------------------------------------------------------------
