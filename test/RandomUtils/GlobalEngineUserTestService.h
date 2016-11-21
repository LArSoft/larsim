/**
 * @file    GlobalEngineUserTestService.h
 * @brief   Test service registering its own engine
 * @author  Gianluca Petrillo
 * @see     GlobalEngineUserService_service.cc
 * 
 */

#ifndef ARTEXT_SEEDSERVICE_TEST_GLOBALENGINEUSERSERVICE_H
#define ARTEXT_SEEDSERVICE_TEST_GLOBALENGINEUSERSERVICE_H


// NuRandomService
#define NUTOOLS_RANDOMUTILS_NuRandomService_USEROOT 1 // get testing::TRandomSeeder (it's a test...)
#include "nutools/RandomUtils/NuRandomService.h"

// From art and its tool chain.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

// ROOT libraries
#include "TRandom.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <memory> // std::unique_ptr<>


// Forward declarations
namespace art {
  class ActivityRegistry;
  class ModuleDescription;
}

namespace testing {
   
   /**
    * @brief Test service registering its own ßs
    * 
    * The service owns all its engines and is asks NuRandomService for their seeds.
    * NuRandomService will use global identifiers for them, that is identifiers not
    * bound to the art execution context (that is, not bound to a specific
    * module).
    * 
    * 
    * Configuration parameters
    * -------------------------
    * 
    * * *instances* (list of strings, default: empty): list of names of the
    *   engines this service is going to own; if it results in an empty list,
    *   a single instance will be created.
    * * *perEventSeeds* (boolean, default: false): don't check that seeds
    *   match expectations on construction; this is needed in per-event policies
    *   where the seed is not known until the events are available
    * 
    */
   class GlobalEngineUserTestService {
         public:
      GlobalEngineUserTestService
        (fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      
         private:
      
      std::vector<std::string> instanceNames; ///< name of engine instances
      
      bool perEventSeeds; ///< whether to skip seed check in constructor
      
      /// our random generator engines
      std::vector<std::unique_ptr<TRandom>> engines;
      
      /// Checks all engines by CheckSeed() calls
      void CheckAllSeeds() const;
      
      /// Throws an exception if the seed is not expected
      static void CheckSeed(TRandom const& engine);
      
      /// Name used for the global engine instance
      static const std::string GlobalInstanceName;
      
      
      // callbacks
      void preModuleConstruction (art::ModuleDescription const&);
      void postModuleConstruction(art::ModuleDescription const&);
      void preModuleBeginRun     (art::ModuleDescription const&);
      void postModuleBeginRun    (art::ModuleDescription const&);
      void preProcessEvent       (art::Event const& evt);
      void preModule             (art::ModuleDescription const&);
      void postModule            (art::ModuleDescription const&);
      void postProcessEvent      (art::Event const&);
      void preModuleEndJob       (art::ModuleDescription const&);
      void postModuleEndJob      (art::ModuleDescription const&);
      
   }; // class GlobalEngineUserTestService
   
} // namespace testing


DECLARE_ART_SERVICE(testing::GlobalEngineUserTestService, LEGACY)


#endif // ARTEXT_SEEDSERVICE_TEST_GLOBALENGINEUSERSERVICE_H
