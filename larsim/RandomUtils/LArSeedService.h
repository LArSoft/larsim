/**
 * @file LArSeedService.h
 * @brief An art service to assist in the distribution of guaranteed unique seeds to all engines within an art job.
 * @author Rob Kutschke (kutschke@fnal.gov)
 * @date   2013/03/14 19:54:49
 * @see SeedService_service.cc
 * 
 */

/******************************************************************************
 *** TODO QA
 *** 
 *** - the original code had copy constructor made private to prevent it from
 ***   being called (now it is deleted); is this a prescription for a service?
 *** 
 ******************************************************************************/

#ifndef SeedService_SeedService_h
#define SeedService_SeedService_h

// C/C++ standard libraries
#include <string>
#include <map>
#include <iomanip>
#include <initializer_list>

// Some helper classes.
#include "larsim/RandomUtils/ArtState.h"
#include "larsim/RandomUtils/Providers/SeedMaster.h"

// CLHEP libraries
#ifdef SEEDSERVICE_USECLHEP
#	include "CLHEP/Random/RandomEngine.h"
#endif // SEEDSERVICE_USECLHEP

// ROOT libraries
#ifdef SEEDSERVICE_USEROOT
#	include "TRandom.h"
#endif // SEEDSERVICE_USEROOT

// From art and its tool chain.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EngineCreator.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

// Forward declarations
namespace art {
  class ActivityRegistry;
  class ModuleDescription;
  class Run;
  class SubRun;
}

namespace sim {
  
  /**
   * @brief An art service to assist in the distribution of guaranteed unique
   * seeds to all engines within an art job.
   * @see SeedMaster
   *
   * The LArSeedService acts as an interface between art framework and the
   * SeedMaster class.
   * 
   * The documentation is mantained in the SeedMaster class.
   * The configuration of LArSeedService is exactly the same as SeedMaster's,
   * and in art it's read from `services.user.LArSeedService`.
   * In a art module, code requests a seed by making one of the following two
   * calls:
   *     
   *     art::ServiceHandle<LArSeedService>->getSeed();
   *     art::ServiceHandle<LArSeedService>->getSeed("instanceName");
   *     
   * This must happen in the constructor of the module.
   * The reason is that we don't want engines to be initialized in the middle
   * of a job.
   * 
   * 
   * Interface with art's RandomNumberGenerator service
   * ---------------------------------------------------
   * 
   * It is possible to use the interface of this service to fully manager the
   * random engines. The internally preferred engine manager is still
   * RandomNumberGenerator, but other engines can be integrated to be managed
   * by LArSeedService as well.
   * 
   * The normal way to obtain random number engines is to ask
   * RandomNumberGenerator service during module's construction:
   *     
   *     // within module construction:
   *     createEngine(seed, "HepJamesRandom", "instanceA");
   *     
   * that will create an engine associated to the current module and seed it.
   * The engine will be of type CLHEP::HepRandomEngine.
   * The equivalent line with LArSeedService interface is:
   *     
   *     // within module construction:
   *     art::ServiceHandle<sim::LArSeedService>()->createEngine
   *       (*this, "HepJamesRandom", "instanceA");
   *     
   * This will in fact internally call the same function above.
   * The gain is automatic seeding of the engine by LArSeedService and the fact
   * that from now on LArSeedService knows about the engine. That is traded with
   * a more cumbersome syntax and the loss of control on the seed.
   * The syntax is the result of the need to talk to a service and the fact that
   * RandomNumberGenerator has a very restricted interface.
   * About the second point, a bit more control is offered by an extended
   * syntax:
   *     
   *     // within module construction:
   *     art::ServiceHandle<sim::LArSeedService>()->createEngine
   *       (*this, "HepJamesRandom", "instanceA", pset, { "Seed", "MySeed" });
   *     
   * Here LArSeedService will provide a seed only if it does not find in the
   * parameter set `pset` any of the specified configuration parameters (first
   * "Seed", then "MySeed" in this example). In other words, if "Seed" exists,
   * its value is used as seed, otherwise if "MySeed" exists, its value is used
   * instead, and otherwise LArSeedService is given control on that seed.
   * The exception is that if the specified seed is a "magic value", the
   * `InvalidSeed` (`0`), it is interpreted as a request to ignore the parameter
   * and use the service to get the seed. This is made as a quick way to remove
   * the seed override from an existing FHiCL file with one line.
   * Note that if LArSeedService does not get the control, also the policies that
   * reseed on event-by-event basis will not do that.
   * 
   * There are a few cases where RandomNumberGenerator is not controlling the
   * random engine. In that case, LArSeedService can still interact with the engine
   * (that has to have been created already):
   *     
   *     // within module construction:
   *     art::ServiceHandle<sim::LArSeedService>()->registerEngine(
   *       sim::LArSeedService::TRandomSeeder(pRandom),
   *       "instanceA", pset, "Seed"
   *       );
   *     
   * The first parameter is a function or functor needed to set the seed of the
   * engine (in the example, it's a optional functor LArSeedService provides to
   * deal with ROOT's TRandom objects, but a custom function can be used if
   * needed).
   * If that function is not provided, LArSeedService will not be able to set the
   * seed of the engine, neither at the beginning of the job nor on every new
   * event, as some policies require.
   */
  class LArSeedService {
      public:
    using seed_t = art::RandomNumberGenerator::seed_t;
    
    using SeedMaster_t = SeedMaster<seed_t>; ///< type of object providing seeds
    
    using EngineId = SeedMaster_t::EngineId; ///< type of random engine ID
    
    /// An invalid seed
    static constexpr seed_t InvalidSeed = SeedMaster_t::InvalidSeed;
    
    LArSeedService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    
    // Accept compiler written d'tor.  Not copyable or assignable.
    // This class is not copyable or assignable: these methods are not implemented.
    LArSeedService(LArSeedService const&) = delete;
    LArSeedService const& operator=(LArSeedService const&) = delete;
    LArSeedService(LArSeedService&&) = delete;
    LArSeedService const& operator=(LArSeedService&&) = delete;
    ~LArSeedService() = default;
    
    /// Returns the seed value for this module label (instance name).
    seed_t getSeed();
    seed_t getSeed(std::string instanceName);
    
    
    /// @{
    /// @name RandomNumberGenerator interface replacement
    /// 
    /// These functions allow to interface with the RandomNumberGenerator.
    /// In this way, LArSeedService will take care of most of the needs for
    /// seeding.
    /// 
    
    //@{
    /**
     * @brief Creates an engine with RandomNumberGenerator service
     * @param module module who will own the new engine
     * @param type the type of engine
     * @param instance the name of the engine instance
     * @return the seed used
     * @see registerEngine()
     *
     * This method creates a new engine by calling
     * RandomNumberGenerator::createEngine() with the seed from getSeed().
     * The meaning of the two parameters is the same as in that function,
     * including the somehow inconvenient order of the arguments
     * 
     * @note The module parameter is needed since the interface to create
     * an engine by RandomNumberGenerator service is private and open only
     * to friends.
     */
    seed_t createEngine
      (art::EngineCreator& module, std::string type, std::string instance = "");
    seed_t createEngine(art::EngineCreator& module);
    //@}
    
    //@{
    /**
     * @brief Creates an engine with RandomNumberGenerator service
     * @param module module who will own the new engine
     * @param type the type of engine
     * @param instance the name of the engine instance
     * @param pset parameter set to read parameters from
     * @param pname name or names of the seed parameters
     * @return the seed used
     * @see registerEngine()
     *
     * This method creates a new engine by calling
     * RandomNumberGenerator::createEngine() with the seed from getSeed().
     * The meaning of the two parameters is the same as in that function,
     * including the somehow inconvenient order of the arguments
     * 
     * The engine seed is set. First, the seed is retrieved from the specified
     * configuration, looking for the first of the parameters in pname that is
     * avaialble. If no parameter is found, the seed is obtained from
     * LArSeedService.
     * 
     * @note The module parameter is needed since the interface to create
     * an engine by RandomNumberGenerator service is private and open only
     * to friends.
     */
    seed_t createEngine(
      art::EngineCreator& module, std::string type, std::string instance,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return createEngine(module, type, instance, pset, { pname }); }
    seed_t createEngine(
      art::EngineCreator& module, std::string type, std::string instance,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
      );
    seed_t createEngine(
      art::EngineCreator& module, std::string type,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return createEngine(module, type, "", pset, pname); }
    seed_t createEngine(
      art::EngineCreator& module, std::string type,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
      )
      { return createEngine(module, type, "", pset, pname); }
    seed_t createEngine(
      art::EngineCreator& module,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return createEngine(module, pset, { pname }); }
    seed_t createEngine(
      art::EngineCreator& module,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
      );
    //@}
    
    /**
     * @brief Registers an existing engine with LArSeedService
     * @param seeder function used to set the seed of the existing engine
     * @param instance name of the engine
     * @see createEngine()
     *
     * This function works similarly to createEngine(), but it uses an existing
     * engine instead of creating a new one by RandomNumberGenerator service.
     * The seeder function must be provided for the service to be of any use:
     * registerEngine() will set the seed immediately, and the seeder function
     * will be used to set the seed for policies that do that on each event.
     * The instance name must also be unique, since for LArSeedService purposes
     * the registered engine is no different from any other, created by
     * RandomNumberGenerator or not.
     * 
     * Three standard functions are provided as seeders, for use with
     * RandomNumberGenerator engines (RandomNumberGeneratorSeeder()),
     * with a CLHEP::HepRandomEngine (CLHEPengineSeeder class) and with ROOT's
     * TRandom (TRandomSeeder class). Note that CLHEP and ROOT classes are not
     * compiled in LArSeedService by default, and the recommendation is to take
     * their implementation as an example and create your own after them).
     * Any seeder function with the prototype of LArSeedService::Seeder_t:
     *     
     *     void Seeder(EngineId const&, seed_t);
     *     
     * or a functor with
     *     
     *     void operator() (EngineId const&, seed_t);
     *     
     * can be used as seeder.
     */
    seed_t registerEngine
      (SeedMaster_t::Seeder_t seeder, std::string instance = "" );
    
    //@{
    /**
     * @brief Registers an existing engine with LArSeedService
     * @param seeder function used to set the seed of the existing engine
     * @param instance name of the engine
     * @param pset parameter set to read parameters from
     * @param pname name or names of the seed parameters
     * @see createEngine()
     *
     * These functions work similarly to registerEngine(), but the preferred
     * way to obtain the seed is from configuration.
     * First, the seed is retrieved from the specified
     * configuration, looking for the first of the parameters in pname that is
     * avaialble. If no parameter is found, the seed is obtained from
     * LArSeedService.
     */
    seed_t registerEngine(
      SeedMaster_t::Seeder_t seeder, std::string instance,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return registerEngine(seeder, instance, pset, { pname }); }
    seed_t registerEngine(
      SeedMaster_t::Seeder_t seeder, std::string instance,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
      );
    seed_t registerEngine(
      SeedMaster_t::Seeder_t seeder,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return registerEngine(seeder, "", pset, pname); }
    seed_t registerEngine(
      SeedMaster_t::Seeder_t seeder,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
      )
      { return registerEngine(seeder, "", pset, pname); }
    //@}
    
    /// @}
    
    /// Prints known (EngineId,seed) pairs.
    template<class Stream> void print(Stream&) const;
    
    /// Prints to the framework Info logger
    void print() const;
    
    // Call backs that will be called by art.
    void preModuleConstruction (art::ModuleDescription const& md);
    void postModuleConstruction(art::ModuleDescription const&);
    void preModuleBeginRun     (art::ModuleDescription const& md);
    void postModuleBeginRun    (art::ModuleDescription const&);
    void preProcessEvent       (art::Event const& evt);
    void preModule             (art::ModuleDescription const& md);
    void postModule            (art::ModuleDescription const&);
    void postProcessEvent      (art::Event const&);
    void postEndJob            ();
    
    
    /// Seeder_t function setting the seed of an engine in RandomNumberGenerator
    static void RandomNumberGeneratorSeeder(EngineId const& id, seed_t seed);
    
#ifdef SEEDSERVICE_USEROOT
    /// Seeder_t functor setting the seed of a ROOT TRandom engine (untested!)
    class TRandomSeeder {
        public:
      TRandomSeeder(TRandom* engine): pRandom(engine) {}
      void operator() (EngineId const&, seed_t seed)
        { if (pRandom) pRandom->SetSeed(seed); }
        protected:
      TRandom* pRandom = nullptr;
    }; // class TRandomSeeder
#endif // SEEDSERVICE_USEROOT
    
#ifdef SEEDSERVICE_USECLHEP
    /// Seeder_t functor setting the seed of a CLHEP::HepRandomEngine engine (untested!)
    class CLHEPengineSeeder {
        public:
      CLHEPengineSeeder(CLHEP::HepRandomEngine* engine): pRandom(engine) {}
      void operator() (EngineId const&, seed_t seed)
        { if (pRandom) pRandom->setSeed(seed, 0); }
        protected:
      CLHEP::HepRandomEngine* pRandom = nullptr;
    }; // class CLHEPengineSeeder
#endif // SEEDSERVICE_USECLHEP
    
  private:
    
    
    /// Class managing the seeds
    SeedMaster_t seeds;
    
    /// Helper to track state of art. It is legal to call getSeed only in certain states.
    SeedServiceHelper::ArtState state;
    
    /// Control the level of information messages.
    int verbosity = 0;
    bool bPrintEndOfJobSummary = false; ///< print a summary at the end of job
    
    
    /// Returns whether the specified engine is already registered
    bool hasEngine(EngineId const& id) const { return seeds.hasEngine(id); }
    
    // Main logic for computing and validating a seed.
    seed_t getSeed(EngineId const&);
    
    // Main logic for computing and validating a seed.
    seed_t getEventSeed(EngineId const&);
    
    /**
     * @brief Reseeds the specified engine instance in the current module
     * @param instance the name of the engine instance
     * @return the seed set, or InvalidSeed if no reseeding happened
     */
    seed_t reseedInstance(EngineId const& id);
    
    //@{
    /// Reseeds all the engines in the current module
    void reseedModule(std::string currentModule);
    void reseedModule();
    //@}
    
    
    /// Registers the engine ID into SeedMaster
    seed_t prepareEngine(EngineId const& id, SeedMaster_t::Seeder_t seeder);
    
    // Helper functions for all policies
    void ensureValidState() const;
    
    //@{
    /// Returns a fully qualified EngineId
    EngineId qualify_engine_label
      (std::string moduleLabel, std::string instanceName) const;
    EngineId qualify_engine_label(std::string instanceName = "") const;
    //@}
    
    //@{
    /// Reads the seed from the first of the specified parameters available
    /// @return whether any parameter was found
    static bool readSeedParameter
      (seed_t& seed, fhicl::ParameterSet const& pset, std::string pname)
      { return readSeedParameter(seed, pset, { pname }); }
    static bool readSeedParameter(
      seed_t& seed, fhicl::ParameterSet const& pset,
      std::initializer_list<std::string> pnames
      );
    //@}
    
    /// Helper to retrieve a seed including configuration
    /// @return the seed, and whether it is fixed (that is, from configuration)
    std::pair<seed_t, bool> getSeed(
      EngineId const&id,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
      );
    
    /// Forces LArSeedService not to change the seed of the specified engine
    void freezeSeed(EngineId const& id);
    
    /// Registers an engine and its seeder
    void registerEngineAndSeeder
      (EngineId const& id, SeedMaster_t::Seeder_t seeder);
    
 }; // class LArSeedService
  
  
  
  //================================================================
  template<class Stream>
  inline void LArSeedService::print(Stream& log) const {
    seeds.print(log);
  } // print(Stream)
  
} // namespace sim

DECLARE_ART_SERVICE(sim::LArSeedService, LEGACY)

#endif /* SeedService_SeedService_h */
