/**
 * @file LArSeedService.h
 * @brief An art service to assist in the distribution of guaranteed unique seeds to all engines within an art job.
 * @author Gianluca Petrillo (petrillo@fnal.gov), Rob Kutschke (kutschke@fnal.gov)
 * @date   2013/03/14 19:54:49
 * @see SeedService_service.cc
 * 
 */


#ifndef LARSIM_RANDOMUTILS_LARSEEDSERVICE_H
#define LARSIM_RANDOMUTILS_LARSEEDSERVICE_H 1

#ifndef LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP
/// Define to non-zero to include special CLHEP random engine support
#  define LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP 0
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP

// ROOT libraries
#ifndef LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT
/// Define to non-zero to include special ROOT random generator support
#  define LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT 0
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT


// C/C++ standard libraries
#include <string>
#include <utility> // std::forward()
#include <initializer_list>

// Some helper classes.
#include "larsim/RandomUtils/ArtState.h"
#include "larsim/RandomUtils/Providers/SeedMaster.h"

// CLHEP libraries
#if (LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP)
#  include "CLHEP/Random/RandomEngine.h"
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP

// ROOT libraries
#if (LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT)
#  include "TRandom.h"
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT

// From art and its tool chain.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

// Forward declarations
namespace art {
  class ActivityRegistry;
  class ModuleDescription;
  class Run;
  class SubRun;
  class EngineCreator;
}

namespace sim {
  
  /**
   * @brief An art service to assist in the distribution of guaranteed unique
   * seeds to all engines within an art job.
   * @see SeedMaster
   *
   * LArSeedService centrally manages seeds for random generator engines.
   * 
   * The LArSeedService acts as an interface between art framework and the
   * SeedMaster class.
   * 
   * The documentation is mantained in the SeedMaster class.
   * The configuration of LArSeedService is exactly the same as SeedMaster's,
   * and in art it's read from `services.LArSeedService`.
   * The following documentation describes features of LArSeedService that are
   * built on top of SeedMaster to have a more convenient interaction within
   * the art framework.
   * 
   * Before asking LArSeedService for its seed, an engine must be in some way
   * registered. Once the engine is registered, its original seed can be queried
   * again by calling `getSeed()` methods.
   * 
   * 
   * Glossary
   * ---------
   * 
   * Here "engine" means a class that is able to generate random numbers
   * according to a flat distribution.
   * Both art and LArSeedService are module-based, that means that the engines
   * are in the context of a specific module instance, and different module
   * instances have independent engines.
   * That is the reason why you don't need to specify anything about the module
   * when creating or obtaining a random engine, and it is also the reason why
   * engines outside module context are not supported by the framework.
   * 
   * Each module can need more than one engine. A module can have any number of
   * engines, and each of them is identified by an "instance name" that is
   * unique within the module.
   * Nonetheless, most modules need just one engine. In that case, a default
   * instance name can be used (that is an empty string).
   * 
   * A "seeder" is a callable object (e.g. a function) that sets the seed of
   * a certain engine. The seeder is expected to find out by its own which
   * engine it has to seed, and fot that it is provided an engine ID.
   * 
   * 
   * Registration of a random generator engine
   * ------------------------------------------
   * 
   * Registration must happen in art module constructor, in one of the following
   * ways:
   *  * by asking this service to create an engine via RandomNumberGenerator
   *    (see `createEngine()` methods) [discouraged]
   *  * by registering an existing engine and its seeding function
   *    (see `registerEngine()` methods)
   *  * by just declaring that an engine exists
   *    (see `declareEngine()` and `getSeed()` methods)
   * The first method and, when a seeder or an engine is provided, also the
   * second method, set the seed of the engine they register (see below).
   * In the third case, it is generally the caller's responsibility to seed the
   * engine. The registration of an engine which has been only declared can be
   * "completed" by calling `defineEngine()` to provide the actual seeder for
   * that engine. The pair of calls `declareEngine()`/`defineEngine()` (or
   * `getSeed()`/`defineEngine()`) is equivalent to a single call to
   * `registerEngine()`, with the added flexibility of having the seed for the
   * engine already available before the registration is completed.
   * 
   * The use of createEngine() class of function is discouraged and this
   * function might be removed in the future, because of its non-clean use of
   * module interfaces.
   * The recommended approach is more prolix; here an example for an engine with
   * a non-default instance name:
   *     
   *     std::string const instanceName = "instanceName";
   *     auto& Seeds = *(art::ServiceHandle<sim::LArSeedService>());
   *     
   *     // declare an engine; LArSeedService associates an (unknown) engine, in
   *     // the current module and an instance name, with a seed (returned)
   *     auto const seed = Seeds.declareEngine(instanceName);
   *     
   *     // now create the engine (for example, use art); seed will be set
   *     createEngine(seed, "HepJamesRandom", instanceName);
   *     
   *     // finally, complete the registration; seed will be set again
   *     art::ServiceHandle<art::RandomNumberGenerator> RNG;
   *     Seeds.defineEngine(RNG->getEngine(instanceName));
   *     
   * This is equivalent to the discouraged call
   *     
   *     auto& Seeds = *(art::ServiceHandle<sim::LArSeedService>());
   *     Seeds.createEngine(*this, "HepJamesRandom", "instanceName");
   *     
   * Please read carefully the documentation of the method of your choice, since
   * they have different requirements and apply to different usage patterns.
   * 
   * The registration must happen in the constructor of the module.
   * That is because we don't want engines to be initialized in the middle of a
   * job.
   * 
   * 
   * Setting the seed of an engine
   * ------------------------------
   * 
   * LArSeedService is able to set the seed of an engine when the engine is
   * registered via either:
   *  * `createEngine()` (creation of a new CLHEP engine)
   *  * `registerEngine()` (registration of an engine or a seeder function),
   *    if the registered seeder function is valid (non-null) or if a CLHEP
   *    engine is being registered (in which case the seeder is automatically
   *    created valid)
   *  * `defineEngine()` (registration of a seeder for an engine that was
   *    already declared), again if the seed is valid
   * LArSeedService is *not* able to automatically set the seed of an engine if
   * it was registered via either:
   *  * `declareEngine()` (declaration of the existence of an engine),
   *    that does not even require the engine to exist
   *  * `getSeed()` (query of a seed), when it (implicitly) declares an engine
   *    which had not been declared yet
   * 
   * If LArSeedService is able to set the seed, it will do so only once, as soon as
   * it can.
   * This means that if the policy allows the seed to be known immediately,
   * the seed will be set on registration. In the case of a per-event policy
   * that requires the presence of an event, the seed can be known only 
   * when the event is available, and LArSeedService will set the seed before the
   * module the engine is associated with starts its main processing method
   * (produce(), filter() or analyze()).
   * 
   * 
   * Changing the seeder
   * --------------------
   * 
   * Currently, changing the seeder of an engine after the engine has been
   * fully registered is not supported. As a consequence, changing the engine
   * is also not supported.
   * 
   * Since only the seeder function is registered in LArSeedService, a seeder
   * function that is flexible enough to change the engine it seeds may work
   * around this limitation.
   * 
   * 
   * Querying the seed of an engine
   * -------------------------------
   * 
   * If necessary, the seed that LArSeedService has assigned to an engine can be
   * requested by one of the following two calls:
   *     
   *     art::ServiceHandle<LArSeedService>->getSeed();
   *     art::ServiceHandle<LArSeedService>->getSeed("instanceName");
   *     
   * depending on whether the engine has a non-empty instance name,
   * Note that this call implicitly "declares" the engine it refers to.
   * 
   * For most policies, the seed is set according to the configuration, once
   * for all. In those cases, `getSeed()` will always return the same value.
   * If the policy prescribes different seeds at different times, the method
   * returns the seed that is assigned to the engine at the time of the call.
   * 
   * Also note that the seed assigned by LArSeedService might not match the current
   * seed of the engine, if:
   *  * LArSeedService is not in charge of setting the seed of the engine,
   *    and the engine seed has not been set yet
   *  * the seed was reset directly after LArSeedService set the engine seed
   * Both circumstances should be avoided.
   * 
   * 
   * Creating the engines independently of LArSeedService
   * --------------------------------------------------
   * 
   * A number of things must happen for an engine to correctly work with
   * LArSeedService:
   *  * the engine instance needs to exist or be created 
   *  * the engine must be "registered" into LArSeedService
   *  * the seed must be obtained from LArSeedService
   *  * the seed must be provided to the engine
   * 
   * A recipe for creating a ROOT engine and let LArSeedService take care of its
   * seeds is:
   *     
   *     // create the engine; ROOT will set some (temporary) seed
   *     fRandom = std::make_unique<TRandom3>();
   *     
   *     // declare the engine; LArSeedService associates its seeder, in
   *     // the current module and an instance name, with a seed (returned);
   *     // the seed is also set in the engine
   *     auto& Seeds = *(art::ServiceHandle<sim::LArSeedService>());
   *     Seeds.registerEngine(TRandomSeeder(fRandom), "instanceName");
   *     
   * Here `fRandom` is supposed to be a member function of the module.
   * The `TRandomSeeder` object is an object that knows how to set the seed
   * of the TRandom-derived ROOT generator passed as constructor argument.
   * For an example of implementation, see the source code of LArSeedService.
   * 
   * 
   * Overriding the seed from LArSeedService at run time
   * -------------------------------------------------
   * 
   * LArSeedService (and SeedMaster, which the former relies upon) will decide
   * which seed to give to each registered engine and, when possible
   * (see above), will set that seed too.
   * 
   * All registration functions offer an extended signature to tell LArSeedService
   * that if there is an explicitly configured seed, that should take presedence
   * over the one automatically assigned by SeedMaster policy.
   * This extended signature includes:
   *  * a FHiCL parameter set
   *  * a configuration parameter name, or a list of them
   * LArSeedService will look in the specified parameter set and if it finds a
   * value corresponding to any of the specified parameter names, will set
   * the seed of the engine to that value, and it will mark the engine as
   * "frozen" (meaning that LArSeedService will not ever set a seed again on that
   * engine). [see also the exception to this rule below]
   * 
   * The typical use of this function is to read a parameter ("Seed" is the
   * suggested name) from the configuration of the module. Note that this is in
   * contrast with the location where LArSeedService seeds are normally configured,
   * that is in the configuration of LArSeedService service itself.
   * For example:
   *     
   *     // within module construction:
   *     art::ServiceHandle<sim::LArSeedService>()->registerEngine
   *       (engine, "instanceName", pset, { "Seed", "MySeed" });
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
   * Note that if LArSeedService does not get the control on the seed,
   * policies that reseed on event-by-event basis will not act on the engine.
   * 
   * 
   * Engines outside of the module context
   * ======================================
   * 
   * It is possible to have engines that are not associated with any module.
   * If no module is current, an engine will be registered in a "global"
   * context. This happens during service construction, and once the service
   * construction phase is completed, no more global engines can be registered.
   * Would one ever need to access the seed of such engine, a specific
   * interface needs to be used: `getGlobalSeed()` to get the configured seed
   * at the beginning of the job, or `getGlobalCurrentSeed()` to get the seed
   * specific for the current event, if any. These are equivalent to the module
   * context methods `getSeed()` and `getCurrentSeed()`.
   * 
   * The art service RandomNumberGenerator does not support the creation of
   * module-independent engines. The ownership of each global engine is by the
   * using service, as well as the management of engine's lifetime.
   * After an engine has been instantiated, it can be registered with
   * `registerEngine()`. Likewise, a global engine can be declared first with
   * `declareEngine()`, then instantiated, and then optionally defined with
   * `defineEngine()`. This is completely analogous to the module-context
   * engines. The only difference is that no `createEngine()` interface is
   * available for global engines.
   * Whether these methods create a global or module engine depends only
   * on the context they are called in.
   * 
   * LArSeedService does not manage the engine life in any way. If a service owns
   * an engine, it also needs a way to access it, as nothing equivalent to
   * RandomNumberGenerator's `getEngine()` is provided.
   * LArSeedService does manage the seeding of the registered engines, even the
   * global ones. If the seed policy involves a event-dependent seed, all global
   * engines are seeded together at the beginning of the event, before any
   * module is executed. This mechanism relies on the fact that LArSeedService gets
   * its preProcessEvent callback executed before the ones of the services that
   * own any engine. This is guaranteed if the service constructors invoke
   * LArSeedService before they register their callbacks.
   * 
   * For an example of usage, see GlobalEngineUserTestService service in the
   * test suite. Here is an excerpt using a ROOT TRandom3 engine, that can
   * be constructed without seed. In the service constructor, the code:
   *     
   *     // create a new engine, that we own (engine be a class data member)
   *     engine = std::make_unique<TRandom3>();
   *     auto seed = Seeds.registerEngine
   *       (sim::LArSeedService::TRandomSeeder(engine.get()), "MyService");
   *     
   *     // MyService callback registrations in the ActivityRegistry
   *     // should follow the first call to LArSeedService!
   *     
   * where sim::LArSeedService::TRandomSeeder is a seeder class for TRandom
   * engines (optionally provided in LArSeedService.h, and very easy to write).
   * This excerpt creates and owns a TRandom3 engine, and then it registers it
   * to LArSeedService, which immediately sets its seed and will take care of
   * reseeding the engine at each event, should the policy require it.
   * The service will access the engine as a data member, and should it need the
   * seed it will use:
   *     
   *     art::ServiceHandle<sim::LArSeedService>()->getGlobalSeed("MyService");
   *     
   * Note that is a good idea to give the engine an instance name after the
   * service name itself, since the instance name is global and it's the only
   * thing distinguishing global engines, and name conflicts between different
   * services may easily arise.
   * 
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
    
    /// Returns whether the specified seed is valid
    static constexpr bool isSeedValid(seed_t seed)
      { return seed != InvalidSeed; }
    
    
    /**
     * @brief Returns a seed for the engine with specified instance name
     * @param instanceName name of the engine instance
     * @return a seed for the engine with specified instance name
     * @see getGlobalSeed()
     * 
     * The seed for an engine in the context of the current module is returned.
     * If you need the seed for an engine outside that context, use
     * `getGlobalSeed()` instead.
     * 
     * The engine needs to have been registered before, in any of the supported
     * ways.
     * If it has not, this call will declare it with declareEngine()
     * and no further registration will be allowed.
     * 
     * While this method can be called at any time, the registration of an
     * engine can happen only at construction time, and it will make the call to
     * this method fail if it is called at any other time.
     */
    seed_t getSeed(std::string instanceName);
    
    /**
     * @brief Returns a seed for the engine with default instance name
     * @return a seed for the engine with default instance name
     * 
     * This method is equivalent to getSeed(std::string) with an empty instance
     * name.
     */
    seed_t getSeed();
    
    
    /**
     * @brief Returns a seed for the global engine with specified instance name
     * @param instanceName name of the engine instance
     * @return a seed for the global engine with specified instance name
     * @see getSeed()
     * 
     * A "global" engine is not bound to a specific execution context.
     * The only context LArSeedService is aware of is the module, so this
     * translates into engines that are not bound to any module.
     * To instruct LArSeedService to ignore the current context (that may be
     * a running module, or no running module at all), `getGlobalSeed()` is
     * used instead of `getSeed()`, that will consider the context and in fact
     * consider the absence of context an error.
     * 
     * The engine needs to have been registered before, in any of the supported
     * ways.
     * If it has not, this call will declare it with declareEngine()
     * and no further registration will be allowed.
     * 
     * While this method can be called at any time, the registration of an
     * engine can happen only at construction time, and it will make the call to
     * this method fail if it is called at any other time.
     */
    seed_t getGlobalSeed(std::string instanceName);
    
    
    /// Returns the last computed seed for specified engine of current module
    seed_t getCurrentSeed(std::string instanceName) const
      { return seeds.getCurrentSeed(qualify_engine_label(instanceName)); }
    
    /// Returns the last computed seed for the default engine of current module
    seed_t getCurrentSeed() const
      { return seeds.getCurrentSeed(qualify_engine_label()); }
    
    /// Returns the last computed seed for the specified global engine
    seed_t getGlobalCurrentSeed(std::string instanceName) const
      { return seeds.getCurrentSeed(qualify_global_engine(instanceName)); }
    
    
    /// @{
    /**
     * @name Create and register an engine
     * 
     * The life time of the engine is managed by art::RandomNumberGenerator,
     * while the seeding is managed by this service.
     * 
     * This is a replacement of art::EngineCreator-derived classes.
     * The use of createEngine() class of function is discouraged and this
     * function might be removed in the future, because of its non-clean use of
     * module interfaces.
     * The recommended approach is more verbose;
     * here an example for an engine with a non-default instance name:
     *     
     *     std::string const instanceName = "instanceName";
     *     auto& Seeds = *(art::ServiceHandle<sim::LArSeedService>());
     *     
     *     // declare an engine; LArSeedService associates an (unknown) engine, in
     *     // the current module and an instance name, with a seed (returned)
     *     auto const seed = Seeds.declareEngine(instanceName);
     *     
     *     // now create the engine (for example, use art); seed will be set
     *     createEngine(seed, "HepJamesRandom", instanceName);
     *     
     *     // finally, complete the registration; seed will be set again
     *     art::ServiceHandle<art::RandomNumberGenerator> RNG;
     *     Seeds.defineEngine(RNG->getEngine(instanceName));
     *     
     * This is equivalent to the discouraged call
     *     
     *     auto& Seeds = *(art::ServiceHandle<sim::LArSeedService>());
     *     Seeds.createEngine(*this, "HepJamesRandom", "instanceName");
     *     
     */
    
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
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      );
    seed_t createEngine(
      art::EngineCreator& module, std::string type,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return createEngine(module, type, "", pset, pname); }
    seed_t createEngine(
      art::EngineCreator& module, std::string type,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      )
      { return createEngine(module, type, "", pset, pnames); }
    seed_t createEngine(
      art::EngineCreator& module,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return createEngine(module, pset, { pname }); }
    seed_t createEngine(
      art::EngineCreator& module,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      );
    //@}
    
    /// @}
    
    /// @{
    /// @name Register an existing engine
    /// 
    /// The life time of the engine is under user's control, while the seeding
    /// is managed by this service.
    /// 
    
    /**
     * @brief Registers an existing engine with LArSeedService
     * @param seeder function used to set the seed of the existing engine
     * @param instance name of the engine
     * @return the seed assigned to the engine (may be invalid)
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
     * @return the seed assigned to the engine (may be invalid)
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
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      );
    seed_t registerEngine(
      SeedMaster_t::Seeder_t seeder,
      fhicl::ParameterSet const& pset, std::string pname
      )
      { return registerEngine(seeder, "", pset, pname); }
    seed_t registerEngine(
      SeedMaster_t::Seeder_t seeder,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      )
      { return registerEngine(seeder, "", pset, pnames); }
    //@}
    
#if (LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP)
    /**
     * @brief Registers an existing CLHEP engine with LArSeedService
     * @param engine a reference to the CLHEP random generator engine
     * @param instance name of the engine
     * @param pset parameter set to read parameters from
     * @param pname name or names of the seed parameters
     * @return the seed assigned to the engine (may be invalid)
     *
     * The specified engine is not managed.
     * It may be owned by RandomNumberGenerator service.
     * 
     * The engine is expected to be valid as long as this service performs
     * reseeding.
     */
    seed_t registerEngine
      (CLHEP::HepRandomEngine& engine, std::string instance = "")
      { return registerEngine(CLHEPengineSeeder(engine), instance); }
    seed_t registerEngine(
      CLHEP::HepRandomEngine& engine, std::string instance,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      )
      {
        return registerEngine
          (CLHEPengineSeeder(engine), instance, pset, pnames);
      }
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP
    /// @}
    
    
    /// @{
    /// @name Declare the presence of an engine
    
    /**
     * @brief Declares the presence of an engine with a given instance name
     * @param instance name of the instance of the engine (empty by default)
     * @return the seed assigned to the engine (may be invalid)
     * 
     * The existence of an engine with the specified instance name is recorded,
     * and a seed is assigned to it. The engine will be identified by the
     * instance name and by context information (the current module).
     * 
     * Differently from createEngine() and registerEngine(), the actual
     * existence of a engine is not required. It is up to the user to manage
     * the engine, if any at all, including the seeding.
     */
    seed_t declareEngine(std::string instance = "");
    
    /**
     * @brief Declares the presence of an engine with a given instance name
     * @param instance name of the instance of the engine
     * @param pset parameter set where to find a possible fixed seed request
     * @param pname the name of the parameter for the fixed seed request
     * @return the seed assigned to the engine (may be invalid)
     * 
     * The existence of an engine with the specified instance name is recorded,
     * and a seed is assigned to it. The engine will be identified by the
     * instance name and by context information (the current module).
     * 
     * The preferred way to obtain the seed is from configuration.
     * First, the seed is retrieved from the specified configuration,
     * looking for the first of the parameters in pname that is available.
     * If no parameter is found, the seed is obtained from LArSeedService.
     * 
     * Differently from createEngine() and registerEngine(), the actual
     * existence of a engine is not required. It is up to the user to manage
     * the engine, if any at all, including the seeding.
     */
    seed_t declareEngine
      (std::string instance, fhicl::ParameterSet const& pset, std::string pname)
      { return declareEngine(instance, pset, { pname }); }
    
    /**
     * @brief Declares the presence of an engine with a given instance name
     * @param instance name of the instance of the engine
     * @param pset parameter set where to find a possible fixed seed request
     * @param pnames name of the parameters for the fixed seed request
     * @return the seed assigned to the engine (may be invalid)
     * @see declareEngine(std::string, fhicl::ParameterSet const&, std::string)
     * 
     * This method provides the same function as 
     * declareEngine(std::string, fhicl::ParameterSet const&, std::string),
     * but it can pick the seed from the first parameter among the ones in pset
     * whose name is in pnames.
     */
    seed_t declareEngine(
      std::string instance,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      );
    
    /**
     * @brief Declares the presence of an engine with a default instance name
     * @param pset parameter set where to find a possible fixed seed request
     * @param pname the name of the parameter for the fixed seed request
     * @return the seed assigned to the engine (may be invalid)
     * @see declareEngine(fhicl::ParameterSet const&, std::string)
     * 
     * This method provides the same function as 
     * declareEngine(std::string, fhicl::ParameterSet const&, std::string),
     * but it gives the engine a empty instance name.
     */
    seed_t declareEngine(fhicl::ParameterSet const& pset, std::string pname)
      { return declareEngine("", pset, pname); }
    
    /**
     * @brief Declares the presence of an engine with a default instance name
     * @param pset parameter set where to find a possible fixed seed request
     * @param pnames name of the parameters for the fixed seed request
     * @return the seed assigned to the engine (may be invalid)
     * @see declareEngine(std::string, fhicl::ParameterSet const&, std::string)
     * 
     * This method provides the same function as 
     * declareEngine(std::string, fhicl::ParameterSet const&, std::string),
     * but it can pick the seed from the first parameter among the ones in pset
     * whose name is in pnames.
     * Also, it gives the engine a empty instance name.
     */
    seed_t declareEngine(
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      )
      { return declareEngine("", pset, pnames); }
    /// @}
    
    /// @{
    /**
     * @brief Defines a seeder for a previously declared engine
     * @param seeder seeder associated to the engine
     * @param instance name of engine instance (default: empty)
     * @return the seed assigned to the engine (may be invalid)
     * @see declareEngine()
     * 
     * The seeder is the same object as in registerEngine().
     * This function can be used to finalise the declaration of an engine.
     * If the engine was just declared with declareEngine() (as opposed to
     * registered with registerEngine() or created with createEngine()),
     * "defining" the engine will hook it to LArSeedService, that will take care
     * of setting seeds automatically when needed.
     * This step is not mandatory, but no automatic seeding will happen if it is
     * omitted.
     */
    seed_t defineEngine
      (SeedMaster_t::Seeder_t seeder, std::string instance = {});
    
    #if (LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP)
    /**
     * @brief Defines a seeder for a previously declared engine
     * @param instance name of engine instance
     * @param engine CLHEP engine to be associated to the instance
     * @return the seed assigned to the engine (may be invalid)
     * @see declareEngine()
     * 
     * This method operates on the default engine instance and performs the
     * same operations as defineEngine(std::string, Seeder_t).
     * A seeder is internally created for the CLHEP random engine.
     */
    seed_t defineEngine
      (CLHEP::HepRandomEngine& engine, std::string instance = {})
      { return defineEngine(CLHEPengineSeeder(engine), instance); }
    #endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP
    /// @}
    
    
    /// Prints known (EngineId,seed) pairs.
    template<class Stream>
    void print(Stream&& out) const
      { seeds.print(std::forward<Stream>(out)); }
    
    /// Prints to the framework Info logger
    void print() const { print(mf::LogInfo("LArSeedService")); }
    
    /// Seeder_t function setting the seed of an engine in RandomNumberGenerator
    static void RandomNumberGeneratorSeeder(EngineId const& id, seed_t seed);
    
#if (LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT)
    /// Seeder_t functor setting the seed of a ROOT TRandom engine (untested!)
    class TRandomSeeder {
        public:
      TRandomSeeder(TRandom* engine): pRandom(engine) {}
      void operator() (EngineId const&, seed_t seed)
        { if (pRandom) pRandom->SetSeed(seed); }
        protected:
      TRandom* pRandom = nullptr;
    }; // class TRandomSeeder
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USEROOT
    
#if (LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP)
    /// Seeder_t functor setting the seed of a CLHEP::HepRandomEngine engine (untested!)
    class CLHEPengineSeeder {
        public:
      CLHEPengineSeeder(CLHEP::HepRandomEngine& e): engine(e) {}
      CLHEPengineSeeder(CLHEP::HepRandomEngine* e): engine(*e) {}
      void operator() (EngineId const&, seed_t seed)
        {
          engine.setSeed(seed, 0);
          LOG_DEBUG("CLHEPengineSeeder")
            << "CLHEP engine: '" << engine.name() << "'[" << ((void*) &engine)
            << "].setSeed(" << seed << ", 0)";
        }
        protected:
      CLHEP::HepRandomEngine& engine;
    }; // class CLHEPengineSeeder
#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP
    
  private:
    
    
    /// Class managing the seeds
    SeedMaster_t seeds;
    
    /**
     * Helper to track state of art.
     * 
     * The state is updated by LArSeedService itself, and therefore knows only
     * about what it is notified about, when it is notified about.
     * For example, service construction phase may start before the service
     * was even constructed, but the state will be updated only on LArSeedService
     * construction.
     */
    SeedServiceHelper::ArtState state;
    
    /// Control the level of information messages.
    int verbosity = 0;
    bool bPrintEndOfJobSummary = false; ///< print a summary at the end of job
    
    /// Register an engine and seeds it with the seed from the master
    seed_t registerEngineID(
      EngineId const& id,
      SeedMaster_t::Seeder_t seeder = SeedMaster_t::Seeder_t()
      );
    
    /// Set the seeder of an existing engine
    seed_t defineEngineID(EngineId const& id, SeedMaster_t::Seeder_t seeder);


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
    
    /// Reseed all the global engines
    void reseedGlobal();
    
    /// Registers the engine ID into SeedMaster
    seed_t prepareEngine(EngineId const& id, SeedMaster_t::Seeder_t seeder);
    
    // Helper functions for all policies
    void ensureValidState(bool bGlobal = false) const;
    
    //@{
    /// Returns a fully qualified EngineId
    EngineId qualify_engine_label
      (std::string moduleLabel, std::string instanceName) const;
    EngineId qualify_engine_label(std::string instanceName = "") const;
    EngineId qualify_global_engine(std::string instanceName = "") const
      { return EngineId(instanceName, EngineId::global); }
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
    
    /// Query a seed from the seed master
    seed_t querySeed(EngineId const& id);
    
    /// Helper to retrieve a seed including configuration
    /// @return the seed, and whether it is fixed (that is, from configuration)
    std::pair<seed_t, bool> findSeed(
      EngineId const& id,
      fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
      );
    
    /// Forces LArSeedService not to change the seed of the specified engine
    void freezeSeed(EngineId const& id, seed_t frozen_seed);
    
    /// Registers an engine and its seeder
    void registerEngineAndSeeder
      (EngineId const& id, SeedMaster_t::Seeder_t seeder);
    
    /// Calls the seeder with the specified seed and engine ID
    seed_t seedEngine(EngineId const& id) { return seeds.reseed(id); }
    
    // Call backs that will be called by art.
    void preModuleConstruction (art::ModuleDescription const& md);
    void postModuleConstruction(art::ModuleDescription const&);
    void preModuleBeginRun     (art::ModuleDescription const& md);
    void postModuleBeginRun    (art::ModuleDescription const&);
    void preProcessEvent       (art::Event const& evt);
    void preModule             (art::ModuleDescription const& md);
    void postModule            (art::ModuleDescription const&);
    void postProcessEvent      (art::Event const&);
    void preModuleEndJob       (art::ModuleDescription const& md);
    void postModuleEndJob      (art::ModuleDescription const&);
    void postEndJob            ();
    
    
 }; // class LArSeedService
  
  
  
} // namespace sim

DECLARE_ART_SERVICE(sim::LArSeedService, LEGACY)

#endif // LARSIM_RANDOMUTILS_LARSEEDSERVICE_H
