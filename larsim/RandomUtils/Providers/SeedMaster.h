/**
 * @file   SeedMaster.h
 * @brief  A class to assist in the distribution of guaranteed unique seeds
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   20141111
 * @see    LArSeedService.h SeedMaster.cc
 */

#ifndef SeedService_SeedMaster_h
#define SeedService_SeedMaster_h

// C/C++ standard libraries
#include <vector>
#include <string>
#include <map>
#include <memory> // std::unique_ptr<>

// From art and its tool chain
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/Exception.h"

// Some helper classes
#include "larsim/RandomUtils/Providers/MapKeyIterator.h"
#include "larsim/RandomUtils/Providers/EngineId.h"
#include "larsim/RandomUtils/Providers/Policies.h"
#include "larsim/RandomUtils/Providers/EventSeedInputData.h"

// more headers included in the implementation section below


// The master list of all the policies
#define SEED_SERVICE_POLICIES                   \
  SEED_SERVICE_POLICY(unDefined)                \
  SEED_SERVICE_POLICY(autoIncrement)            \
  SEED_SERVICE_POLICY(linearMapping)            \
  SEED_SERVICE_POLICY(preDefinedOffset)         \
  SEED_SERVICE_POLICY(preDefinedSeed)           \
  SEED_SERVICE_POLICY(random)                   \
  SEED_SERVICE_POLICY(perEvent)                 \
  /**/


namespace sim {
  
  /**
   * @brief A class to assist in the distribution of guaranteed unique seeds to all engine IDs.
   * @tparam SEED type of random engine seed
   * @see LArSeedService
   *
   * @attention direct use of this class is limited to art-less contexts;
   * within art, use `sim::LArSeedService` instead.
   * 
   * This class is configured from a FHiCL parameter set.
   * The complete configuration depends on the policy chosen; the following
   * parameters are common to all the policies:
   *     
   *     LArSeedService : {
   *        policy           : "autoIncrement" // Required: Other legal value are listed in SEED_SERVICE_POLICIES
   *        verbosity        : 0               // Optional: default=0, no informational printout
   *        endOfJobSummary  : false           // Optional: print list of all managed seeds at end of job.
   *     }
   *     
   * The policy parameter tells the service to which algorithm to use.
   * If the value of the policy parameter is not one of the known policies, the code will
   * throw an exception.
   *
   * Code instanciating a SeedMaster can request a seed by making one of the
   * following two calls:
   *     
   *     SeedMasterInstance.getSeed("moduleLabel");
   *     SeedMasterInstance.getSeed("moduleLabel", "instanceName");
   *     
   * art module code requests a seed directly to LArSeedService service.
   *
   * It is the callers responsibility to use the appropriate form.
   *
   * When `getSeed` is called with a particular module label and instance name, it computes
   * a seed value, saves it and returns it.  If there is a subsequent call to `getSeed` with
   * the same module label and instance name, the class will return the saved value of
   * the seed.  The following text will use the phrase "unique calls to `getSeed`"; two
   * calls with the same module label and instance names are not considered unique.
   *
   * If the policy is defined as `autoIncrement`, the additional configurable
   * items are:
   *     
   *     LArSeedService : {
   *        policy           : "autoIncrement"
   *        // ... and all the common ones, plus:
   *        baseSeed         : 0     // Required: An integer >= 0.
   *        checkRange       : true  // Optional: legal values true (default) or false
   *        maxUniqueEngines : 20    // Required iff checkRange is true.
   *     }
   *     
   * In this policy, the seed is set to `baseSeed`+`offset`, where
   * on the first unique call to `getSeed` the `offset` is set to 0; on the
   * second unique call to `getSeed` it is set to 1, and so on.
   *
   * If the policy is defined as `linearMapping`, the additional configurable
   * items are:
   *     
   *     LArSeedService : {
   *        policy           : "linearMapping"
   *        // ... and all the common ones, plus:
   *        nJob             : 0     // Required: An integer >= 0.
   *        checkRange       : true  // Optional: legal values true (default) or false
   *        maxUniqueEngines : 20    // Required iff checkRange is true.
   *     }
   *     
   * In this policy, the seed is set to
   * `maxUniqueEngines`*`nJob`+`offset`, where on the first unique call to `getSeed`
   * the `offset` is set to 0; on the second unique call to `getSeed` it is set to 1,
   * and so on.
   *
   * If the policy is defined as `preDefinedOffset`, the additional configurable
   * items are:
   *     
   *     LArSeedService : {
   *        policy           : "preDefinedOffset"
   *        // ... and all the common ones, plus:
   *        baseSeed         : 0     // Required: An integer >= 0.
   *        checkRange       : true  // Optional: legal values true (default) or false
   *        maxUniqueEngines : 20    // Required iff checkRange is true
   *        
   *        module_label1: offset1      // for each module with a nameless engine
   *        module_label2: {            // for each module with nemed engine instances
   *          instance_name1: offset21  //   ... one entry for each instance name
   *          instance_name2: offset22
   *        }
   *     }
   *     
   * In this policy, when `getSeed` is called, the class will look
   * into the parameter set to find a defined offset for the specified module label
   * and instance name.  The returned value of the seed will be `baseSeed`+`offset`.
   *
   * If the policy is defined as `preDefinedSeed`, the additional configurable
   * items are:
   *     
   *     LArSeedService : {
   *        policy           : "preDefinedSeed"
   *        // ... and all the common ones, plus:
   *        
   *        module_label1: seed1      // for each module with a nameless engine
   *        module_label2: {          // for each module with nemed engine instances
   *          instance_name1: seed21  //   ... one entry for each instance name
   *          instance_name2: seed22
   *        }
   *     }
   *     
   * This policy allows to specify the actual seed to be used.
   * Note that the policy
   * does not impose any constraint on the user-provided set of seeds.
   * In particular, the uniqueness of the seeds is not
   * enforced. Intended for debugging and special tests, use with care.
   *
   * If the policy is defined as `random`, the additional configurable
   * items are:
   *     
   *     LArSeedService : {
   *        policy           : "random"
   *        // ... and all the common ones, plus:
   *        masterSeed: master_seed // optional: an integer >= 0
   *     }
   *     
   * With this policy, the seed is extracted from a local random number
   * generator.
   * The seed used to initialize this additional random number generatot is
   * taken from the clock, unless the `masterSeed` parameter is set to specify
   * the actual seed.
   * This policy is meant as a quick way to disentangle the code from the
   * random seed policy used, and it's meant for special needs only and
   * definitely not for production.
   * You can enable this policy instead of whatever is already in your
   * configuration by adding at the end of your configuration:
   *     
   *     services.LArSeedService.policy: "random"
   *     
   * (this assumes that the configuration of the SeedMaster is read from
   * `services.LArSeedService`, that is the case in the art framework).
   * 
   * The FHiCL grammar to specify the offsets takes two forms.  If no instance name
   * is given, the offset is given by:
   *     
   *     moduleLabel : offset
   *     
   * When a module has multiple instances, the offsets are given by:
   *     
   *     moduleLabel : {
   *        instanceName1 : offset1
   *        instanceName2 : offset2
   *     }
   *     
   * `SeedMaster` does several additional checks, except for the `preDefinedSeed` policy.
   *
   * If one (module label, instance name) has the same seed as another (module label, instance name),
   * the class will throw an exception.
   *
   * If the `checkRange` parameter is set to `true`, and if an offset is generated with a value outside
   * the allowed range (typically `0<= offset < maxUniqueEngines-1`) then the class will also throw an exception.
   *
   * It is the responsibility of the user to ensure that the parameters
   * (e.g. `nJob` and `maxUniqueEngines`) are chosen
   * it a way that ensures the required level of uniqueness of seeds.  The example grid jobs have
   * a single point of maintenance to achieve this: the user must specify the starting job number
   * for each grid submission.
   */
  template <typename SEED>
  class SeedMaster {
      public:
    using seed_t = SEED;                   ///< type of served seeds
    using SeedMaster_t = SeedMaster<SEED>; ///< type of this class
      
    using EngineId = SeedMasterHelper::EngineId; ///< type of engine ID
    
    /// type of a function setting a seed
    using Seeder_t = std::function<void(EngineId const&, seed_t)>;
    
      private:
    /// Type of abstract class for policy implementation
    using PolicyImpl_t = details::RandomSeedPolicyBase<seed_t>;
    
    /// Type for seed data base
    using map_type = std::map<EngineId, seed_t>;
    
    /// type of map of seeders associated with the engines
    using Seeders_t = std::map<EngineId, Seeder_t>;
    
      public:
    /// type of data used for event seeds
    using EventData_t = typename PolicyImpl_t::EventData_t;
    
    /// An invalid seed
    static constexpr seed_t InvalidSeed = PolicyImpl_t::InvalidSeed;
    
    enum Policy {
#define SEED_SERVICE_POLICY(x) x,
      SEED_SERVICE_POLICIES
#undef SEED_SERVICE_POLICY
    };
    
    static const std::vector<std::string>& policyNames();
    
    /// An iterator to the configured engine IDs
    using EnginesIteratorBox
      = SeedServiceHelper::MapKeyConstIteratorBox<Seeders_t>;
    
    SeedMaster(const fhicl::ParameterSet&);
    
    // Accept compiler written c'tors, assignments and d'tor.
    
    /// Returns whether the specified engine is already registered
    bool hasEngine(EngineId const& id) const { return seeders.count(id) > 0; }
    
    /// Returns the seed value for this module label
    seed_t getSeed(std::string moduleLabel);
    
    /// Returns the seed value for this module label and instance name
    seed_t getSeed(std::string moduleLabel, std::string instanceName);
    
    /// Returns the seed value for the engine with the specified ID
    seed_t getSeed(EngineId const&);
    
    //@{
    /// Returns the seed value for the event with specified data
    seed_t getEventSeed(EventData_t const& data, std::string instanceName);
    seed_t getEventSeed(EventData_t const& data, EngineId const& id);
    //@}
    
    
    /**
     * @brief Register the specified function to reseed the engine id
     * @param id ID of the engine to be associated to the seeder
     * @param seeder function to be used for seeding the engine
     *
     * SeedMaster keeps a list of functions that can be used to reseed an
     * existing engine.
     * When reseedEvent() (or reseed()) is called, these functions are invoked
     * to set the seed to the engine.
     */
    void registerSeeder(EngineId const& id, Seeder_t seeder);
    
    
    /// Forces SeedMaster not to change the seed of a registered engine
    void freezeSeed(EngineId const& id) { seeders.at(id) = nullptr; }
    
    
    /// Reseeds the specified engine with a global seed (if any)
    seed_t reseed(EngineId const& id);
    
    /// Reseeds the specified engine with an event seed (if any)
    seed_t reseedEvent(EngineId const& id, EventData_t const& data);
    
    /// Prints known (EngineId,seed) pairs
    template<typename Stream> void print(Stream&&) const;
    
    /// Returns an object to iterate in range-for through configured engine IDs
    EnginesIteratorBox engineIDsRange() const { return { seeders }; }
    
    /// Prepares for a new event
    void onNewEvent();
    
    /// Prints to the framework Info logger
    void print() const { print(mf::LogVerbatim("SEEDS")); }    
    
      private:
    /// Control the level of information messages
    int verbosity;
    
    /// Which of the supported policies to use?
    Policy policy;
    
    /// List of seeds already computed. Seeds are unique per EngineID
    map_type knownSeeds;
    
    /// List of event seeds already computed. Seeds are unique per EngineID
    map_type knownEventSeeds;
    
    /// List of seeder functors
    Seeders_t seeders;

    /// Helper function to parse the policy name
    void setPolicy(std::string policyName);
    
    /// Throws if the seed has already been used
    void ensureUnique
      (EngineId const& id, seed_t seed, map_type const& map) const;
    void ensureUnique(EngineId const& id, seed_t seed) const
      { return ensureUnique(id, seed, knownSeeds); }
    
    /// the instance of the random policy
    std::unique_ptr<PolicyImpl_t> policy_impl;
    
  }; // class SeedMaster
  
} // namespace sim


//==============================================================================
//===  Template implementation
//===


// C++ include files
#include <ostream>
#include <iomanip> // std::setw()
#include <ostream> // std::endl
#include <algorithm> // std::find(), std::copy()
#include <iterator> // std::ostream_iterator<>, std::distance()

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art include files
#include "art/Utilities/Exception.h"


//----------------------------------------------------------------------------
template <typename SEED>
std::vector<std::string> const& sim::SeedMaster<SEED>::policyNames() {
  static std::vector<std::string> names;
  if(names.empty()) {
    const char *cnames[] = {
#define SEED_SERVICE_POLICY(x) #x,
      SEED_SERVICE_POLICIES
#undef SEED_SERVICE_POLICY
    };
    names = std::vector<std::string>
      (cnames, cnames + sizeof(cnames)/sizeof(cnames[0]));
  }

  return names;
} // SeedMaster<SEED>::policyNames()



//----------------------------------------------------------------------------
template <typename SEED>
sim::SeedMaster<SEED>::SeedMaster(fhicl::ParameterSet const& pSet):
  verbosity(pSet.get<int>("verbosity",0)),
  policy(unDefined),
  knownSeeds(),
  knownEventSeeds(),
  seeders()
{
  
  // Throw if policy is not recognized.
  const std::string strPolicy = pSet.get<std::string>("policy");
  setPolicy(strPolicy);
  
  // Finish parsing the parameter set, as required by the selected policy
  switch(policy) {
  case autoIncrement:
    policy_impl.reset(new details::AutoIncrementPolicy<seed_t>(pSet));
    break;
  case linearMapping:
    policy_impl.reset(new details::LinearMappingPolicy<seed_t>(pSet));
    break;
  case preDefinedOffset:
    policy_impl.reset(new details::PredefinedOffsetPolicy<seed_t>(pSet));
    break;
  case preDefinedSeed:
    policy_impl.reset(new details::PredefinedSeedPolicy<seed_t>(pSet));
    break;
  case random:
    policy_impl.reset(new details::RandomPolicy<seed_t>(pSet));
    break;
  case perEvent:
    policy_impl.reset(new details::PerEventPolicy<seed_t>(pSet));
    break;
  default:
    throw art::Exception(art::errors::LogicError)
      << "SeedMaster(): Internal error: unknown policy_ value";
  } // switch
  
  if ( verbosity > 0 )
    print(mf::LogVerbatim("SeedMaster"));
  
} // SeedMaster<SEED>::SeedMaster()



//----------------------------------------------------------------------------
template <typename SEED>
inline typename sim::SeedMaster<SEED>::seed_t
sim::SeedMaster<SEED>::getSeed
  (std::string moduleLabel)
{
  return getSeed(EngineId(moduleLabel));
} // SeedMaster<SEED>::getSeed(string)


//----------------------------------------------------------------------------
template <typename SEED>
typename sim::SeedMaster<SEED>::seed_t sim::SeedMaster<SEED>::getSeed
  (std::string moduleLabel, std::string instanceName)
{
  return getSeed(EngineId(moduleLabel,instanceName));
} // SeedMaster<SEED>::getSeed(string, string)


//----------------------------------------------------------------------------
template <typename SEED>
void sim::SeedMaster<SEED>::registerSeeder
  (EngineId const& id, Seeder_t seeder)
{
  seeders[id] = seeder;
} // SeedMaster<SEED>::registerSeeder()


//----------------------------------------------------------------------------
template <typename SEED>
typename sim::SeedMaster<SEED>::seed_t sim::SeedMaster<SEED>::reseed
  (EngineId const& id)
{
  seed_t seed = getSeed(id);
  if (seed != InvalidSeed) { // reseed
    Seeder_t seeder = seeders[id]; // it might create a new (empty) entry
    if (seeder) seeder(id, seed);
  }
  return seed;
} // SeedMaster<SEED>::reseed()


template <typename SEED>
typename sim::SeedMaster<SEED>::seed_t sim::SeedMaster<SEED>::reseedEvent
  (EngineId const& id, EventData_t const& data)
{
  seed_t seed = getEventSeed(data, id);
  if (seed != InvalidSeed) { // reseed
    Seeder_t seeder = seeders[id]; // it might create a new (empty) entry
    if (seeder) seeder(id, seed);
  }
  return seed;
} // SeedMaster<SEED>::reseedEvent()



//----------------------------------------------------------------------------
template <typename SEED> template <typename Stream>
void sim::SeedMaster<SEED>::print(Stream&& log) const {
  log << "\nSummary of seeds computed by the LArSeedService";
  
  // allow the policy implementation to print whatever it feels to
  std::ostringstream sstr;
  policy_impl->print(sstr);
  if (!sstr.str().empty()) log << '\n' << sstr.str();
  
  if ( !knownSeeds.empty() ) {
    log << "\n Seed Value     ModuleLabel.InstanceName";
    for (auto const& p: knownSeeds)
      log << "\n" << std::setw(10) << p.second << "      " << p.first;
  } // if
  log << '\n' << std::endl;
} // SeedMaster<SEED>::print(Stream)


//----------------------------------------------------------------------------
template <typename SEED>
typename sim::SeedMaster<SEED>::seed_t sim::SeedMaster<SEED>::getSeed
  (EngineId const& id)
{
  // Check for an already computed seed.
  typename map_type::iterator iSeed = knownSeeds.find(id);
  if (iSeed != knownSeeds.end()) return iSeed->second;
  
  // Compute the seed.
  seed_t seed = policy_impl->getSeed(id);
  if (policy_impl->yieldsUniqueSeeds()) ensureUnique(id, seed);
  
  // Save the result.
  knownSeeds.emplace(id, seed);
  
  return seed;
} // SeedMaster<SEED>::getSeed()


//----------------------------------------------------------------------------
template <typename SEED>
typename sim::SeedMaster<SEED>::seed_t sim::SeedMaster<SEED>::getEventSeed
  (EventData_t const& data, EngineId const& id)
{
  // Check for an already computed seed.
  typename map_type::iterator iSeed = knownEventSeeds.find(id);
  if (iSeed != knownEventSeeds.end()) return iSeed->second;
  
  // Compute the seed.
  seed_t seed = policy_impl->getEventSeed(id, data);
  if ((seed != InvalidSeed) && policy_impl->yieldsUniqueSeeds())
    ensureUnique(id, seed, knownEventSeeds);
  
  // Save the result.
  knownEventSeeds.emplace(id, seed);
  
  return seed;
} // SeedMaster<SEED>::getEventSeed(EngineId)


template <typename SEED>
typename sim::SeedMaster<SEED>::seed_t sim::SeedMaster<SEED>::getEventSeed
  (EventData_t const& data, std::string instanceName)
{
  return getEventSeed(data, EngineId(data.moduleLabel, instanceName));
} // SeedMaster<SEED>::getEventSeed(string)



//----------------------------------------------------------------------------
template <typename SEED>
inline void sim::SeedMaster<SEED>::onNewEvent() {
  // forget all we know about the current event
  knownEventSeeds.clear();
} // SeedMaster<SEED>::onNewEvent()


//----------------------------------------------------------------------------
template <typename SEED>
void sim::SeedMaster<SEED>::setPolicy(std::string policyName) {
  
  std::vector<std::string>::const_iterator iter
    = std::find(policyNames().begin(), policyNames().end(), policyName);
  if (iter != policyNames().end()) {
    policy = Policy(std::distance(policyNames().begin(), iter));
  }
  
  if (policy == unDefined) {
    std::ostringstream os;
    os<< "LArSeedService::setPolicy(): Unrecognized policy: "
      << policyName
      << "\n Known policies are: ";
    
    std::copy(policyNames().begin(), policyNames().end(),
      std::ostream_iterator<std::string>(os, ", "));
    
    throw art::Exception(art::errors::Configuration) << os.str();
  }
} // SeedMaster<SEED>::setPolicy()


//----------------------------------------------------------------------------
template <typename SEED>
void sim::SeedMaster<SEED>::ensureUnique
  (EngineId const& id, seed_t seed, map_type const& seeds) const
{
  
  for (auto const& p: seeds) {
    
    // Do not compare to self
    if ( p.first == id ) continue;
    
    if ( p.second == seed ){
      throw art::Exception(art::errors::LogicError)
        << "LArSeedService::ensureUnique() seed: "<<seed
        << " already used by module.instance: " << p.first << "\n"
        << "May not be reused by module.instance: " << id;
    }
  } // for
} // SeedMaster<SEED>::ensureUnique()


#endif // SeedService_SeedMaster_h
