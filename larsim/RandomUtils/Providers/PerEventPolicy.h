/**
 * @file   PerEventPolicy.h
 * @brief  Implementation of the per-event random seed assignment policy
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   20150211
 * @see    SeedMaster.h
 * 
 * No code in this files is directly serviceable.
 * Documentation is up to date though.
 */

#ifndef LARSIM_RANDOMUTILS_PROVIDERS_PEREVENTPOLICY_H
#define LARSIM_RANDOMUTILS_PROVIDERS_PEREVENTPOLICY_H 1

// C/C++ standard libraries
#include <string>
#include <sstream>
#include <ostream> // std::endl
#include <type_traits> // std::make_signed<>

// From art and its tool chain
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/Exception.h"

// Some helper classes
#include "larsim/RandomUtils/Providers/RandomSeedPolicyBase.h"
#include "larsim/RandomUtils/Providers/EngineId.h"
#include "larsim/RandomUtils/Providers/EventSeedInputData.h"


namespace sim {
  
  namespace details {
    
    /** ************************************************************************
     * @brief Implementation of the "perEvent" policy
     *
     * This policy extracts seeds depending on contextual information from the
     * event.
     * The information that enters the seed is the event ID
     * (run, subrun, event), the process name, and the engine ID.
     */
    template <typename SEED>
    class PerEventPolicy: public RandomSeedPolicyBase<SEED> {
        public:
      using base_t = RandomSeedPolicyBase<SEED>;
      using this_t = PerEventPolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// type for contextual event information
      using EventData_t = SeedServiceHelper::EventSeedInputData;
      
      typedef enum {
        saEventTimestamp_v1,             ///< event timestamp algorithm (v1)
        NAlgos,                          ///< total number of seed algorithms
        saUndefined,                     ///< algorithm not defined
        saDefault = saEventTimestamp_v1  ///< default algorithm
      } SeedAlgo_t; ///< seed algorithms; see reseed() documentation for details
      
      /// Configures from a parameter set
      /// @see configure()
      PerEventPolicy(fhicl::ParameterSet const& pset): base_t("perEvent")
        { this_t::configure(pset); }
      
      /// Returns whether the returned seed should be unique: for us it "no".
      virtual bool yieldsUniqueSeeds() const { return false; }
      
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       * Parameters:
       * - *algorithm* (string, default: "EventTimestamp_v1"): the name of the
       *   algorithm to get the seed
       * - *offset* (integer, optional): if specified, all the final seeds are
       *   incremented by this value; the resulting seed is not checked, it
       *   might even be invalid. This is considered an emergency hack when
       *   one absolutely needs to have a different seed than the one assigned
       *   to the event. This also defies the purpose of the policy, since after
       *   this, to reproduce the random sequences the additional knowledge of
       *   which offset was used is necessary.
       */
      virtual void configure(fhicl::ParameterSet const& pset) override;
      
      /// Prints the details of the configuration of the random generator
      virtual void print(std::ostream& out) const override;
      
      
      /// Default algorithm version
      static constexpr const char* DefaultVersion = "v1";
      
      
      /// Converts some information into a valid seed by means of hash values
      template <typename Hashable>
      static seed_t SeedFromHash(Hashable const& info)
        { return makeValid(std::hash<Hashable>()(info)); }
      
      /// Converts run, subrun and event numbers into a string
      static std::string UniqueEventIDString(EventData_t const& info);
      
      /// Converts event ID and timestamp information into a string
      static std::string UniqueEventString(EventData_t const& info);
      
      
        private:
      /// type for seed offset
      using SeedOffset_t = typename std::make_signed<seed_t>::type;
      
      SeedAlgo_t algo; ///< the algorithm to extract the seed
      
      SeedOffset_t offset; ///< offset added to all the seeds
      
      /// Per-job seed: invalid seeds are returned
      virtual seed_t createSeed(SeedMasterHelper::EngineId const&) override
        { return base_t::InvalidSeed; }
      
      /**
       * @brief Returns a seed proper for the specified event information
       * @param id random number engine ID (moule label and instance name)
       * @param info event information
       * @return a seed specific to the information provided
       *
       * The algorithm used to combine the provided information into a seed is
       * defined by the configuration. The following algorithms are supported:
       * - *EventTimestamp_v1*: includes event ID (run, subrun and event
       *   numbers), event timestamp, process name and engine ID into a hash
       *   value, used for the seed
       */
      virtual seed_t createEventSeed
        (SeedMasterHelper::EngineId const& id, EventData_t const& info)
        override;
      
      /// Renders a seed valid
      template <typename T>
      static seed_t makeValid(T value) { return ValidSeed<seed_t>(value); }
      
      /// Implementation of the EventTimestamp_v1 algorithm
      static seed_t EventTimestamp_v1
        (SeedMasterHelper::EngineId const& id, EventData_t const& info);
      
      
      //@{
      /// Algorithm name (manual) handling
      static const std::vector<std::string> algoNames;
      
      static std::vector<std::string> InitAlgoNames();
      //@}
    }; // class PerEventPolicy<>
    
    
    
    //--------------------------------------------------------------------------
    //---  PerEventPolicy template implementation
    //---
    template <typename SEED>
    std::vector<std::string> PerEventPolicy<SEED>::InitAlgoNames() {
      
      std::vector<std::string> names((size_t) NAlgos);
      
      names[saEventTimestamp_v1] = "EventTimestamp_v1";
      
      return names;
    } // PerEventPolicy<SEED>::InitAlgoNames()
    
    
    template <typename SEED>
    const std::vector<std::string> PerEventPolicy<SEED>::algoNames
      = PerEventPolicy<SEED>::InitAlgoNames();
    
    
    //--------------------------------------------------------------------------
    template <typename SEED>
    std::string PerEventPolicy<SEED>::UniqueEventIDString
      (EventData_t const& info)
    {
      std::ostringstream sstr;
      sstr << "Run: " << info.runNumber << " Subrun: " << info.subRunNumber
        << " Event: " << info.eventNumber;
      return sstr.str();
    } // PerEventPolicy<SEED>::UniqueEventIDString()
    
    
    template <typename SEED>
    std::string PerEventPolicy<SEED>::UniqueEventString
      (EventData_t const& info)
    {
      std::ostringstream sstr;
      sstr << " Timestamp: " << info.time;
      return UniqueEventIDString(info) + sstr.str();
    } // PerEventPolicy<SEED>::UniqueEventString()
    
    
    template <typename SEED>
    auto PerEventPolicy<SEED>::EventTimestamp_v1
      (SeedMasterHelper::EngineId const& id, EventData_t const& info)
      -> seed_t
    {
      if (!info.isTimeValid) {
        throw art::Exception(art::errors::InvalidNumber)
          << "Input event has an invalid timestamp,"
          " random seed per-event policy EventTimestamp_v1 can't be used.\n";
      }
      std::string s = UniqueEventString(info)
        + " Process: " + info.processName
        + " Module: " + id.moduleLabel;
      if (!id.instanceName.empty())
        s.append(" Instance: ").append(id.instanceName);
      seed_t seed = SeedFromHash(s);
      LOG_DEBUG("PerEventPolicy") << "Seed from: '" << s << "': " << seed;
      return seed;
    } // PerEventPolicy<SEED>::EventTimestamp_v1()
    
    
    //--------------------------------------------------------------------------
    template <typename SEED>
    void PerEventPolicy<SEED>::configure(fhicl::ParameterSet const& pset) {
      // set the algorithm
      algo = saUndefined;
      std::string algorithm_name
         = pset.get<std::string>("algorithm", "default");
      
      if (algorithm_name == "default") algo = saDefault;
      else {
        for (size_t iAlgo = 0; iAlgo < (size_t) NAlgos; ++iAlgo) {
          if (algorithm_name != algoNames[iAlgo]) continue;
          algo = (SeedAlgo_t) iAlgo;
          break;
        } // for
      }
      if (algo == saUndefined) {
        throw art::Exception(art::errors::Configuration)
          << "No valid event random seed algorithm specified!\n";
      }
      
      // read an optional overall offset
      offset = pset.get<SeedOffset_t>("offset", 0);
      
      // EventTimestamp_v1 does not require specific configuration
      
    } // PerEventPolicy<SEED>::configure()
    
    
    //--------------------------------------------------------------------------
    /// Prints the details of the configuration of the random generator
    template <typename SEED>
    void PerEventPolicy<SEED>::print(std::ostream& out) const {
      base_t::print(out);
      out
        << "\n  algorithm version: " << algoNames[algo];
      if (offset != 0)
        out << "\n  constant offset:   " << offset;
    } // PerEventPolicy<SEED>::print()
    
    
    //--------------------------------------------------------------------------
    template <typename SEED>
    typename PerEventPolicy<SEED>::seed_t PerEventPolicy<SEED>::createEventSeed
        (SeedMasterHelper::EngineId const& id, EventData_t const& info)
    {
      seed_t seed = base_t::InvalidSeed;
      switch (algo) {
        case saEventTimestamp_v1:
          seed = EventTimestamp_v1(id, info);
          break;
        case saUndefined:
          throw art::Exception(art::errors::Configuration)
            << "Per-event random number seeder not configured!\n";
        default:
          throw art::Exception(art::errors::LogicError)
            << "Unsupported per-event random number seeder (#"
            << ((int) algo) << ")\n";
      } // switch
      return seed + offset;
    } // PerEventPolicy<SEED>::createEventSeed()
    
    
  } // namespace details
  
} // namespace sim


#endif // LARSIM_RANDOMUTILS_PROVIDERS_PEREVENTPOLICY_H
