/**
 * @file   RandomSeedPolicyBase.h
 * @brief  Defines an interface for random seed assignment policies
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   20150211
 * @see    SeedMaster.h
 * 
 * No code in this files is directly serviceable.
 * Documentation is up to date though.
 * 
 * Also a few additional feature extensions are implemented here.
 */

#ifndef LARSIM_RANDOMUTILS_PROVIDERS_RANDOMSEEDPOLICYBASE_H
#define LARSIM_RANDOMUTILS_PROVIDERS_RANDOMSEEDPOLICYBASE_H 1

// C/C++ standard libraries
#include <vector>
#include <array>
#include <string>
#include <bitset>
#include <algorithm> // std::find()
#include <sstream>
#include <ostream> // std::endl

// From art and its tool chain
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// Art include files
#include "canvas/Utilities/Exception.h"

// Some helper classes
#include "larsim/RandomUtils/Providers/EngineId.h"
#include "larsim/RandomUtils/Providers/EventSeedInputData.h"


namespace sim {
  
  namespace details {
    
    //--------------------------------------------------------------------------
    /// Class storing a seed in the valid range
    template <typename SEED>
    class ValidSeed {
        public:
      using seed_t = SEED; ///< type of random seed
      
      static constexpr seed_t Min =         1; ///< Smallest allowed seed
      static constexpr seed_t Max = 900000000; ///< Largest allowed seed
      
      /// Forces the specified value into the allowed seed range
      template <typename T>
      static constexpr seed_t MakeValid(T s)
        { return BounceUp(Min + seed_t(s) % (Max - Min + 1), Min, Max); }
      
      /// Returns whether the specified seed is valid
      template <typename T>
      static constexpr bool isValid(T s) { return (s >= Min) && (s <= Max); }
      
      
      /// Constructor: converts from a value
      template <typename T>
      ValidSeed(T s): seed(MakeValid(s)) {}
      
      ValidSeed(const ValidSeed&) = delete;
      ValidSeed(ValidSeed&&) = default;
      ValidSeed& operator= (const ValidSeed&) = delete;
      ValidSeed& operator= (ValidSeed&&) = default;
      
      /// Return the converted seeda
      operator seed_t() const { return seed; }
      
        protected:
      seed_t seed; ///< the converted seed
      
        private:
      /// May be merged into MakeValid() with C++14
      template <typename T>
      static constexpr seed_t BounceUp(T s, T min, T max)
        { return (s < min)? s + (max - min + 1): s; }

    }; // class ValidSeed
    
    //--------------------------------------------------------------------------
    /// Interface for a policy implementation
    template <typename SEED>
    class RandomSeedPolicyBase {
        public:
      using seed_t = SEED; ///< type of the random seed
      
      /// type of data used for event seeds
      using EventData_t = SeedServiceHelper::EventSeedInputData;
      
      /// An invalid seed
      static constexpr seed_t InvalidSeed = 0;
      
      /// Constructor; requires the policy name
      RandomSeedPolicyBase(std::string policy_name):
        name(policy_name) {}
      
      // Virtual destructor
      virtual ~RandomSeedPolicyBase() {}
      
      /// Configure this policy
      virtual void configure(fhicl::ParameterSet const&) {}
      
      /// Returns the next random number
      virtual seed_t getSeed(SeedMasterHelper::EngineId const& id)
        { return createSeed(id); }
      
      /// Returns a random number specific to an event
      virtual seed_t getEventSeed
        (SeedMasterHelper::EngineId const& id, EventData_t const& eventInfo)
        { return createEventSeed(id, eventInfo); }
      
      /// Returns the given name of the policy
      std::string getName() const { return name; }
      
      /// Prints information on the configuration of this policy
      virtual void print(std::ostream& out) const
        { out << "Random policy: '" << getName() << "'"; }
      
      /// Returns whether the returned seed should be unique
      virtual bool yieldsUniqueSeeds() const { return true; }
      
        protected:
      std::string name; ///< name of the policy
      
      /// Extracts the next random number seed
      virtual seed_t createSeed(SeedMasterHelper::EngineId const&) = 0;
      
      /// Extracts a seed for specified event information; returns InvalidSeed
      virtual seed_t createEventSeed
        (SeedMasterHelper::EngineId const&, EventData_t const&)
        { return InvalidSeed; }
      
    }; // class RandomSeedPolicyBase
    
    
    
    /// Helper class to support range checking
    template <typename SEED>
    class RangeCheckHelper {
        public:
      using seed_t = SEED;
      
      /// Constructor; specify configuration labels
      RangeCheckHelper(
        std::string maxSeedsLabel = "maxUniqueEngines",
        std::string baseSeedLabel = "baseSeed",
        std::string checkRangeLabel = "checkRange"
        )
        { SetConfigLabels(maxSeedsLabel, baseSeedLabel, checkRangeLabel); }
      
      
      void SetConfigLabels(
        std::string maxSeedsLabel = "maxUniqueEngines",
        std::string baseSeedLabel = "baseSeed",
        std::string checkRangeLabel = "checkRange"
        );
      
      /// Configures from a parameter set
      /// @return if the configuration after the call is complete
      bool configure(fhicl::ParameterSet const& pset);
      
      
      /// Sets whether to perform the check or not
      void SetCheck(bool doCheck = true)
        { bCheck = doCheck; hasParameters.set(pmDoCheck); }
      
      /// Sets the base seed directly
      void SetBaseSeed(seed_t base_seed)
        { BaseSeed = base_seed; hasParameters.set(pmBaseSeed); }
      
      /// Sets the number of seeds directly
      void SetNSeeds(seed_t nSeeds)
        { MaxSeeds = nSeeds; hasParameters.set(pmMaxSeeds); }
      
      /// Performs the check on the specified seed
      bool operator() (seed_t seed) const
        {
          return
            !bCheck || ((seed >= BaseSeed) && (seed < BaseSeed + MaxSeeds));
        } // operator()
      
      /// Throws an exception if the range check on seed fails
      void EnsureRange
        (std::string policy, SeedMasterHelper::EngineId const& id, seed_t seed)
        const;

      /// Returns whether all the parameters are configured
      bool isConfigured() const;
      
      /// Returns the items currently not configured
      std::vector<std::string> missingConfig() const;
      
      /// Prints the configuration int the specified stream
      template <typename STREAM>
      void print(STREAM& out, std::string indent = std::string()) const;
        
        protected:
      typedef enum {
        pmMaxSeeds,
        pmBaseSeed,
        pmDoCheck,
        pmNParams
      } Params_t;
      static constexpr unsigned int NParams = (unsigned int) pmNParams;
      
      std::array<std::string, NParams> paramLabels;
      std::bitset<NParams> hasParameters;
      
      bool bCheck = true; ///< should we perform the check?
      seed_t BaseSeed; ///< minimum valid seed
      seed_t MaxSeeds; ///< number of valid seeds
      
    }; // class RangeCheckHelper
    
    
    template <typename SEED>
    void RangeCheckHelper<SEED>::SetConfigLabels(
      std::string maxSeedsLabel   /* = "maxUniqueEngines" */,
      std::string baseSeedLabel   /* = "baseSeed" */,
      std::string checkRangeLabel /* = "checkRange" */
    ) {
      paramLabels[pmMaxSeeds] = maxSeedsLabel;
      paramLabels[pmBaseSeed] = baseSeedLabel;
      paramLabels[pmDoCheck]  = checkRangeLabel;
    } // RangeCheckHelper::SetConfigLabels()
    
    
    template <typename SEED>
    bool RangeCheckHelper<SEED>::configure(fhicl::ParameterSet const& pset) {
      if (!paramLabels[pmMaxSeeds].empty()) {
        seed_t seed;
        if (pset.get_if_present(paramLabels[pmMaxSeeds], seed)) SetNSeeds(seed);
      }
      if (!paramLabels[pmBaseSeed].empty()) {
        seed_t seed;
        if (pset.get_if_present(paramLabels[pmBaseSeed], seed))
          SetBaseSeed(seed);
      }
      if (!paramLabels[pmDoCheck].empty()) {
        bool flag;
        if (pset.get_if_present(paramLabels[pmDoCheck], flag)) SetCheck(flag);
      }
      return isConfigured();
    } // RangeCheckHelper<SEED>::configure()
    
    
    template <typename SEED>
    bool RangeCheckHelper<SEED>::isConfigured() const {
      return (hasParameters.test(pmDoCheck) && !bCheck) || hasParameters.all();
    } // RangeCheckHelper<SEED>::isConfigured()
    
    
    template <typename SEED>
    void RangeCheckHelper<SEED>::EnsureRange(
      std::string policy,
      SeedMasterHelper::EngineId const& id, seed_t seed
    ) const {
      if (operator()(seed)) return;
      seed_t offset = seed - BaseSeed;
      throw art::Exception(art::errors::LogicError)
        << "LArSeedService (policy: " << policy << ") for engine: "
        << id << " the offset of seed " << seed << " is: " << offset << "."
        "\nAllowed seed offsets are in the range 0....(N-1) where N is: "
        << MaxSeeds << " (as configured in maxUniqueEngines)";
    } // RangeCheckHelper<SEED>::EnsureRange()
    
    
    template <typename SEED> template <typename STREAM>
    void RangeCheckHelper<SEED>::print
      (STREAM& out, std::string indent /* = "" */) const
    {
      if (!isConfigured())
        out << indent << "seed range checker not configured!";
      else if (bCheck)
        out << indent << "maximum number of seeds: " << MaxSeeds;
      else
        out << indent << "no limit on number of seeds.";
    } // RangeCheckHelper<SEED>::print()
    
    
    template <typename SEED>
    std::vector<std::string> RangeCheckHelper<SEED>::missingConfig() const {
      if (hasParameters.test(pmDoCheck) && !bCheck) return {};
      std::vector<std::string> missing;
      for (unsigned int i = 0; i < NParams; ++i)
        if (!hasParameters.test(i)) missing.push_back(paramLabels[i]);
      return missing;
    } // RangeCheckHelper<SEED>::missingConfig()
    
    
    
    /// Range-checked policy (abstract)
    template <typename SEED>
    class CheckedRangePolicy: public RandomSeedPolicyBase<SEED> {
        public:
      using base_t = RandomSeedPolicyBase<SEED>;
      using this_t = CheckedRangePolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Constructor; requires the policy name
      CheckedRangePolicy
        (std::string policy_name, fhicl::ParameterSet const& pset):
        base_t(policy_name)
        { this_t::configure(pset); }
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       * Parameters:
       * - *baseSeed* (unsigned integer): the first seed to be delivered
       * - *checkRange* (boolean, default: true): whether to verify that each
       *   seed is within the expected range
       * - *maxUniqueEngines* (unsigned integer, mandatory if /checkRange/ is
       *   true) the maximum number on seeds we expect to create
       */
      virtual void configure(fhicl::ParameterSet const& pset) override
        { base_t::configure(pset); static_configure(pset); }
      
      /// Returns the next random number
      virtual seed_t getSeed(SeedMasterHelper::EngineId const& id)
        {
          seed_t seed = this->createSeed(id);
          ensureRange(id, seed);
          return seed;
        } // getSeed()
      
      /// Prints information on the configuration of this policy
      virtual void print(std::ostream& out) const override
        {
          base_t::print(out);
          range_check.print(out << "\n", "  ");
        } // print()
      
        protected:
      RangeCheckHelper<seed_t> range_check;
      
      /// Lets the derived class to start configuration
      CheckedRangePolicy(std::string policy_name): base_t(policy_name) {}
      
      /// Performs the range checks and complains if needed
      virtual void ensureRange
        (SeedMasterHelper::EngineId const& id, seed_t seed) const
        { range_check.EnsureRange(this->getName(), id, seed); }
      
      /// Check that the configuration is complete
      void CheckRangeConfiguration() const;
      
      /// Local configuration; does not require the range config to be complete
      void static_configure(fhicl::ParameterSet const& pset);

    }; // class CheckedRangePolicy
    
    
    template <typename SEED>
    void CheckedRangePolicy<SEED>::static_configure
      (fhicl::ParameterSet const& pset)
    {
      range_check.configure(pset);
    } // CheckedRangePolicy<SEED>::static_configure()
    
    template <typename SEED>
    void CheckedRangePolicy<SEED>::CheckRangeConfiguration() const {
      if (!range_check.isConfigured()) {
        art::Exception e(art::errors::Configuration);
        e << "configuration of policy '" << this->getName()
          << "' incomplete:";
        for (std::string const& name: range_check.missingConfig())
          e << " " << name;
        throw e;
      }
    } // CheckedRangePolicy<SEED>::CheckRangeConfiguration()
    
    
    
    /** ************************************************************************
     * @brief Base class for policies reacting at engine instance level
     * @see CheckedRangePolicy
     *
     */
    template <typename SEED>
    class PerInstancePolicy: public CheckedRangePolicy<SEED> {
        public:
      using base_t = CheckedRangePolicy<SEED>;
      using this_t = PerInstancePolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Configures from a parameter set
      /// @see configure()
      PerInstancePolicy(std::string name, fhicl::ParameterSet const& pset):
        base_t(name)
        { this_t::configure(pset); }
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       */
      virtual void configure(fhicl::ParameterSet const& pset) override
        { base_t::configure(pset); static_configure(pset); }
      
        protected:
      fhicl::ParameterSet parameters; ///< configuration parameters
      
      /// Internal constructor: does not configure. For use in derived classes
      PerInstancePolicy(std::string name): base_t(name) {}
      
      /// Retrieves the parameter (seed) for the specified engine ID
      seed_t getInstanceSeed(SeedMasterHelper::EngineId const& id) const
        { return getInstanceParameter<seed_t>(parameters, id); }
      
      
      /// Retrieves the parameter (seed) for the specified engine ID
      template <typename T>
      static T getInstanceParameter
        (fhicl::ParameterSet const& pset, SeedMasterHelper::EngineId const& id);
      
      
      void static_configure(fhicl::ParameterSet const& pset);
      
    }; // class PerInstancePolicy<>
    
    
    template <typename SEED>
    void PerInstancePolicy<SEED>::static_configure
      (fhicl::ParameterSet const& pset)
    {
      parameters = pset; // copy the parameters locally
    } // PerInstancePolicy<SEED>::configure()
    
    
    // This method reads a element of type T from the instance configuration
    template <typename SEED> template <typename T>
    T PerInstancePolicy<SEED>::getInstanceParameter(
      fhicl::ParameterSet const& pset, SeedMasterHelper::EngineId const& id
    ) {
      
      // first check if the instance is actually global;
      // if so, look for it directly in the parameter set
      if (id.isGlobal()) {
        // We expect the element to be just an item.
        if (pset.is_key_to_table(id.instanceName)) {
          // this is mostly a limitation of the FHiCL syntax,
          // that we can overcome with some cumbersomeness if we need to.
          throw art::Exception(art::errors::Configuration)
            << "A seed for the global instance '" << id
            << "' was requested, but the configuration sets named instances ("
            << pset.get<fhicl::ParameterSet>(id.instanceName).to_compact_string()
            << ").\n";
        }
        T param;
        if (!pset.get_if_present(id.instanceName, param)) {
          throw art::Exception(art::errors::Configuration)
            << "LArSeedService: unable to find the parameter for global instance'"
            << id << "'\n";
        }
        return param;
      } // if global
      
      // there must be /some/ configuration for the module
      if (!pset.has_key(id.moduleLabel)) {
        throw art::Exception(art::errors::Configuration)
          << "A seed for the instance '" << id
          << "' was requested, but there is no configuration at all for '"
          << id.moduleLabel << "' module label.";
      }
      
      T param;
      if (!id.hasInstanceName()) { // Case 1: no instance name.
        // We expect the element to be just an item.
        if (pset.is_key_to_table(id.moduleLabel)) {
          // this is mostly a limitation of the FHiCL syntax,
          // that we can overcome with some cumbersomeness if we need to.
          throw art::Exception(art::errors::Configuration)
            << "A seed for the nameless instance '" << id
            << "' was requested, but the configuration sets named instances ("
            << pset.get<fhicl::ParameterSet>(id.moduleLabel).to_compact_string()
            << ").\nNameless and named engine instances can't coexist.";
        }
        if (!pset.get_if_present(id.moduleLabel, param)) {
          throw art::Exception(art::errors::Configuration)
            << "LArSeedService: unable to find the parameter for '" << id << "'";
        }
      } // if no instance name
      else { // Case 2: instance name is given.
        
        if (pset.is_key_to_atom(id.moduleLabel)) {
          // see above
          throw art::Exception(art::errors::Configuration)
            << "A seed for '" << std::string(id) << "' was requested,"
               " but the configuration sets a nameless instance of '"
            << id.moduleLabel << "'.\n"
            << "Nameless and named engine instances can't coexist.";
        }
        fhicl::ParameterSet subSet;
        if (!pset.get_if_present(id.moduleLabel, subSet)) {
          throw art::Exception(art::errors::Configuration)
            << "LArSeedService: unable to find the parameter block for: '"
            << id << "'";
        }
        
        if (!subSet.get_if_present(id.instanceName, param)) {
          throw art::Exception(art::errors::Configuration)
            << "LArSeedService: unable to find the parameter value for: '"
            << id << "'";
        }
      } // if instance name
      
      return param;
    } // PerInstancePolicy<SEED>::getInstanceParameter<>()
    
    
  } // namespace details
  
} // namespace sim


#endif // LARSIM_RANDOMUTILS_PROVIDERS_RANDOMSEEDPOLICYBASE_H
