/**
 * @file   StandardPolicies.h
 * @brief  Implementation of some standard and historical random seed assignment
 *         policies
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   20150211
 * @see    SeedMaster.h
 * 
 * No code in this files is directly serviceable.
 * Documentation is up to date though.
 */

#ifndef LARSIM_RANDOMUTILS_PROVIDERS_STANDARDPOLICIES_H
#define LARSIM_RANDOMUTILS_PROVIDERS_STANDARDPOLICIES_H 1

// C/C++ standard libraries
#include <string>
#include <ostream> // std::endl

// From art and its tool chain
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// Some helper classes
#include "larsim/RandomUtils/Providers/RandomSeedPolicyBase.h"
#include "larsim/RandomUtils/Providers/EngineId.h"


namespace sim {
  
  namespace details {
    
    /** ************************************************************************
     * @brief Implementation of the "autoIncrement" policy
     * @see CheckedRangePolicy
     *
     * This is heavily based on CheckedRangePolicy.
     */
    template <typename SEED>
    class AutoIncrementPolicy: public CheckedRangePolicy<SEED> {
        public:
      using base_t = CheckedRangePolicy<SEED>;
      using this_t = AutoIncrementPolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Configures from a parameter set
      /// @see configure()
      AutoIncrementPolicy(fhicl::ParameterSet const& pset):
        base_t("autoIncrement")
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
        {
          base_t::range_check.SetConfigLabels
            ("maxUniqueEngines", "", "checkRange");
          base_t::configure(pset);
          static_configure(pset);
        } // configure()
      
      /// Prints the configuration of this policy
      virtual void print(std::ostream& out) const override;
      
        protected:
      seed_t first_seed;
      seed_t next_seed; ///< next seed delivered
      
      /// Returns the next random number
      virtual seed_t createSeed(SeedMasterHelper::EngineId const&) override
        { return next_seed++; }
      
      void static_configure(fhicl::ParameterSet const&);
    }; // class AutoIncrementPolicy<>
    
    
    template <typename SEED>
    void AutoIncrementPolicy<SEED>::static_configure
      (fhicl::ParameterSet const& pset)
    {
      first_seed = pset.get<seed_t>("baseSeed");
      base_t::range_check.SetBaseSeed(first_seed);
      base_t::CheckRangeConfiguration();
      next_seed = first_seed;
    } // AutoIncrementPolicy<SEED>::configure()
    
    
    template <typename SEED>
    void AutoIncrementPolicy<SEED>::print(std::ostream& out) const {
      base_t::print(out);
      out << "\n  first seed: " << first_seed;
    } // AutoIncrementPolicy<SEED>::print()
    
    
    /** ************************************************************************
     * @brief Implementation of the "linearMapping" policy
     * @see CheckedRangePolicy
     *
     * This is heavily based on CheckedRangePolicy.
     */
    template <typename SEED>
    class LinearMappingPolicy: public CheckedRangePolicy<SEED> {
        public:
      using base_t = CheckedRangePolicy<SEED>;
      using this_t = LinearMappingPolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Configures from a parameter set
      /// @see configure()
      LinearMappingPolicy(fhicl::ParameterSet const& pset):
        base_t("linearMapping")
        { this_t::configure(pset); }
      
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       * Parameters:
       * - *nJob* (unsigned integer): the number of this job; the first seed
       * - *checkRange* (boolean, default: true): whether to verify that each
       *   seed is within the expected range
       * - *maxUniqueEngines* (unsigned integer, mandatory) the maximum number
       *   on seeds we expect to create
       */
      virtual void configure(fhicl::ParameterSet const& pset) override
        {
          base_t::range_check.SetConfigLabels("", "", "checkRange");
          base_t::configure(pset);
          static_configure(pset);
        }
      
      /// Prints the configuration of this policy
      virtual void print(std::ostream& out) const override;
      
        protected:
      seed_t first_seed; ///< base seed
      seed_t next_seed; ///< next seed delivered
      unsigned int nSeedsPerJob;
      
      /// Returns the next random number
      virtual seed_t createSeed(SeedMasterHelper::EngineId const&) override
        { return next_seed++; }
      
      void static_configure(fhicl::ParameterSet const& pset);
      
    }; // class LinearMappingPolicy<>
    
    
    template <typename SEED>
    void LinearMappingPolicy<SEED>::static_configure
      (fhicl::ParameterSet const& pset)
    {
      // this code is for legacy support, and it could disappear in the future
      if (!pset.get_if_present<seed_t>("nJob", first_seed)) {
        if (!pset.get_if_present<seed_t>("baseSeed", first_seed)) {
          // this is going to fail; I am doing this just to get
          // the more appropriate error message possible
          first_seed = pset.get<seed_t>("nJob");
        }
        else {
          mf::LogWarning("SeedMaster") <<
            std::string(80, '*') <<
            "\nDEPRECATION WARNING: 'baseSeed' parameter has been deprecated"
              " for linearMapping policy, in favour of 'nJob'."
            "\nPlease update your configuration accordingly."
            << "\n" << std::string(80, '*');
        }
      }
    //  first_seed = pset.get<seed_t>("nJob");
      nSeedsPerJob = pset.get<seed_t>("maxUniqueEngines");
      first_seed *= nSeedsPerJob;
      ++first_seed; // we don't want 0 as a seed
      next_seed = first_seed;
      base_t::range_check.SetBaseSeed(next_seed);
      base_t::range_check.SetNSeeds(nSeedsPerJob);
      base_t::CheckRangeConfiguration();
    } // LinearMappingPolicy<SEED>::configure()
    
    
    template <typename SEED>
    void LinearMappingPolicy<SEED>::print(std::ostream& out) const {
      base_t::print(out);
      out
        << "\n  first seed:    " << first_seed
        << "\n  seeds per job: " << nSeedsPerJob;
    } // LinearMappingPolicy<SEED>::print()
    
    
    
    /** ************************************************************************
     * @brief Implementation of the "preDefinedSeed" policy
     *
     */
    template <typename SEED>
    class PredefinedSeedPolicy: public PerInstancePolicy<SEED> {
        public:
      using base_t = PerInstancePolicy<SEED>;
      using this_t = PredefinedSeedPolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Configures from a parameter set
      /// @see configure()
      PredefinedSeedPolicy(fhicl::ParameterSet const& pset):
        base_t("preDefinedSeed")
        { this_t::configure(pset); }
      
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       * Parameters: one entry per engine.
       * The FHiCL grammar to specify the seeds takes two forms.
       * If no instance name is given, the seed is given by:
       *
       *     moduleLabel : seed
       *
       * When a module has multiple instances, the seeds are given by:
       *
       *     moduleLabel : {
       *        instanceName1 : seed1
       *        instanceName2 : seed2
       *     }
       * 
       */
      virtual void configure(fhicl::ParameterSet const& pset) override
        { base_t::configure(pset); static_configure(pset); }
      
      /// Prints the configuration of this policy
      virtual void print(std::ostream& out) const override;
      
      
      /// Returns whether the returned seed should be unique: for us it "no".
      virtual bool yieldsUniqueSeeds() const { return false; }
      
        protected:
      
      /// Returns the seed stored in the parameter set
      virtual seed_t createSeed(SeedMasterHelper::EngineId const& id) override
        { return base_t::getInstanceSeed(id); }
      
      void static_configure(fhicl::ParameterSet const&)
        { base_t::range_check.SetCheck(false); }
      
    }; // class PredefinedSeedPolicy<>
    
    
    template <typename SEED>
    void PredefinedSeedPolicy<SEED>::print(std::ostream& out) const {
      base_t::print(out);
      out << "\n  seeds directly from the configuration";
    } // PredefinedSeedPolicy<SEED>::print()
    
    
    
    /** ************************************************************************
     * @brief Implementation of the "preDefinedOffset" policy
     * @see CheckedRangePolicy
     *
     * This is heavily based on CheckedRangePolicy.
     */
    template <typename SEED>
    class PredefinedOffsetPolicy: public PerInstancePolicy<SEED> {
        public:
      using base_t = PerInstancePolicy<SEED>;
      using this_t = PredefinedOffsetPolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Configures from a parameter set
      /// @see configure()
      PredefinedOffsetPolicy(fhicl::ParameterSet const& pset):
        base_t("preDefinedOffset")
        { this_t::configure(pset); }
      
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       * Parameters:
       * - *baseSeed* (unsigned integer): the base seed
       * - *checkRange* (boolean, default: true): whether to verify that each
       *   seed is within the expected range
       * - *maxUniqueEngines* (unsigned integer, mandatory if /checkRange/ is
       *   true) the maximum number on seeds we expect to create
       * - in addition, one entry per engine (see below)
       * 
       * The FHiCL grammar to specify the offsets takes two forms.
       * If no instance name is given, the offset is given by:
       *
       *     moduleLabel : offset
       *
       * When a module has multiple instances, the offsets are given by:
       *
       *     moduleLabel : {
       *        instanceName1 : offset1
       *        instanceName2 : offset2
       *     }
       */
      virtual void configure(fhicl::ParameterSet const& pset) override
        {
          base_t::range_check.SetConfigLabels
            ("maxUniqueEngines", "", "checkRange");
          base_t::configure(pset);
          static_configure(pset);
        }
      
      /// Prints the configuration of this policy
      virtual void print(std::ostream& out) const override;
      
        protected:
      seed_t base_seed;
      
      /// Returns the seed stored in the parameter set
      virtual seed_t createSeed(SeedMasterHelper::EngineId const& id) override
        { return base_seed + base_t::getInstanceSeed(id); }
      
      void static_configure(fhicl::ParameterSet const&);
      
    }; // class PredefinedOffsetPolicy<>
    
    
    template <typename SEED>
    void PredefinedOffsetPolicy<SEED>::static_configure
      (fhicl::ParameterSet const& pset)
    {
      base_seed = pset.get<seed_t>("baseSeed");
      base_t::range_check.SetBaseSeed(base_seed);
      base_t::CheckRangeConfiguration();
    } // PredefinedOffsetPolicy<SEED>::configure()
    
    
    template <typename SEED>
    void PredefinedOffsetPolicy<SEED>::print(std::ostream& out) const {
      base_t::print(out);
      out << "\n  base seed: " << base_seed;
    } // PredefinedOffsetPolicy<SEED>::print()
    
    
  } // namespace details
  
} // namespace sim


#endif // LARSIM_RANDOMUTILS_PROVIDERS_STANDARDPOLICIES_H
