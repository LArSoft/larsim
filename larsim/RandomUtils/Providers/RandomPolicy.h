/**
 * @file   RandomPolicy.h
 * @brief  Implementation of the random seed assignment policy "random"
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   20150211
 * @see    SeedMaster.h
 * 
 * No code in this files is directly serviceable.
 * Documentation is up to date though.
 */

#ifndef SeedService_RandomPolicy_h
#define SeedService_RandomPolicy_h

// C/C++ standard libraries
#include <ostream> // std::endl
#include <memory> // std::unique_ptr<>
#include <random> // std::uniform_int_distribution, std::default_random_engine
#include <chrono> // std::system_clock

// From art and its tool chain
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// Some helper classes
#include "larsim/RandomUtils/Providers/RandomSeedPolicyBase.h"


namespace sim {
  
  namespace details {
    
    /** ************************************************************************
     * @brief Implementation of the "random" policy
     *
     * This policy extracts seeds randomly. Each seed is between 1 and the
     * maximum seed value specified at construction, both extremes included.
     * The sequence of these seeds is initialized by a random seed.
     * The random generator used is potentially low quality -- it does not
     * matter for this application.
     */
    template <typename SEED>
    class RandomPolicy: public RandomSeedPolicyBase<SEED> {
        public:
      using base_t = RandomSeedPolicyBase<SEED>;
      using this_t = RandomPolicy<SEED>;
      using seed_t = typename base_t::seed_t;
      
      /// Configures from a parameter set
      /// @see configure()
      RandomPolicy(fhicl::ParameterSet const& pset): base_t("random")
        { this_t::configure(pset); }
      
      
      /**
       * @brief Configure this policy
       * @param pset the parameter set for the configuration
       * 
       * Parameters:
       * - *masterSeed* (unsigned integer, optional): the seed of the seed
       *   generator; by default, it's taken from the system clock
       */
      virtual void configure(fhicl::ParameterSet const& pset) override;
      
      /// Prints the details of the configuration of the random generator
      virtual void print(std::ostream& out) const override;
      
      
        private:
      class RandomImpl {
          public:
        RandomImpl(seed_t master_seed, seed_t min_seed, seed_t max_seed):
          seed(master_seed),
          generator(master_seed), distribution(min_seed, max_seed)
          {}
        
        seed_t master_seed() const { return seed; }
        seed_t min() const { return distribution.min(); }
        seed_t max() const { return distribution.max(); }
        seed_t operator() () { return distribution(generator); }
        
          private:
        seed_t seed; ///< seed given at construction, for the record
        std::default_random_engine generator; ///< random engine
        std::uniform_int_distribution<seed_t> distribution; ///< flat
      }; // RandomImpl
      
      std::unique_ptr<RandomImpl> random_seed;
      
      /// Extracts a random seed
      virtual seed_t createSeed(SeedMasterHelper::EngineId const&) override
        { return (*random_seed)(); }
      
    }; // class RandomPolicy<>
    
    
    template <typename SEED>
    void RandomPolicy<SEED>::configure(fhicl::ParameterSet const& pset) {
      constexpr seed_t MagicMaxSeed = 900000000;
      seed_t master_seed = 0;
      if (!pset.get_if_present("masterSeed", master_seed)) {
        // get the base seed randomly too, from the clock,
        // and within [1; MagicMaxSeed]
        master_seed = 1 +
          std::chrono::system_clock::now().time_since_epoch().count()
          % MagicMaxSeed;
      }
      random_seed.reset(new RandomImpl(master_seed, 1, MagicMaxSeed));
    } // RandomPolicy<SEED>::configure()
    
    
    /// Prints the details of the configuration of the random generator
    template <typename SEED>
    void RandomPolicy<SEED>::print(std::ostream& out) const {
      base_t::print(out);
      out
        << "\n  master seed: " << random_seed->master_seed()
        << "\n  seed within: [ " << random_seed->min()
                << " ; " << random_seed->max() << " ]"
        ;
    } // RandomPolicy<SEED>::print()
    
    
  } // namespace details
  
} // namespace sim


#endif // SeedService_RandomPolicy_h
