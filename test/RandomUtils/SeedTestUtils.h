/**
 * @file   SeedTestUtils.h
 * @brief  Functions used in LArSeedService tests
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 17, 2016
 * @see    SeedTestUtils.cxx
 */

#ifndef ARTEXTENSIONS_TEST_SEEDSERVICE_SEEDTESTUTILS_H
#define ARTEXTENSIONS_TEST_SEEDSERVICE_SEEDTESTUTILS_H 1

// C++ includes.
#include <exception>
#include <string>

// Framework includes.
#include "art/Utilities/Exception.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// forward declarations
namespace CLHEP {
  class HepRandomEngine;
  class RandFlat;
} // namespace CLHEP


namespace testing {
  namespace LArSeedService {
    
    /// Type of seed
    using seed_t = art::RandomNumberGenerator::seed_t;
    
    /// Returns whether the exception looks to be from LArSeedService
    bool isSeedServiceException(std::exception const& e);
    
    
    /// Returns whether the exception looks to be from LArSeedService
    /// and has specified code
    bool isSeedServiceException
      (std::exception const& e, art::errors::ErrorCodes code);
    
    /// Extracts "statistics" from a flat random number distribution
    unsigned short RollStat(CLHEP::RandFlat& rand);
    
    /// Creates a "character statistics" using the specified random engine
    std::string CreateCharacter(CLHEP::HepRandomEngine& engine);
    
    
    /// @{
    /// Returns the seed of the specified engine
    seed_t readSeed(CLHEP::HepRandomEngine const& engine);
    /// @}
    
    
  } // namespace LArSeedService
} // namespace testing
  
#endif // ARTEXTENSIONS_TEST_SEEDSERVICE_SEEDTESTUTILS_H
