/**
 * @file   SeedTestUtils.h
 * @brief  Functions used in NuRandomService tests
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
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// forward declarations
namespace CLHEP {
  class HepRandomEngine;
  class RandFlat;
} // namespace CLHEP


namespace testing {
  namespace NuRandomService {
    
    /// Type of seed
    using seed_t = art::RandomNumberGenerator::seed_t;
    
    /// Returns whether the exception looks to be from NuRandomService
    bool isSeedServiceException(std::exception const& e);
    
    
    /// Returns whether the exception looks to be from NuRandomService
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
    
    
  } // namespace NuRandomService
} // namespace testing
  
#endif // ARTEXTENSIONS_TEST_SEEDSERVICE_SEEDTESTUTILS_H
