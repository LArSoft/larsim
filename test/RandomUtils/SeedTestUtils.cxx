/**
 * @file   SeedTestUtils.cxx
 * @brief  Functions used in NuRandomService tests (implementation file)
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 17, 2016
 * @see    SeedTestUtils.h
 */

// library header
#include "SeedTestUtils.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// C++ includes.
#include <string>
#include <sstream>
#include <iomanip> // std::setw()

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine
#include "CLHEP/Random/RandFlat.h"

// Framework includes.
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace testing {
  namespace NuRandomService {
    
    /// Returns whether the exception looks to be from NuRandomService
    bool isSeedServiceException(std::exception const& e) {
      
      // is it an art::Exception? (NuRandomService exceptions are)
      art::Exception const* art_exc = dynamic_cast<art::Exception const*>(&e);
      if (!art_exc) return false;
      
      // All configuration exceptions are from NuRandomService.
      // Why? just because.
      if (art_exc->categoryCode() == art::errors::Configuration)
        return true;
      
      // does it have "NuRandomService" in the message? (NuRandomService exceptions should)
      if ( (art_exc->explain_self().find("NuRandomService") == std::string::npos)
        && (art_exc->explain_self().find("SeedMaster") == std::string::npos)
        )
        return false;
      
      // we can't do better than this
      return true;
    } // isSeedServiceException(exception)
    
    
    /// Returns whether the exception looks to be from NuRandomService and has specified code
    bool isSeedServiceException
      (std::exception const& e, art::errors::ErrorCodes code)
    {
      if (!isSeedServiceException(e)) return false;
      // by now we know "e" is a art::Exception
      return static_cast<art::Exception const&>(e).categoryCode() == code;
    } // isSeedServiceException(exception, ErrorCodes)
    
    
    /// Extracts "statistics" from a flat random number distribution
    unsigned short RollStat(CLHEP::RandFlat& rand) {
    //  std::array<unsigned short, 4> rolls;
      unsigned int min = 6, total = 0;
      for (int i = 0; i < 4; ++i) {
        unsigned int roll = 1 + rand.fireInt(0, 6);
      //  rolls[i] = roll;
        if (min > roll) min = roll;
        total += roll;
      }
      total -= min;
      return total;
    } // RollStat()
    
    /// Creates a "character statistics" using the specified random engine
    std::string CreateCharacter(CLHEP::HepRandomEngine& engine) {
      CLHEP::RandFlat flat(engine);
      constexpr size_t NStats = 6;
      static const std::array<std::string, NStats> statNames
        = { "STR", "DEX", "CON", "INT", "WIS", "CHA" };
      std::array<unsigned short int, NStats> stats;
      std::generate(stats.begin(), stats.end(), [&flat]{ return RollStat(flat); });
      
      short int bonus = 0;
      std::ostringstream sstr;
      for (size_t iStat = 0; iStat < NStats; ++iStat) {
        sstr << "  " << statNames[iStat] << "=" << std::setw(2) << stats[iStat];
        bonus += stats[iStat] / 2 - 5;
      } // for
      sstr << "  [bonus: " << std::setw(3) << std::showpos << bonus << "]";
      return sstr.str();
    } // CreateCharacter()
    
    
    /// Returns the seed of the specified engine (CLHEP overload)
    seed_t readSeed(CLHEP::HepRandomEngine const& engine)
      { return engine.getSeed(); }
    
    
  } // namespace NuRandomService
} // namespace testing
  
