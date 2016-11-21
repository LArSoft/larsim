/**
 * @file   TestGeneratedEventTimestamp_module.cc
 * @brief  Test for GeneratedEventTimestamp
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 23, 2015
 *
 * Generated using artmod from cetpkgsupport v1_08_02.
 */

// C/C++ standard libraries
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <algorithm> // std::generate()
#include <iomanip> // std::setw()

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// artextensions libraries
#include "nutools/RandomUtils/NuRandomService.h"


class TestGeneratedEventTimestamp: public art::EDAnalyzer {
    public:
  explicit TestGeneratedEventTimestamp(fhicl::ParameterSet const& pset);
  
  void analyze(art::Event const & e) override;
  
    protected:
  static unsigned short RollStat(CLHEP::RandFlat& rand);
  static std::string CreateCharacter(CLHEP::HepRandomEngine& engine);
  
}; // class TestGeneratedEventTimestamp


//------------------------------------------------------------------------------
//---  TestGeneratedEventTimestamp implementation
//---
TestGeneratedEventTimestamp::TestGeneratedEventTimestamp
  (fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{
  // create two default random engines; obtain the random seed from NuRandomService,
  // unless overridden in configuration with key "Seed" and "AuxSeed"
  art::ServiceHandle<rndm::NuRandomService> Seeds;
  Seeds->createEngine(*this, pset, "Seed");
  Seeds->createEngine(*this, "HepJamesRandom", "aux", pset, "AuxSeed");
} // TestGeneratedEventTimestamp::TestGeneratedEventTimestamp()


//------------------------------------------------------------------------------
void TestGeneratedEventTimestamp::analyze(art::Event const& event) {
  
  mf::LogInfo log("TestGeneratedEventTimestamp");
  
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine const& MainEngine = rng->getEngine();
  CLHEP::HepRandomEngine const& AuxEngine = rng->getEngine("aux");
  
  rndm::NuRandomService::seed_t seed = MainEngine.getSeed(),
    aux_seed = AuxEngine.getSeed();
  
  log
    <<   "Event time stamp: " << event.time().value()
    << "\nRandom seeds: " << seed << " (main), " << aux_seed << " (aux)";
  
} // TestGeneratedEventTimestamp::analyze()


//------------------------------------------------------------------------------
unsigned short TestGeneratedEventTimestamp::RollStat(CLHEP::RandFlat& rand) {
  std::vector<unsigned short> rolls(4);
  unsigned int min = 6, total = 0;
  for (int i = 0; i < 4; ++i) {
    unsigned int roll = 1 + rand.fireInt(0, 6);
    rolls[i] = roll;
    if (min > roll) min = roll;
    total += roll;
  }
  total -= min;
  return total;
} // TestGeneratedEventTimestamp::RollStat()


//------------------------------------------------------------------------------
std::string TestGeneratedEventTimestamp::CreateCharacter
  (CLHEP::HepRandomEngine& engine)
{
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
} // TestGeneratedEventTimestamp::CreateCharacter()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(TestGeneratedEventTimestamp)

//------------------------------------------------------------------------------
