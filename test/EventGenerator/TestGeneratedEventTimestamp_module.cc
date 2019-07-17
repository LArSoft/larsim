/**
 * @file   TestGeneratedEventTimestamp_module.cc
 * @brief  Test for GeneratedEventTimestamp
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 23, 2015
 *
 * Generated using artmod from cetpkgsupport v1_08_02.
 */

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
namespace fhicl { class ParameterSet; }

// artextensions libraries
#include "nurandom/RandomUtils/NuRandomService.h"

class TestGeneratedEventTimestamp: public art::EDAnalyzer {
    public:
  explicit TestGeneratedEventTimestamp(fhicl::ParameterSet const& pset);

private:
  void analyze(art::Event const & e) override;
  CLHEP::HepRandomEngine& fEngine;
  CLHEP::HepRandomEngine& fAuxEngine;
}; // class TestGeneratedEventTimestamp


//------------------------------------------------------------------------------
TestGeneratedEventTimestamp::TestGeneratedEventTimestamp
  (fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
  // create two random engines; obtain the random seed from NuRandomService,
  // unless overridden in configuration with key "Seed" and "AuxSeed"
  , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, pset, "Seed"))
  , fAuxEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "aux", pset, "AuxSeed"))
{}


//------------------------------------------------------------------------------
void TestGeneratedEventTimestamp::analyze(art::Event const& event)
{
  mf::LogInfo("TestGeneratedEventTimestamp")
    <<   "Event time stamp: " << event.time().value()
    << "\nRandom seeds: " << fEngine.getSeed() << " (main), "
    << fAuxEngine.getSeed() << " (aux)";
  }


DEFINE_ART_MODULE(TestGeneratedEventTimestamp)
