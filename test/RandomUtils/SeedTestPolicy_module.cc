/**
 * @file SeedTestPolicy_module.cc
 * @brief Test the LArSeedService.
 * @author Rob Kutschke (kutschke@fnal.gov), Gianluca Petrillo (petrillo@fnal.gov)
 * @see LArSeedService.hh
 */

// test library
#include "SeedTestUtils.h"

// LArSoft libraries
#define LARSIM_RANDOMUTILS_LARSEEDSERVICE_USECLHEP 1
#include "larsim/RandomUtils/LArSeedService.h"

// framework
#include "art/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/CurrentModule.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// supporting libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine
#include "CLHEP/Random/JamesRandom.h" // CLHEP::HepJamesRandom
#include "CLHEP/Random/RandFlat.h"

// C/C++ standard library
#include <string>
#include <vector>
#include <sstream>
#include <iomanip> // std::setw()
#include <memory> // std::unique_ptr<>


namespace testing {
  
  /**
   * @brief Test module for LArSeedService
   * 
   * The test writes on screen the random seeds it gets.
   * 
   * Note that the test does not actually get any random number, unless the
   * *useGenerators* option is turned on.
   * 
   * Configuration parameters:
   * - <b>instanceNames</b> (string list): use one random number generator
   *   for each instance name here specified; if not specified, a single
   *   default instance is used
   * - <b>expectedErrors</b> (unsigned integer, default: 0): expect this number
   *   of errors from LArSeedService, and throw an exception if we get a different
   *   number
   * - <b>useGenerators</b> (boolean, default: true): uses
   *   art RandomGeneratorService with the seeds from LArSeedService
   * - <b>perEventSeeds</b> (boolean, default: false): set it to true if the
   *   selected policy gives per-event seeds; in this case, the check of
   *   seed always being the same is skipped
   * 
   */
  class SeedTestPolicy: public art::EDAnalyzer {
    
      public:
    
    explicit SeedTestPolicy(fhicl::ParameterSet const& pset);
    
    virtual void analyze(art::Event const& event) override;
    
    virtual void endJob() override;
    
      private:
    using seed_t = testing::LArSeedService::seed_t;
    
    std::vector<std::string> instanceNames;
    std::vector<seed_t> startSeeds; ///< seeds after the constructor
    unsigned int nExpectedErrors; ///< number of expected errors
    bool useGenerators; ///< instanciate and use random number generators
    bool perEventSeeds; ///< whether we expect different seeds on each event
    
    unsigned int nErrors = 0; ///< Number of errors detected so far
    
    std::unique_ptr<CLHEP::HepRandomEngine> localEngine; ///< self-managed 
    
    /// Returns whether the engine associated with the specified index is local
    bool isLocalEngine(size_t iEngine) const;
    
    /// Returns the engine associated with the given instance index
    /// (may throw if none is present)
    CLHEP::HepRandomEngine& getRandomEngine(size_t iEngine);
  
    /// Verifies that the engine with the specified index has the expected seed
    seed_t verifySeed(size_t iEngine);
    
    seed_t obtainSeed(std::string instanceName = "");
    
    /// Returns whether e is an exception we can handle (and, if so, it handles)
    bool handleSeedServiceException(art::Exception& e);
    
  }; // class SeedTestPolicy
  
} // namespace testing



//------------------------------------------------------------------------------
//--- implementation
//--- 

testing::SeedTestPolicy::SeedTestPolicy(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)
  , instanceNames  (pset.get<std::vector<std::string>>("instanceNames", {}))
  , nExpectedErrors(pset.get<unsigned int>            ("expectedErrors", 0U))
  , useGenerators  (pset.get<bool>                    ("useGenerators", true))
  , perEventSeeds  (pset.get<bool>                    ("perEventSeeds", false))
{
  
  //
  // print some configuration information
  //
  { // anonymous block
    mf::LogInfo log("SeedTestPolicy");
    log << "Construct SeedTestPolicy with "
      << instanceNames.size() <<  " engine instances:";
    for (auto const& instanceName: instanceNames)
      log << " " << instanceName;
    
  } // anonymous block
  
  auto* Seeds = &*(art::ServiceHandle<sim::LArSeedService>());
  
  // by default, have at least one, default engine instance
  if (instanceNames.empty()) instanceNames.push_back("");
  
  mf::LogInfo log("SeedTestPolicy"); // cumulative log
  
  //
  // register all the engines, and store their seeds
  //
  for (std::string const& instanceName: instanceNames) {
    seed_t const seed = obtainSeed(instanceName);
    log << "\nSeed for '" << instanceName << "' is: " << seed;
    startSeeds.push_back(seed);
  } // for first loop (declaration)
  
  
  //
  // verify the seed of each instance
  //
  for (size_t iEngine = 0; iEngine < instanceNames.size(); ++iEngine) {
    std::string const& instanceName = instanceNames[iEngine];
    
    // This involved condition tree ensures that SeedMaster is queried
    // for a seed exactly once per instance, no matter what.
    // This is relevant for the error count.
    // Out of it, a seed is returned.
    seed_t seed = sim::LArSeedService::InvalidSeed;
    if (isLocalEngine(iEngine)) {
      if (useGenerators) {
        localEngine = std::make_unique<CLHEP::HepJamesRandom>();
        try {
          seed = Seeds->defineEngine(*localEngine, instanceName);
        }
        catch(art::Exception& e) {
          if (!handleSeedServiceException(e)) throw;
        }
        mf::LogInfo("SeedTestConstruct")
          << "Engine instance '" << instanceName
          << "' will be owned by the test module.";
      } // if use generators
      else seed = obtainSeed(instanceName);
    } // if local
    else { // if managed by art
      seed = obtainSeed(instanceName);
      if (useGenerators) {
        createEngine(seed, "HepJamesRandom", instanceName);
        // registration still matters for per-event policies
        Seeds->defineEngine(getRandomEngine(iEngine), instanceName);
        verifySeed(iEngine);
      }
    } // if ... else
    
    // check that the seed returned by the service is still the same
    seed_t const expectedSeed = startSeeds.at(iEngine);
    if (seed != expectedSeed) {
      throw art::Exception(art::errors::LogicError)
        << "LArSeedService returned different seed values for engine instance '"
        << instanceName << "': first " << expectedSeed << ", now " << seed
        << "\n";
    } // if unexpected seed
  } // for second loop

  //
  // An engine with the following label has already been registered
  // (it's the one managed by RandomGeneratorService).
  // Registering another should raise an exception
  // (incidentally, if useGenerators is false we are trying to register nullptr)
  // 
  bool bBug = false;
  try {
    Seeds->declareEngine(instanceNames.front());
    bBug = true;
  }
  catch(std::exception const& e) {
    if (!testing::LArSeedService::isSeedServiceException(e)) throw;
  }
  if (bBug) {
    throw art::Exception(art::errors::LogicError)
      << "Registration of local engine with duplicate label"
         " did not throw an exception";
  }
  
} // testing::SeedTestPolicy::SeedTestPolicy()


//------------------------------------------------------------------------------
void testing::SeedTestPolicy::analyze(art::Event const& event) {

  mf::LogVerbatim("SeedTestPolicy")
    << "SeedTestPolicy::analyze() " << event.id() << " with "
    << instanceNames.size() << " random engines";
  
  if (useGenerators) {
    for (size_t iEngine = 0; iEngine < instanceNames.size(); ++iEngine) {
      //
      // collect information and resources
      //
      std::string const& instanceName = instanceNames[iEngine];
      seed_t const startSeed = startSeeds.at(iEngine);
      CLHEP::HepRandomEngine& engine = getRandomEngine(iEngine);
      
      //
      // check seed (if per event, it should be the opposite)
      //
      seed_t const actualSeed = testing::LArSeedService::readSeed(engine);
      if (perEventSeeds) {
        if (actualSeed == startSeed) {
          // this has a ridiculously low chance of begin fortuitous
          throw art::Exception(art::errors::LogicError)
            << "per event seed " << actualSeed << " of engine '" << instanceName
            << "' is the same as at beginning!\n";
        }
      }
      else {
        if (actualSeed != startSeed) {
          throw art::Exception(art::errors::LogicError)
            << "expected seed " << startSeed << " for engine '" << instanceName
            << "', got " << actualSeed << " instead!\n";
        }
      }
      
      //
      // print character statistics
      //
      mf::LogVerbatim("SeedTestPolicy")
        << std::setw(12) << (instanceName.empty()? "<default>": instanceName)
        << ": " << testing::LArSeedService::CreateCharacter(engine)
        << "   (seed: " << actualSeed << ")";
      
    } // for
  } // if use generators
  
} // testing::SeedTestPolicy::analyze()


//------------------------------------------------------------------------------
void testing::SeedTestPolicy::endJob() {
  // if we have an unexpected amount of errors, bail out
  if (nExpectedErrors != nErrors) {
    art::Exception e(art::errors::Configuration);
    e << "SeedTestPolicy: detected " << nErrors << " errors";
    if (nExpectedErrors) e << ", " << nExpectedErrors << " expected";
    throw e << "!\n";
  }
} // testing::SeedTestPolicy::endJob()


//------------------------------------------------------------------------------

testing::SeedTestPolicy::seed_t testing::SeedTestPolicy::obtainSeed
  (std::string instanceName /* = "" */)
{
  // Returns the seed for the specified engine instance, or 0 in case of
  // configuration error (in which case, an error counter is increased)
  seed_t seed = sim::LArSeedService::InvalidSeed;
  try {
    art::ServiceHandle<sim::LArSeedService> seeds;
    // currently (v0_00_03), the two calls are actually equivalent
    seed
      = instanceName.empty()? seeds->getSeed(): seeds->getSeed(instanceName);
  }
  catch(art::Exception& e) {
    if (!testing::LArSeedService::isSeedServiceException(e)) throw;
    
    ++nErrors;
    mf::LogError log("SeedTestPolicy");
    log << "Detected";
    if (nErrors > nExpectedErrors) log << " UNEXPECTED";
    log << " error #" << nErrors << ":\n" << e;
  }
  return seed;
} // testing::SeedTestPolicy::obtainSeed()


bool testing::SeedTestPolicy::handleSeedServiceException(art::Exception& e) {
  if (!testing::LArSeedService::isSeedServiceException(e)) return false;
  
  ++nErrors;
  mf::LogError log("SeedTest01");
  log << "Detected";
  if (nErrors > nExpectedErrors) log << " UNEXPECTED";
  log << " error #" << nErrors << ":\n" << e << "\n";
  return true;
} // testing::SeedTestPolicy::handleSeedServiceException()


bool testing::SeedTestPolicy::isLocalEngine(size_t i) const {
  return (i == 0) && (instanceNames.size() != 1);
} // testing::SeedTestPolicy::isLocalEngine()


CLHEP::HepRandomEngine& testing::SeedTestPolicy::getRandomEngine
  (size_t iEngine)
{
  // do we have a local engine?
  // it's the first one, only if there are more than one engine
  return isLocalEngine(iEngine)
    ? (*localEngine)
    : art::ServiceHandle<art::RandomNumberGenerator>()
        ->getEngine(instanceNames.at(iEngine))
    ;
} // testing::SeedTestPolicy::getRandomEngine()


testing::SeedTestPolicy::seed_t testing::SeedTestPolicy::verifySeed
  (size_t iEngine)
{
  CLHEP::HepRandomEngine const& engine = getRandomEngine(iEngine);
  
  seed_t const actualSeed = testing::LArSeedService::readSeed(engine);
  seed_t const expectedSeed = startSeeds.at(iEngine);
  // if the expected seed is invalid, we are not even sure it was ever set;
  // the engine is in an invalid state and that's it
  if (!sim::LArSeedService::isSeedValid(expectedSeed)) return actualSeed;
  
  if (actualSeed != expectedSeed) {
    std::string const& instanceName = instanceNames[iEngine];
    throw art::Exception(art::errors::LogicError)
      << "expected seed " << expectedSeed << " for engine '" << instanceName
      << "', got " << actualSeed << " instead!";
  }
  return actualSeed;
} // testing::SeedTestPolicy::verifySeed()  
  
  
//------------------------------------------------------------------------------
DEFINE_ART_MODULE(testing::SeedTestPolicy)

//------------------------------------------------------------------------------
