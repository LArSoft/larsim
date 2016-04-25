/**
 * @file SeedTest01_module.cc
 * @brief Test the LArSeedService.
 * @author Rob Kutschke (kutschke@fnal.gov)
 * @see LArSeedService.hh
 */


// art extensions
#include "larsim/RandomUtils/LArSeedService.h"

// C++ includes.
#include <string>
#include <vector>
#include <sstream>
#include <iomanip> // std::setw()

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine
#include "CLHEP/Random/RandFlat.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// Framework includes.
#include "art/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/CurrentModule.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace testing {
  
  /// Returns whether the exception looks to be from LArSeedService
  bool isSeedServiceException(std::exception const& e) {
    
    // is it an art::Exception? (LArSeedService exceptions are)
    art::Exception const* art_exc = dynamic_cast<art::Exception const*>(&e);
    if (!art_exc) return false;
    
    // All configuration exceptions are from LArSeedService.
    // Why? just because.
    if (art_exc->categoryCode() == art::errors::Configuration)
      return true;
    
    // does it have "LArSeedService" in the message? (LArSeedService exceptions should)
    if ( (art_exc->explain_self().find("LArSeedService") == std::string::npos)
      && (art_exc->explain_self().find("SeedMaster") == std::string::npos)
      )
      return false;
    
    // we can't do better than this
    return true;
  } // isSeedServiceException(exception)
  
  
  /// Returns whether the exception looks to be from LArSeedService and has specified code
  bool isSeedServiceException
    (std::exception const& e, art::errors::ErrorCodes code)
  {
    if (!isSeedServiceException(e)) return false;
    // by now we know "e" is a art::Exception
    return static_cast<art::Exception const&>(e).categoryCode() == code;
  } // isSeedServiceException(exception, ErrorCodes)
  
  
  /**
   * @brief Test module for LArSeedService
   * 
   * The test writes on screen the random seeds it gets.
   * 
   * Note that the test does not actually get any random number, unless the
   * *useGenerators* option is turned on.
   * 
   * Configuration parameters:
   * - <b>testMode</b> (integer, required):
   *   - 0: tries to get a seed on construction (twice)
   *   - 1: tries to get a seed for each instance on construction (twice)
   *   - 2: tries to get a seed on each run
   *   - 3: tries to get a seed on each subrun (should fail)
   *   - 4: tries to get a seed on each event (should fail)
   * - <b>instanceNames</b> (string list, required): use one random number
   *   generator for each instance name here specified
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
  class SeedTest01 : public art::EDAnalyzer {

  public:

    typedef art::RandomNumberGenerator::seed_t seed_t;

    explicit SeedTest01(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& event) override;

    void beginRun   ( const art::Run&    run) override;
    void beginSubRun( const art::SubRun& subRun) override;
    
    void endJob() override;
    
  private:
    int testMode_;
    std::vector<std::string> instanceNames_;
    std::vector<seed_t> startSeeds_; ///< seeds after the constructor
    unsigned int nExpectedErrors; ///< number of expected errors
    bool useGenerators; ///< instanciate and use random number generators
    bool bPerEventSeeds; ///< whether we expect different seeds on each event
    
    unsigned int nErrors = 0; ///< Number of errors detected so far
    
    std::string moduleName() const;
    
    seed_t ObtainSeed(std::string instanceName = "");
    
    
    static unsigned short RollStat(CLHEP::RandFlat& rand);
    static std::string CreateCharacter(std::string instance);
    
    static seed_t readSeed(CLHEP::HepRandomEngine const& engine)
      { return engine.getSeed(); }
    
  }; // class SeedTest01

  
  SeedTest01::SeedTest01(fhicl::ParameterSet const& pSet):
    art::EDAnalyzer(pSet),
    testMode_(pSet.get<int>("testMode")),
    instanceNames_(),
    startSeeds_(),
    nExpectedErrors(pSet.get<unsigned int>("expectedErrors", 0U)),
    useGenerators(pSet.get<bool>("useGenerators", true)),
    bPerEventSeeds(pSet.get<bool>("perEventSeeds", false))
  {
    pSet.get_if_present("instanceNames", instanceNames_ );

    { // anonymous block
      mf::LogInfo log("SeedTestConstruct");
      log << "Construct SeedTest01[" << moduleName() << "]. Mode:  "
        << testMode_ <<  " instances:";
      for ( size_t i=0; i< instanceNames_.size(); ++i){
        log << " " << instanceNames_[i];
      }
    }
    
    switch (testMode_) {
      case 0: {
        if (!instanceNames_.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "SeedTest01: Test mode 0 does not use instances.";
        }
        instanceNames_.push_back("");
        
        seed_t seed = ObtainSeed();
        seed_t seed2 = ObtainSeed();
        mf::LogVerbatim("SeedTestConstruct")
          <<   "Seed is :  " << seed
          << "\nSeed2 is : " << seed2;
        
        startSeeds_.resize(instanceNames_.size());
        startSeeds_[0] = seed;
        
        if (seed != seed2) {
          throw art::Exception(art::errors::LogicError)
            << "SeedTest01: Seed has changed: " << seed << " => " << seed2;
        }
        if (useGenerators) createEngine(seed);
        break;
      } // test mode 0
      
      case 1: {
        if (instanceNames_.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "SeedTest01: Instance names are required!";
        }
        
        mf::LogVerbatim log("SeedTestConstruct");
        log << "Module '" << moduleName() << "' has " << instanceNames_.size()
          << " seed instances:";
        for ( size_t i=0; i< instanceNames_.size(); ++i){
          seed_t seed = ObtainSeed( instanceNames_[i] );
          log << "\nSeed for '"
              << instanceNames_[i] << "' is: "
              << seed;
          startSeeds_.push_back(seed);
        } // for first loop
        for ( size_t i=0; i< instanceNames_.size(); ++i){
          seed_t seed = ObtainSeed( instanceNames_[i] );
          log << "\nSeed2 for '"
              << instanceNames_[i] << "' is: "
              << seed;
          if (seed != startSeeds_[i]) {
            throw art::Exception(art::errors::LogicError)
              << "SeedTest01: Seed has changed for instance "
              << instanceNames_[i] << ": " << startSeeds_[i] << " => " << seed;
          }
          if (useGenerators)
            createEngine(seed, "HepJamesRandom", instanceNames_[i]);
        } // for second loop
        break;
      } // if test mode 1
      
      case 2: case 3: case 4: {
        if (!instanceNames_.empty()) {
          throw art::Exception(art::errors::Configuration)
            << "SeedTest01: Test mode " << testMode_ << " does not use instances.";
        }
        instanceNames_.push_back("");
        if (useGenerators) {
          // the point of these test modes is not to make known to LArSeedService
          // that we have random engine instances; therefore we initialize
          // with a default value
          seed_t seed = 0;
          createEngine(seed);
          startSeeds_.push_back(seed);
        }
        break;
      } // if test mode 2, 3 or 4
      
    } // switch (testMode_)
  } // SeedTest01::SeedTest01()

  void SeedTest01::analyze(const art::Event& event){
    { // anonymous block
      mf::LogVerbatim log("SeedTest01");
      log << "SeedTest01[" << moduleName() << "]::analyze "
          << event.id();
      if ( testMode_ == 4 ){
        seed_t seed = ObtainSeed();
        log << "\nSeed from analyze is :  " << seed;
      }
    } // anonymous block
    
    if (useGenerators) {
      for(size_t iInst = 0; iInst < instanceNames_.size(); ++iInst) {
        std::string const& instanceName = instanceNames_[iInst];
        seed_t startSeed = startSeeds_.at(iInst);
        art::ServiceHandle<art::RandomNumberGenerator> rng;
        CLHEP::HepRandomEngine& engine = rng->getEngine(instanceName);
        seed_t actualSeed = readSeed(engine);
        if ((actualSeed != startSeed) && !bPerEventSeeds) {
          throw art::Exception(art::errors::LogicError)
            << "expected seed " << startSeed << " for engine '" << instanceName
            << "', got " << actualSeed << " instead!";
        }
        mf::LogVerbatim("SeedTest01")
          << std::setw(12) << (instanceName.empty()? "<default>": instanceName)
          << ": " << CreateCharacter(instanceName)
          << "   (seed: " << actualSeed << ")";
      } // for
    } // if use generators
  } // SeedTest01::analyze()

  void SeedTest01::beginRun( const art::Run& run){
    mf::LogVerbatim log("SeedTest01");
    log << "SeedTest01[" << moduleName() << "]::beginRun " << run.id();
    if ( testMode_ == 2 ){
      seed_t seed = ObtainSeed();
      log << "\nSeed from beginRun is :  " << seed;
    }
  }

  void SeedTest01::beginSubRun( const art::SubRun& subRun){
    mf::LogVerbatim log("SeedTest01");
    log << "SeedTest01[" << moduleName() << "]::beginSubRun " << subRun.id();
    if ( testMode_ == 3 ){
      seed_t seed = ObtainSeed();
      log << "Seed from beginSubRun is :  " << seed;
    }
  }
  
  void SeedTest01::endJob() {
    // if we have an unexpected amount of errors, bail out
    if (nExpectedErrors != nErrors) {
      art::Exception e(art::errors::Configuration);
      e << "SeedTest01: detected " << nErrors << " errors";
      if (nExpectedErrors) e << ", " << nExpectedErrors << " expected";
      throw e << "!";
    }
  } // SeedTest01::endJob()
  
  
  inline std::string SeedTest01::moduleName() const
    { return art::ServiceHandle<art::CurrentModule>()->label(); }
  
  
  SeedTest01::seed_t SeedTest01::ObtainSeed(std::string instanceName /* = "" */)
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
      if (!isSeedServiceException(e)) throw;
      
      ++nErrors;
      mf::LogError log("SeedTest01");
      log << "Detected";
      if (nErrors > nExpectedErrors) log << " UNEXPECTED";
      log << " error #" << nErrors << ":\n" << e;
    }
    return seed;
  } // SeedTest01::ObtainSeed()
  
  
  unsigned short SeedTest01::RollStat(CLHEP::RandFlat& rand) {
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
  } // SeedTest01::RollStat()
  
  
  std::string SeedTest01::CreateCharacter(std::string instance) {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(instance);
    CLHEP::RandFlat flat(engine);
    
    std::ostringstream sstr;
    sstr
      << "  STR=" << std::setw(2) << RollStat(flat)
      << "  DEX=" << std::setw(2) << RollStat(flat)
      << "  CON=" << std::setw(2) << RollStat(flat)
      << "  INT=" << std::setw(2) << RollStat(flat)
      << "  WIS=" << std::setw(2) << RollStat(flat)
      << "  CHA=" << std::setw(2) << RollStat(flat)
      ;
    return sstr.str();
  } // SeedTest01::CreateCharacter()

} // namespace testing

DEFINE_ART_MODULE(testing::SeedTest01)
