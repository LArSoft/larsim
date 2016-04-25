/**
 * @file SeedService_service.cc
 * @brief Assists in the distribution of guaranteed unique seeds to all engines within a job.
 * @author Rob Kutschke (kutschke@fnal.gov)
 * @see LArSeedService.h SeedMaster.h
 */

// LArSeedService header
#include "larsim/RandomUtils/LArSeedService.h"

// Art include files
#include "art/Utilities/Exception.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Persistency/Provenance/EventID.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine

// C++ include files
#include <iostream>
#include <iomanip>

namespace sim {

  //----------------------------------------------------------------------------
  LArSeedService::LArSeedService
    (fhicl::ParameterSet const& paramSet, art::ActivityRegistry& iRegistry)
    : seeds(paramSet)
    , state()
    , verbosity(paramSet.get<int>("verbosity", 0))
    , bPrintEndOfJobSummary(paramSet.get<bool>("endOfJobSummary",false))
  {
    // Register callbacks.
    iRegistry.sPreModuleConstruction.watch  (this, &LArSeedService::preModuleConstruction  );
    iRegistry.sPostModuleConstruction.watch (this, &LArSeedService::postModuleConstruction );
    iRegistry.sPreModuleBeginRun.watch      (this, &LArSeedService::preModuleBeginRun      );
    iRegistry.sPostModuleBeginRun.watch     (this, &LArSeedService::postModuleBeginRun     );
    iRegistry.sPreProcessEvent.watch        (this, &LArSeedService::preProcessEvent        );
    iRegistry.sPreModule.watch              (this, &LArSeedService::preModule              );
    iRegistry.sPostModule.watch             (this, &LArSeedService::postModule             );
    iRegistry.sPostProcessEvent.watch       (this, &LArSeedService::postProcessEvent       );
    iRegistry.sPostEndJob.watch             (this, &LArSeedService::postEndJob             );
    
  } // LArSeedService::LArSeedService()
  
  
  
  //----------------------------------------------------------------------------
  LArSeedService::EngineId LArSeedService::qualify_engine_label
    (std::string moduleLabel, std::string instanceName) const
    { return { moduleLabel, instanceName }; }
  
  LArSeedService::EngineId LArSeedService::qualify_engine_label
    (std::string instanceName /* = "" */) const
    { return qualify_engine_label( state.moduleLabel(), instanceName); }
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::getSeed() {
    return getSeed(qualify_engine_label());
  } // LArSeedService::getSeed()
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::getSeed(std::string instanceName) {
    return getSeed(qualify_engine_label(instanceName));
  } // LArSeedService::getSeed(string)
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::print() const { seeds.print(); }
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::getSeed(EngineId const& id) {
    
    // Are we being called from the right place?
    ensureValidState();
    
    // Ask the seed master about seeds...
    return seeds.getSeed(id);
  } // LArSeedService::getSeed(EngineID)
  
  
  std::pair<LArSeedService::seed_t, bool> LArSeedService::getSeed(
    EngineId const& id,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
  ) {
    seed_t seed = InvalidSeed;
    bool bFrozen = readSeedParameter(seed, pset, pname);
    // if we got a seed, but it is invalid, we pretend we got nothing
    if (bFrozen && (seed == InvalidSeed)) bFrozen = false;
    return { (bFrozen? seed: getSeed(id)), bFrozen };
  } // LArSeedService::getSeed(ParameterSet)
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::createEngine(
    art::EngineCreator& module, std::string type,
    std::string instance /* = "" */
  ) {
    EngineId id = qualify_engine_label(instance);
    const seed_t seed = prepareEngine(id, RandomNumberGeneratorSeeder);
    module.createEngine(seed, type, instance);
    mf::LogInfo("LArSeedService")
      << "Seeding " << type << " engine \"" << id.artName()
      << "\" with seed " << seed << ".";
    return seed;
  } // LArSeedService::createEngine(strings)
  
  
  LArSeedService::seed_t LArSeedService::createEngine(art::EngineCreator& module) {
    EngineId id = qualify_engine_label();
    const seed_t seed = prepareEngine(id, RandomNumberGeneratorSeeder);
    module.createEngine(seed);
    mf::LogInfo("LArSeedService")
      << "Seeding default-type engine \"" << id.artName()
      << "\" with seed " << seed << ".";
    return seed;
  } // LArSeedService::createEngine()
  
  
  LArSeedService::seed_t LArSeedService::createEngine(
    art::EngineCreator& module, std::string type,
    std::string instance,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
  ) {
    EngineId id = qualify_engine_label(instance);
    registerEngineAndSeeder(id, RandomNumberGeneratorSeeder);
    std::pair<seed_t, bool> seedInfo = getSeed(id, pset, pname);
    module.createEngine(seedInfo.first, type, instance);
    mf::LogInfo("LArSeedService")
      << "Seeding " << type << " engine \"" << id.artName()
      << "\" with seed " << seedInfo.first << ".";
    if (seedInfo.second) freezeSeed(id);
    return seedInfo.first;
  } // LArSeedService::createEngine(ParameterSet)
  
  
  LArSeedService::seed_t LArSeedService::createEngine(
    art::EngineCreator& module,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
  ) {
    EngineId id = qualify_engine_label();
    registerEngineAndSeeder(id, RandomNumberGeneratorSeeder);
    std::pair<seed_t, bool> seedInfo = getSeed(id, pset, pname);
    module.createEngine(seedInfo.first);
    mf::LogInfo("LArSeedService")
      << "Seeding default-type engine \"" << id.artName()
      << "\" with seed " << seedInfo.first << ".";
    if (seedInfo.second) freezeSeed(id);
    return seedInfo.first;
  } // LArSeedService::createEngine(ParameterSet)
  
  
  LArSeedService::seed_t LArSeedService::registerEngine
    (SeedMaster_t::Seeder_t seeder, std::string instance /* = "" */)
  {
    EngineId id = qualify_engine_label(instance);
    const seed_t seed = prepareEngine(id, seeder);
    if (seeder) {
      seeder(id, seed);
      mf::LogInfo("LArSeedService")
        << "Seeding registered engine \"" << id.artName()
        << "\" with seed " << seed << ".";
    }
    return seed;
  } // LArSeedService::registerEngine()
  
  
  LArSeedService::seed_t LArSeedService::registerEngine(
    SeedMaster_t::Seeder_t seeder, std::string instance,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pname
  ) {
    EngineId id = qualify_engine_label(instance);
    registerEngineAndSeeder(id, seeder);
    std::pair<seed_t, bool> seedInfo = getSeed(id, pset, pname);
    if (seeder) {
      seeder(id, seedInfo.first);
      mf::LogInfo("LArSeedService")
        << "Seeding registered engine \"" << id.artName()
        << "\" with seed " << seedInfo.first << ".";
    }
    if (seedInfo.second) freezeSeed(id);
    return seedInfo.first;
  } // LArSeedService::registerEngine()
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::ensureValidState() const {
    // getSeed may only be called from a c'tor or from a beginRun method.
    // In all other cases, throw.
    if ( (state.state() != SeedServiceHelper::ArtState::inModuleConstructor)
      && (state.state() != SeedServiceHelper::ArtState::inBeginRun)
      )
    {
      throw art::Exception(art::errors::LogicError)
        << "LArSeedService: not in a module constructor or beginRun method. May not call getSeed.";
    }
  } // LArSeedService::ensureValidState()
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::reseedInstance(EngineId const& id) {
    // get all the information on the current process, event and module from
    // ArtState:
    SeedMaster_t::EventData_t data(state.getEventSeedInputData());
    seed_t seed = seeds.reseedEvent(id, data);
    if (seed == InvalidSeed) {
      mf::LogDebug("LArSeedService")
        << "No random seed specific to this event for engine '" << id << "'";
    }
    else {
      mf::LogInfo("LArSeedService") << "Random seed for this event, engine '"
        << id << "': " << seed;
    }
    return seed;
  } // LArSeedService::reseedInstance()
  
  
  void LArSeedService::reseedModule(std::string currentModule) {
    for (EngineId const& ID: seeds.engineIDsRange()) {
      if (ID.moduleLabel != currentModule) continue; // not our module? neeext!!
      reseedInstance(ID);
    } // for
  } // LArSeedService::reseedModule(string)
  
  void LArSeedService::reseedModule() { reseedModule(state.moduleLabel()); }
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::RandomNumberGeneratorSeeder
    (EngineId const& id, seed_t seed)
  {
    // no check is performed to verify that the current module is the one
    // specified in id.moduleLabel -- but that is required!
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(id.instanceName);
    engine.setSeed(seed, 0); // the 0 is dummy... or so one hopes
  } // LArSeedService::RandomNumberGeneratorSeeder()
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::registerEngineAndSeeder
    (EngineId const& id, SeedMaster_t::Seeder_t seeder)
  {
    if (hasEngine(id)) {
      throw art::Exception(art::errors::LogicError)
        << "LArSeedService: an engine with ID '" << id
        << "' has already been created!";
    }
    seeds.registerSeeder(id, seeder);
  } // LArSeedService::registerEngine()
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::freezeSeed(EngineId const& id)
    { seeds.freezeSeed(id); }
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::prepareEngine
    (EngineId const& id, SeedMaster_t::Seeder_t seeder)
  {
    registerEngineAndSeeder(id, seeder);
    return getSeed(id);
  } // LArSeedService::prepareEngine()
  
  
  //----------------------------------------------------------------------------
  bool LArSeedService::readSeedParameter(
    seed_t& seed, fhicl::ParameterSet const& pset,
    std::initializer_list<std::string> pnames
  ) {
    for (std::string const& key: pnames)
      if (pset.get_if_present(key, seed)) return true;
    seed = InvalidSeed;
    return false;
  } // LArSeedService::readSeedParameter()
  
  
  //----------------------------------------------------------------------------
  // Callbacks called by art.  Used to maintain information about state.
  void LArSeedService::preModuleConstruction(art::ModuleDescription const& md)  {
    state.set_state(SeedServiceHelper::ArtState::inModuleConstructor);
    state.set_module(md);
  } // LArSeedService::preModuleConstruction()
  
  void LArSeedService::postModuleConstruction(art::ModuleDescription const&) {
    state.reset_state();
  } // LArSeedService::postModuleConstruction()
  
  void LArSeedService::preModuleBeginRun(art::ModuleDescription const& md) {
    state.set_state(SeedServiceHelper::ArtState::inBeginRun);
    state.set_module(md);
  } // LArSeedService::preModuleBeginRun()
  
  void LArSeedService::postModuleBeginRun(art::ModuleDescription const&) {
    state.reset_state();
  } // LArSeedService::postModuleBeginRun()
  
  void LArSeedService::preProcessEvent(art::Event const& evt) {
    state.set_state(SeedServiceHelper::ArtState::inEvent);
    state.set_event(evt);
    seeds.onNewEvent(); // inform the seed master that a new event has come
  } // LArSeedService::preProcessEvent()
  
  void LArSeedService::preModule(art::ModuleDescription const& md) {
    state.set_state(SeedServiceHelper::ArtState::inModuleEvent);
    state.set_module(md);
    
    // Reseed all the engine of this module... maybe
    // (that is, if the current policy alows it).
    reseedModule(md.moduleLabel());
  } // LArSeedService::preModule()
  
  void LArSeedService::postModule(art::ModuleDescription const&) {
    state.reset_module();
    state.reset_state();
  } // LArSeedService::postModule()
  
  void LArSeedService::postProcessEvent(art::Event const&) {
    state.reset_event();
    state.reset_state();
  } // LArSeedService::postProcessEvent()
  
  
  void LArSeedService::postEndJob() {
    if ((verbosity > 0) || bPrintEndOfJobSummary)
      print(); // framework logger decides whether and where it shows up
  } // LArSeedService::postEndJob()
  
  //----------------------------------------------------------------------------

} // end namespace sim

DEFINE_ART_SERVICE(sim::LArSeedService)
