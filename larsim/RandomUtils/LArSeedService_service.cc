/**
 * @file SeedService_service.cc
 * @brief Assists in the distribution of guaranteed unique seeds to all engines within a job.
 * @author Rob Kutschke (kutschke@fnal.gov)
 * @see LArSeedService.h SeedMaster.h
 */

// LArSeedService header
#include "larsim/RandomUtils/LArSeedService.h"

// Art include files
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EngineCreator.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Provenance/ModuleDescription.h"
#include "canvas/Persistency/Provenance/EventID.h"

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
    state.transit_to(SeedServiceHelper::ArtState::inServiceConstructor);
    
    // Register callbacks.
    iRegistry.sPreModuleConstruction.watch  (this, &LArSeedService::preModuleConstruction  );
    iRegistry.sPostModuleConstruction.watch (this, &LArSeedService::postModuleConstruction );
    iRegistry.sPreModuleBeginRun.watch      (this, &LArSeedService::preModuleBeginRun      );
    iRegistry.sPostModuleBeginRun.watch     (this, &LArSeedService::postModuleBeginRun     );
    iRegistry.sPreProcessEvent.watch        (this, &LArSeedService::preProcessEvent        );
    iRegistry.sPreModule.watch              (this, &LArSeedService::preModule              );
    iRegistry.sPostModule.watch             (this, &LArSeedService::postModule             );
    iRegistry.sPostProcessEvent.watch       (this, &LArSeedService::postProcessEvent       );
    iRegistry.sPreModuleEndJob.watch        (this, &LArSeedService::preModuleEndJob        );
    iRegistry.sPostModuleEndJob.watch       (this, &LArSeedService::postModuleEndJob       );
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
  LArSeedService::seed_t LArSeedService::getGlobalSeed(std::string instanceName) {
    EngineId ID(instanceName, EngineId::global);
    LOG_DEBUG("LArSeedService")
      << "LArSeedService::getGlobalSeed(\"" << instanceName << "\")";
    return getSeed(ID);
  } // LArSeedService::getGlobalSeed()
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::getSeed(EngineId const& id) {
    
    // We require an engine to have been registered before we yield seeds;
    // this should minimise unexpected conflicts.
    if (hasEngine(id)) return querySeed(id); // ask the seed to seed master
    
    // if it hasn't been declared, we declare it now
    // (this is for backward compatibility with the previous behaviour).
    // registerEngineID() will eventually call this function again to get the
    // seed... so we return it directly.
    // Also note that this effectively "freezes" the engine since no seeder
    // is specified.
    return registerEngineID(id);
        
  } // LArSeedService::getSeed(EngineId)
  
  
  LArSeedService::seed_t LArSeedService::querySeed(EngineId const& id) {
    return seeds.getSeed(id); // ask the seed to seed master
  } // LArSeedService::querySeed()
  
  
  std::pair<LArSeedService::seed_t, bool> LArSeedService::findSeed(
    EngineId const& id,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
  ) {
    seed_t seed = InvalidSeed;
    // try and read the seed from configuration; if succeed, it's "frozen"
    bool const bFrozen = readSeedParameter(seed, pset, pnames);
    
    // if we got a valid seed, use it as frozen
    if (bFrozen && (seed != InvalidSeed))
      return { seed, true };  
    
    // seed was not good enough; get the seed from the master
    return { querySeed(id), false };
    
  } // LArSeedService::findSeed()
  
  
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
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
  ) {
    EngineId id = qualify_engine_label(instance);
    registerEngineAndSeeder(id, RandomNumberGeneratorSeeder);
    std::pair<seed_t, bool> seedInfo = findSeed(id, pset, pnames);
    module.createEngine(seedInfo.first, type, instance);
    mf::LogInfo("LArSeedService")
      << "Seeding " << type << " engine \"" << id.artName()
      << "\" with seed " << seedInfo.first << ".";
    if (seedInfo.second) freezeSeed(id, seedInfo.first);
    return seedInfo.first;
  } // LArSeedService::createEngine(ParameterSet)
  
  
  LArSeedService::seed_t LArSeedService::createEngine(
    art::EngineCreator& module,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
  ) {
    EngineId id = qualify_engine_label();
    registerEngineAndSeeder(id, RandomNumberGeneratorSeeder);
    std::pair<seed_t, bool> seedInfo = findSeed(id, pset, pnames);
    module.createEngine(seedInfo.first);
    mf::LogInfo("LArSeedService")
      << "Seeding default-type engine \"" << id.artName()
      << "\" with seed " << seedInfo.first << ".";
    if (seedInfo.second) freezeSeed(id, seedInfo.first);
    return seedInfo.first;
  } // LArSeedService::createEngine(ParameterSet)
  
  
  LArSeedService::seed_t LArSeedService::registerEngine
    (SeedMaster_t::Seeder_t seeder, std::string instance /* = "" */)
  {
    return registerEngineID(qualify_engine_label(instance), seeder);
  } // LArSeedService::registerEngine(Seeder_t, string)
  
  
  LArSeedService::seed_t LArSeedService::registerEngine(
    SeedMaster_t::Seeder_t seeder, std::string instance,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
  ) {
    EngineId id = qualify_engine_label(instance);
    registerEngineAndSeeder(id, seeder);
    std::pair<seed_t, bool> seedInfo = findSeed(id, pset, pnames);
    seedEngine(id); // seed it before freezing
    if (seedInfo.second) freezeSeed(id, seedInfo.first);
    seed_t const seed = seedInfo.first;
    return seed;
  } // LArSeedService::registerEngine(Seeder_t, string, ParameterSet, init list)
  
    
  LArSeedService::seed_t LArSeedService::declareEngine(std::string instance) {
    return registerEngine(SeedMaster_t::Seeder_t(), instance);
  } // LArSeedService::declareEngine(string)
  
  
  LArSeedService::seed_t LArSeedService::declareEngine(
    std::string instance,
    fhicl::ParameterSet const& pset, std::initializer_list<std::string> pnames
  ) {
    return registerEngine(SeedMaster_t::Seeder_t(), instance, pset, pnames);
  } // LArSeedService::declareEngine(string, ParameterSet, init list)
  
  
  LArSeedService::seed_t LArSeedService::defineEngine
    (SeedMaster_t::Seeder_t seeder, std::string instance /* = {} */)
  {
    return defineEngineID(qualify_engine_label(instance), seeder);
  } // LArSeedService::defineEngine(string, Seeder_t)
  
    
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::registerEngineID(
    EngineId const& id,
    SeedMaster_t::Seeder_t seeder /* = SeedMaster_t::Seeder_t() */
  ) {
    prepareEngine(id, seeder);
    return seedEngine(id);
  } // LArSeedService::registerEngine()
  
  
  LArSeedService::seed_t LArSeedService::defineEngineID
    (EngineId const& id, SeedMaster_t::Seeder_t seeder)
  {
    if (!hasEngine(id)) {
      throw art::Exception(art::errors::LogicError)
        << "Attempted to define engine '" << id.artName()
        << "', that was not declared\n";
    }
    
    if (seeds.hasSeeder(id)) {
      throw art::Exception(art::errors::LogicError)
        << "Attempted to redefine engine '" << id.artName()
        << "', that has already been defined\n";
    }
    
    ensureValidState();
    
    seeds.registerSeeder(id, seeder);
    seed_t const seed = seedEngine(id);
    return seed;
  } // LArSeedService::defineEngineID()
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::ensureValidState(bool bGlobal /* = false */) const {
    if (bGlobal) {
      // registering engines may only happen in a service c'tor
      // In all other cases, throw.
      if ( (state.state() != SeedServiceHelper::ArtState::inServiceConstructor))
      {
        throw art::Exception(art::errors::LogicError)
          << "LArSeedService: not in a service constructor."
          << " May not register \"global\" engines.\n";
      }
    }
    else { // context-aware engine
      // registering engines may only happen in a c'tor
      // (disabling the ability to do that or from a beginRun method)
      // In all other cases, throw.
      if ( (state.state() != SeedServiceHelper::ArtState::inModuleConstructor)
      //  && (state.state() != SeedServiceHelper::ArtState::inModuleBeginRun)
        )
      {
        throw art::Exception(art::errors::LogicError)
          << "LArSeedService: not in a module constructor."
          << " May not register engines.\n";
      }
    } // if
  } // LArSeedService::ensureValidState()
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::reseedInstance(EngineId const& id) {
    // get all the information on the current process, event and module from
    // ArtState:
    SeedMaster_t::EventData_t const data(state.getEventSeedInputData());
    seed_t const seed = seeds.reseedEvent(id, data);
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
  
  
  void LArSeedService::reseedGlobal() {
    for (EngineId const& ID: seeds.engineIDsRange()) {
      if (!ID.isGlobal()) continue; // not global? neeext!!
      reseedInstance(ID);
    } // for
  } // LArSeedService::reseedGlobal()
  
  
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
    // Are we being called from the right place?
    ensureValidState(id.isGlobal());
    
    if (hasEngine(id)) {
      throw art::Exception(art::errors::LogicError)
        << "LArSeedService: an engine with ID '" << id.artName()
        << "' has already been created!\n";
    }
    seeds.registerNewSeeder(id, seeder);
  } // LArSeedService::registerEngineAndSeeder()
  
  
  //----------------------------------------------------------------------------
  void LArSeedService::freezeSeed(EngineId const& id, seed_t frozen_seed)
    { seeds.freezeSeed(id, frozen_seed); }
  
  
  //----------------------------------------------------------------------------
  LArSeedService::seed_t LArSeedService::prepareEngine
    (EngineId const& id, SeedMaster_t::Seeder_t seeder)
  {
    registerEngineAndSeeder(id, seeder);
    return querySeed(id);
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
    state.transit_to(SeedServiceHelper::ArtState::inModuleConstructor);
    state.set_module(md);
  } // LArSeedService::preModuleConstruction()
  
  void LArSeedService::postModuleConstruction(art::ModuleDescription const&) {
    state.reset_state();
  } // LArSeedService::postModuleConstruction()
  
  void LArSeedService::preModuleBeginRun(art::ModuleDescription const& md) {
    state.transit_to(SeedServiceHelper::ArtState::inModuleBeginRun);
    state.set_module(md);
  } // LArSeedService::preModuleBeginRun()
  
  void LArSeedService::postModuleBeginRun(art::ModuleDescription const&) {
    state.reset_state();
  } // LArSeedService::postModuleBeginRun()
  
  void LArSeedService::preProcessEvent(art::Event const& evt) {
    state.transit_to(SeedServiceHelper::ArtState::inEvent);
    state.set_event(evt);
    seeds.onNewEvent(); // inform the seed master that a new event has come
    
    LOG_DEBUG("LArSeedService") << "preProcessEvent(): will reseed global engines";
    reseedGlobal(); // why don't we do them all?!?
    
  } // LArSeedService::preProcessEvent()
  
  void LArSeedService::preModule(art::ModuleDescription const& md) {
    state.transit_to(SeedServiceHelper::ArtState::inModuleEvent);
    state.set_module(md);
    
    // Reseed all the engine of this module... maybe
    // (that is, if the current policy alows it).
    LOG_DEBUG("LArSeedService") << "preModule(): will reseed engines for module '"
      << md.moduleLabel() << "'";
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
  
  void LArSeedService::preModuleEndJob(art::ModuleDescription const& md) {
    state.transit_to(SeedServiceHelper::ArtState::inEndJob);
    state.set_module(md);
  } // LArSeedService::preModuleBeginRun()
  
  void LArSeedService::postModuleEndJob(art::ModuleDescription const&) {
    state.reset_state();
  } // LArSeedService::preModuleBeginRun()
  
  void LArSeedService::postEndJob() {
    if ((verbosity > 0) || bPrintEndOfJobSummary)
      print(); // framework logger decides whether and where it shows up
  } // LArSeedService::postEndJob()
  
  //----------------------------------------------------------------------------

} // end namespace sim

DEFINE_ART_SERVICE(sim::LArSeedService)
