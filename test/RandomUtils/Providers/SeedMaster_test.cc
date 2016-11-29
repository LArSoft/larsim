/**
 * @file   SeedMaster_test.cc
 * @brief  Test the SeedMaster object, core of NuRandomService
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 3rd, 2014
 * @see    NuRandomService.h SeedTest01_module.cc
 *
 * The test runs based on the configuration file specified on the command line.
 * 
 * Usage: SeedMaster_test  ConfigFile [ParameterSetPath [ModuleType ...]]
 * 
 * In general, the configuration files good for SeedTest01 module should be
 * good for this one as well (as long as ParameterSetPath is specified as
 * services.NuRandomService).
 * 
 * The configuration is in FHiCL format; no search of the FHiCL file is
 * performed: the path of the configuration file must be available from the
 * current directory (or you can specify it as absolute path).
 * 
 * The optional parameter set path selects which part of the configuration file
 * should be used to configure the SeedMaster. For a typical art configuration,
 * it would be "services.NuRandomService". An empty path (default) uses the
 * whole configuration directly.
 * 
 * For each specified ModuleType, all the art modules with that type are
 * considered as having test module configuration (in a art-like
 * configuration, "physics.analyzers" and "physics.producers").
 */

// C/C++ standard libraries
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm> // std::find()
#include <iostream>

// CET libraries
#include "cetlib/filepath_maker.h"

// utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/parse.h"
#include "fhiclcpp/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "canvas/Utilities/Exception.h"

// art extensions
#include "nutools/RandomUtils/Providers/SeedMaster.h"


//------------------------------------------------------------------------------
//--- stuff to facilitate the interaction with SeedMaster
//---
using seed_t = unsigned long;
using SeedMaster_t = rndm::SeedMaster<seed_t>;

/// Returns the seed for the specified module/instance, or 0 on error
seed_t ObtainSeed
  (SeedMaster_t& seeds, std::string module_name, std::string instance_name = "")
{
  // Returns the seed for the specified engine instance, or 0 in case of
  // configuration error (in which case, an error counter is increased)
  try {
    return seeds.getSeed(instance_name.empty()?
      SeedMaster_t::EngineId(module_name):
      SeedMaster_t::EngineId(module_name, instance_name)
      );
  }
  catch(art::Exception& e) {
    mf::LogError("SeedMaster") << "Caught an exception while asking seed for '"
      << module_name << "." << instance_name << ":\n"
      << e.what();
    return 0;
  }
} // ObtainSeed()


//------------------------------------------------------------------------------
//--- stuff to facilitate reading configuration
//---
inline fhicl::ParameterSet GetConfiguration
  (const fhicl::ParameterSet& full_pset, std::string key)
{
  return key.empty()? full_pset: full_pset.get<fhicl::ParameterSet>(key);
} // GetConfiguration()


inline fhicl::ParameterSet FindSeedServiceConfiguration(
  const fhicl::ParameterSet& full_pset,
  std::string key = "services.NuRandomService"
) {
  return GetConfiguration(full_pset, key);
} // FindSeedServiceConfiguration()


std::vector<fhicl::ParameterSet> FindModulesConfiguration(
  const fhicl::ParameterSet& full_pset,
  const std::vector<std::string>& module_types,
  std::string base_key = "physics.analyzers"
) {
  fhicl::ParameterSet base_pset = GetConfiguration(full_pset, base_key);
  
  // if this element has "module_type"
  std::string module_type;
  if (base_pset.get_if_present("module_type", module_type)) {
    // so this is a module configuration;
    // is this module among the ones we are looking for?
    if (std::find(module_types.begin(), module_types.end(), module_type)
      != module_types.end()
      )
      return { base_pset };
    else return {};
  } // if module configuration
  
  // not a module; maybe a module is nested somewhere in?
  std::vector<fhicl::ParameterSet> psets;
  for (std::string pset_key: base_pset.get_pset_names()) {
    std::vector<fhicl::ParameterSet> new_psets = FindModulesConfiguration(
      base_pset.get<fhicl::ParameterSet>(pset_key),
      module_types, ""
      );
    psets.insert(psets.end(), new_psets.begin(), new_psets.end());
  } // for
  return psets;
} // FindModulesConfiguration()


/// Returns a string "ModuleType[ModuleName]" out of specified module config
inline std::string GetModuleID(fhicl::ParameterSet const& pset) {
  return pset.get<std::string>("module_type", "Unknown")
    + "[" + pset.get<std::string>("module_name", "unknown") + "]";
} // GetModuleID()


//------------------------------------------------------------------------------
//--- stuff to facilitate the use of message facility
//---
void StartMessageFacility(fhicl::ParameterSet const& pset) {
  fhicl::ParameterSet mf_pset;
  if (!pset.get_if_present("services.message", mf_pset)) {
    // a destination which will react to messages more severe than DEBUG
    std::string MessageFacilityConfiguration = R"(
    destinations : {
      stdout: {
        type:      cout
        threshold: INFO
        categories: {
          default: {
            limit : 10
          }
        } // categories
      } // stdout
    } // destinations
    statistics: cout
    )";
    fhicl::make_ParameterSet(MessageFacilityConfiguration, mf_pset);
    std::cout << "Using default message facility configuration:\n"
      << mf_pset.to_indented_string(1) << std::endl;
  } // if no configuration is available
  
  mf::StartMessageFacility(mf::MessageFacilityService::SingleThread, mf_pset);
  mf::SetApplicationName("SeedMaster_test");
  mf::SetContext("Initialization");
  mf::LogInfo("MessageFacility") << "MessageFacility started.";
} // StartMessageFacility()


//------------------------------------------------------------------------------
//--- stuff to facilitate the tests
//---

/// Tests a "module" with a given configuration
int TestModule(SeedMaster_t& seeds, const fhicl::ParameterSet& pset) {
  
  mf::SetModuleName("TestModule");
  
  using seed_t = SeedMaster_t::seed_t;
  
  std::string module_name
    = pset.get<std::string>("module_name", "unknown");
  std::string module_id = GetModuleID(pset);
  
  // read the instance 
  std::vector<std::string> instance_names;
  pset.get_if_present("instanceNames", instance_names);
  unsigned int nExpectedErrors = pset.get<unsigned int>("expectedErrors", 0);
  
  unsigned int nErrors = 0;
  if (instance_names.empty()) {
    instance_names.push_back("");
    mf::LogInfo("SeedMaster_test") << "Added a default engine instance.";
  }  
  mf::LogVerbatim("SeedMaster_test")
    << module_id << " has " << instance_names.size() << " seed instances:";
  
  // seed test: do we successfully get seeds?
  std::vector<seed_t> our_seeds;
  for (std::string instance_name: instance_names) {
    seed_t seed = ObtainSeed(seeds, module_name, instance_name);
    mf::LogVerbatim("SeedMaster_test")
      << "Seed for '" << instance_name << "' is: " << seed;
    if (seed == 0) {
      LOG_ERROR(module_id)
        << "instance " << instance_name << " got seed 0!";
      if (++nErrors <= nExpectedErrors) {
        mf::LogProblem(module_id) << "  (error #" << nErrors
          << ", " << nExpectedErrors << " expected)";
      }
    }
    our_seeds.push_back(seed);
  } // for declaration loop
  
  // consistency test: are they still the same?
  std::vector<seed_t>::const_iterator iOldSeed = our_seeds.begin();
  for (std::string instance_name: instance_names) {
    seed_t seed = ObtainSeed(seeds, module_name, instance_name);
    mf::LogVerbatim("SeedMaster_test")
      << "Seed for '" << instance_name << "' is: " << seed << " (second query)";
    if (seed != *iOldSeed) {
      LOG_ERROR(module_id)
        << "seed has changed for instance "
        << instance_name << ": " << *iOldSeed << " => " << seed;
      if (++nErrors <= nExpectedErrors) {
        mf::LogProblem(module_id) << " error " << nErrors
          << ", " << nExpectedErrors << " expected";
      }
    } // if different from before
    if (*iOldSeed == 0) {
      LOG_ERROR(module_id)
        << "instance " << instance_name << " got seed 0!";
      if (++nErrors <= nExpectedErrors) {
        mf::LogProblem(module_id) << "  (error #" << nErrors
          << ", " << nExpectedErrors << " expected)";
      }
    } // if seed is 0
    ++iOldSeed;
  } // for consistency check
  
  // as many errors as expected, balance is even
  return (nErrors > nExpectedErrors)?
    nErrors - nExpectedErrors: nExpectedErrors - nErrors;
} // TestModule()



//------------------------------------------------------------------------------
//--- stuff to run the facilitated stuff
//---
int main(int argc, const char** argv) {
  
  //****************************************************************************
  //*** parse command line options
  //***
  std::string config_path;
  std::string parameter_set_name;
  std::vector<std::string> module_types{ "SeedTestPolicy" };
  
  int iParam = 0;
  // configuration file path (mandatory)
  if (++iParam < argc) {
    config_path = argv[iParam];
  }
  else {
    std::cerr
      << "SeedMaster_test: please specify a configuration file." << std::endl;
    return 1;
  }
  
  // configuration key for the SeedMaster (optional)
  if (++iParam < argc) {
    parameter_set_name = argv[iParam];
  }
  
  // module types for modules autodetection (optional)
  if (++iParam < argc) {
    module_types.assign(argv + iParam, argv + argc);
  }
  
  
  //****************************************************************************
  //*** read the configuration from the file specified on the command line
  //***
  // simple file lookup policy: assume the file name specification is complete
  cet::filepath_maker policy;

  // parse a configuration file; obtain intermediate form
  fhicl::intermediate_table table;
  fhicl::parse_document(config_path, policy, table);
  
  // translate into a parameter set
  fhicl::ParameterSet global_pset;
  fhicl::make_ParameterSet(table, global_pset);
  
  // initialize the message facility
  StartMessageFacility(global_pset);
  mf::SetModuleName("main");
  
  // read the configuration of the seed service
  fhicl::ParameterSet pset;
  try {
    pset = FindSeedServiceConfiguration(global_pset, parameter_set_name);
  }
  catch (...) {
    mf::LogError("SeedMaster_t")
      << "Failed to read the '" << parameter_set_name
      << "' parameters set from configuration file '" << config_path << "'";
    throw;
  }
  
  mf::LogInfo("SeedMaster_test") << "Parameters from '" << config_path << "':\n"
    << pset.to_indented_string();
  
  // look for test modules
  std::string module_base_key = "physics.analyzers";
  std::vector<fhicl::ParameterSet> module_psets
    = FindModulesConfiguration(global_pset, module_types, module_base_key);
  if (module_psets.empty()) {
    mf::LogError("SeedMaster_test")
      << "could not find any suitable module configuration in '"
      << module_base_key << "' of '" << config_path
      << "'";
    return 1;
  } // if no modules

  { // anonymous block
    mf::LogInfo log("SeedMaster_test");
    log << "Found configurations for " << module_psets.size() << " modules.";
    for (const fhicl::ParameterSet& pset: module_psets) {
      std::string module_id = GetModuleID(pset);
      log << '\n' << "*** " << module_id << " ***"
        << std::string(72 - 8 - module_id.length(), '*')
        << '\n' << pset.to_indented_string(1);
    } // for
    log << '\n' << std::string(72, '*');
  } // end anonymous block
  
  bool endOfJobSummary = pset.get<bool>("endOfJobSummary", false);
  
  //****************************************************************************
  //*** perform the tests...
  //***
  // create a new SeedMaster with the specified parameters set
  std::unique_ptr<SeedMaster_t> pSeeds;
  try {
    pSeeds.reset(new SeedMaster_t(pset));
  }
  catch (const art::Exception& e) {
    mf::LogError("SeedMaster_test")
      << "Exception caught while initializing SeedMaster:\n"
      << e.what();
    return 1;
  }
  catch (const fhicl::exception& e) {
    mf::LogError("SeedMaster_test")
      << "FHiCL exception caught while initializing SeedMaster:\n"
      << e.what();
    return 1;
  }
  
  unsigned int nErrors = 0;
  for (const fhicl::ParameterSet& module_pset: module_psets) {
    try {
      nErrors += TestModule(*pSeeds, module_pset);
    }
    catch(const art::Exception& e) {
      mf::LogError("SeedMaster_test")
        << "Exception caught while testing module " << GetModuleID(module_pset)
        << ":\n"
        << e.what();
      ++nErrors;
    }
  } // for
  mf::SetModuleName("main");
  
  
  //****************************************************************************
  //*** print the outcome and go
  //***
  
  if (endOfJobSummary) pSeeds->print();
  
  if (nErrors > 0) {
    mf::LogError("SeedMaster_test")
      << "Test terminated with " << nErrors << " errors.";
    return nErrors;
  }
  mf::LogInfo("SeedMaster_test") << "Test successful.";
  return 0;
} // main()
