//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYHelper.cxx
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

// framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larsim/EventGenerator/MARLEY/MARLEYHelper.h"

// ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"

// MARLEY includes
#include "marley/Event.hh"
#include "marley/RootJSONConfig.hh"

namespace {
  // We need to convert from MARLEY's energy units (MeV) to LArSoft's
  // (GeV) using this conversion factor
  constexpr double MeV_to_GeV = 1e-3;
}

//------------------------------------------------------------------------------
evgen::MARLEYHelper::MARLEYHelper(
  const fhicl::Table<evgen::MARLEYHelper::Config>& conf,
  rndm::NuRandomService& rand_service, const std::string& helper_name)
  : fHelperName(helper_name)
{
  // Configure MARLEY using the FHiCL parameters
  this->reconfigure(conf);

  // Register this MARLEY generator with the NuRandomService. For simplicity,
  // we use a lambda as the seeder function (see NuRandomService.h for
  // details). This allows the SeedService to automatically re-seed MARLEY
  // whenever necessary. The user can set an explicit seed for MARLEY in the
  // FHiCL configuration using the "seed" parameter. If you need to get the
  // seed for MARLEY from the SeedService, note that we're using use the value
  // of the input variable helper_name as its generator instance name.
  rndm::NuRandomService::seed_t marley_seed = rand_service.registerEngine(
    [this](rndm::NuRandomService::EngineId const& /* unused */,
      rndm::NuRandomService::seed_t lar_seed) -> void
    {
      if (fMarleyGenerator && fMarleyGenerator.get()) {
        auto seed = static_cast<uint_fast64_t>(lar_seed);
        fMarleyGenerator.get()->reseed(seed);
      }
    },
    fHelperName, conf.get_PSet(), { "seed" }
  );

  // Unless I'm mistaken, the call to registerEngine should seed the generator
  // with the seed from the FHiCL configuration file if one is included, but it
  // doesn't appear to do so (as of 16 Aug 2016, larsoft v06_03_00). As a
  // workaround, I manually reseed the generator (if needed) here using the
  // result of the call to registerEngine, which will be the seed from the
  // FHiCL file if one was given.
  // TODO: figure out what's going on here, and remove this workaround as
  // needed
  uint_fast64_t marley_cast_seed = static_cast<uint_fast64_t>(marley_seed);
  if (marley_cast_seed != fMarleyGenerator->get_seed()) {
    fMarleyGenerator->reseed(marley_cast_seed);
  }

  // Log initialization information from the MARLEY generator
  MF_LOG_INFO(fHelperName) << fMarleyLogStream.str();
  fMarleyLogStream = std::stringstream();

  // Do any needed setup of the MARLEY class dictionaries
  load_marley_dictionaries();
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::add_marley_particles(simb::MCTruth& truth,
  const std::vector<marley::Particle*>& particles,
  const TLorentzVector& vtx_pos, bool track)
{
  // Loop over the vector of MARLEY particles and add simb::MCParticle
  // versions of each of them to the MCTruth object.
  for (const marley::Particle* p : particles) {
    // Treat all of these particles as primaries, which have negative
    // track IDs by convention
    int trackID = -1*(truth.NParticles() + 1);

    int pdg = p->pdg_code();
    double mass = p->mass() * MeV_to_GeV;
    double px = p->px() * MeV_to_GeV;
    double py = p->py() * MeV_to_GeV;
    double pz = p->pz() * MeV_to_GeV;
    double E = p->total_energy() * MeV_to_GeV;
    TLorentzVector mom(px, py, pz, E);

    int status = 0; // don't track the particles in LArG4 by default
    if (track) status = 1;

    simb::MCParticle part(trackID /* trackID to use in Geant4 */, pdg,
      "MARLEY", -1 /* primary particle */, mass, status);

    part.AddTrajectoryPoint(vtx_pos, mom);
    truth.Add(part);
  }
}

//------------------------------------------------------------------------------
simb::MCTruth evgen::MARLEYHelper::create_MCTruth(
  const TLorentzVector& vtx_pos, marley::Event* marley_event)
{
  simb::MCTruth truth;

  truth.SetOrigin(simb::kSuperNovaNeutrino);

  marley::Event event = fMarleyGenerator->create_event();

  // Add the initial and final state particles to the MCTruth object.
  add_marley_particles(truth, event.get_initial_particles(), vtx_pos, false);
  add_marley_particles(truth, event.get_final_particles(), vtx_pos, true);

  // calculate a few parameters for the call to SetNeutrino
  const marley::Particle& nu = event.projectile();
  const marley::Particle& lep = event.ejectile();
  double qx = nu.px() - lep.px();
  double qy = nu.py() - lep.py();
  double qz = nu.pz() - lep.pz();
  double Enu = nu.total_energy();
  double Elep = lep.total_energy();
  double Q2 = qx*qx + qy*qy + qz*qz - std::pow(Enu - Elep, 2);

  // For definitions of Bjorken x, etc., a good reference is Mark Thomson's
  // set of slides on deep inelastic scattering (http://tinyurl.com/hcn5n6l)
  double bjorken_x = Q2 / (2 * event.target().mass() * (Enu - Elep));
  double inelasticity_y = 1. - Elep / Enu;

  // Include the initial excitation energy of the final-state nucleus when
  // calculating W (the final-state invariant mass of the hadronic system)
  // since the other parameters (X, Y) also take into account the 2-2
  // scattering reaction only.
  const marley::Particle& res = event.residue();
  double hadronic_mass_W = res.mass() + event.Ex();

  // TODO: do a more careful job of setting the parameters here
  truth.SetNeutrino(
    simb::kCC, // change when MARLEY can handle NC
    simb::kUnknownInteraction, // not sure what the mode should be
    simb::kUnknownInteraction, // not sure what the interaction type should be
    marley_utils::get_nucleus_pid(18, 40), // Ar-40 PDG code
    marley_utils::NEUTRON, // nucleon PDG
    0, // MARLEY handles low enough energies that we shouldn't need HitQuark
    hadronic_mass_W * MeV_to_GeV,
    bjorken_x, // dimensionless
    inelasticity_y, // dimensionless
    Q2 * std::pow(MeV_to_GeV, 2)
  );

  if (marley_event) *marley_event = event;

  // Process the MARLEY logging messages (if any) captured by our
  // stringstream and forward them to the messagefacility logger
  std::string line;
  while(std::getline(fMarleyLogStream, line)) {
    MF_LOG_INFO(fHelperName) << line;
  }

  // Reset the MARLEY log stream
  fMarleyLogStream = std::stringstream();

  return truth;
}

//------------------------------------------------------------------------------
std::string evgen::MARLEYHelper::find_file(const std::string& fileName,
  const std::string& fileType)
{
  cet::search_path searchPath("FW_SEARCH_PATH");

  std::string fullName;
  searchPath.find_file(fileName, fullName);

  if (fullName.empty())
    throw cet::exception("MARLEYHelper")
      << "Cannot find MARLEY " << fileType << " data file '"
      << fileName << '\'';

  return fullName;
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::load_full_paths_into_json(
  marley::JSON& json, const std::string& key)
{
  if ( json.has_key(key) ) {

    marley::JSON& value = json.at(key);

    if ( value.is_array() ) {
      // Replace each file name (which may appear in the FHiCL configuration
      // without a full path) with the full path found using cetlib
      for (auto& element : value.array_range()) {
        element = find_file(element.to_string(), key);
      }
    }

    else value = find_file(value.to_string(), key);
  }
  else throw cet::exception("MARLEYHelper") << "Missing \"" << key
    << "\" key in the MARLEY parameters.";
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::reconfigure(
  const fhicl::Table<evgen::MARLEYHelper::Config>& conf)
{
  // Convert the FHiCL parameters into a JSON object that MARLEY can understand
  marley::JSON json = fhicl_parameter_to_json(&conf);

  // Update the reaction and structure data file names to the full paths
  // using cetlib to search for them
  load_full_paths_into_json(json, "reactions");
  load_full_paths_into_json(json, "structure");

  // Also update the path for a neutrino source spectrum given in a ROOT
  // TFile
  if ( json.has_key("source") ) {
    marley::JSON& source_object = json.at("source");

    if ( source_object.has_key("tfile") ) {
      load_full_paths_into_json(source_object, "tfile");
    }
  }

  // Convert the neutrino direction 3-vector into the JSON object expected by
  // MARLEY
  if ( json.has_key("direction") ) {
    // Get the vector components from the length-3 array
    marley::JSON& dir_object = json.at("direction");

    double x = dir_object.at(0).to_double();
    double y = dir_object.at(1).to_double();
    double z = dir_object.at(2).to_double();

    // Replace the array with an equivalent JSON object
    dir_object = marley::JSON();
    dir_object["x"] = x;
    dir_object["y"] = y;
    dir_object["z"] = z;
  }

  // Create a new MARLEY configuration based on the JSON parameters
  MF_LOG_INFO("MARLEYHelper " + fHelperName) << "MARLEY will now use"
    " the JSON configuration\n" << json.dump_string() << '\n';
  marley::RootJSONConfig config(json);

  // Create a new marley::Generator object based on the current configuration
  fMarleyGenerator = std::make_unique<marley::Generator>(
    config.create_generator());
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::load_marley_dictionaries()
{

  static bool already_loaded_marley_dict = false;

  if (already_loaded_marley_dict) return;

  // Current (24 July 2016) versions of ROOT 6 require runtime
  // loading of headers for custom classes in order to use
  // dictionaries correctly. If we're running ROOT 6+, do the
  // loading here, and give the user guidance if there are any
  // problems.
  //
  // This is the same technique used in the MARLEY source code
  // for the executable (src/marley.cc). If you change how this
  // code works, please sync changes with the executable as well.
  if (gROOT->GetVersionInt() >= 60000) {
    MF_LOG_INFO("MARLEYHelper " + fHelperName) << "ROOT 6 or greater"
      << " detected. Loading class information\nfrom headers"
      << " \"marley/Particle.hh\" and \"marley/Event.hh\"";
    TInterpreter::EErrorCode* ec = new TInterpreter::EErrorCode();
    gInterpreter->ProcessLine("#include \"marley/Particle.hh\"", ec);
    if (*ec != 0) throw cet::exception("MARLEYHelper " + fHelperName)
      << "Error loading MARLEY header Particle.hh. For MARLEY headers stored"
      << " in /path/to/include/marley/, please add /path/to/include"
      << " to your ROOT_INCLUDE_PATH environment variable and"
      << " try again.";
    gInterpreter->ProcessLine("#include \"marley/Event.hh\"");
    if (*ec != 0) throw cet::exception("MARLEYHelper") << "Error loading"
      << " MARLEY header Event.hh. For MARLEY headers stored in"
      << " /path/to/include/marley/, please add /path/to/include"
      << " to your ROOT_INCLUDE_PATH environment variable and"
      << " try again.";
  }

  // No further action is required for ROOT 5 because the compiled
  // dictionaries (which are linked to this algorithm) contain all of
  // the needed information
  already_loaded_marley_dict = true;
}

//------------------------------------------------------------------------------
marley::JSON evgen::MARLEYHelper::fhicl_parameter_to_json(
  const fhicl::detail::ParameterBase* par)
{
  // Don't bother with the rest of the analysis if we've been handed a
  // nullptr or if the parameter is disabled in the current configuration
  if (par && par->should_use()) {

    auto type = par->parameter_type();
    if (fhicl::is_atom(type)) {

      // This is an ugly hack to load FHiCL atoms into JSON objects. We have
      // to do it type-by-type using the current (12/2016) implementation of
      // the fhicl::Atom and fhicl::OptionalAtom template classes.
      if ( !par->is_optional() ) {
        if (auto atm = dynamic_cast<const fhicl::Atom<double>* >(par))
          return fhicl_atom_to_json<>(atm);
        else if (auto atm = dynamic_cast<const fhicl::Atom<std::string>* >(par))
        {
          return fhicl_atom_to_json<>(atm);
        }
        else if (auto atm = dynamic_cast<const fhicl::Atom<bool>* >(par))
        {
          return fhicl_atom_to_json<>(atm);
        }
        else if (auto atm = dynamic_cast<const fhicl::Atom<int>* >(par))
        {
          return fhicl_atom_to_json<>(atm);
        }
        else throw cet::exception("MarleyGen") << "Failed to deduce"
          << " the type of the FHiCL atom " << par->key();
      }
      else { // optional atoms
        if (auto opt_atom = dynamic_cast<const fhicl::OptionalAtom<double>* >(par))
        {
          return fhicl_optional_atom_to_json<>(opt_atom);
        }
        else if (auto opt_atom = dynamic_cast<const fhicl::OptionalAtom<std::string>* >(par))
        {
          return fhicl_optional_atom_to_json<>(opt_atom);
        }
        else if (auto opt_atom = dynamic_cast<const fhicl::OptionalAtom<bool>* >(par))
        {
          return fhicl_optional_atom_to_json<>(opt_atom);
        }
        else if (auto opt_atom = dynamic_cast<const fhicl::OptionalAtom<int>* >(par))
        {
          return fhicl_optional_atom_to_json<>(opt_atom);
        }
        else throw cet::exception("MarleyGen") << "Failed to deduce"
          << " the type of the FHiCL optional atom " << par->key();
      }
    }
    else if (fhicl::is_sequence(type)) {

      // This is another ugly hack to allow us to fill the JSON array. We
      // have to do it type-by-type using the current (12/2016)
      // implementation of the fhicl::Sequence template class.
      if (auto seq = dynamic_cast<const fhicl::Sequence<double>* >(par))
        return fhicl_sequence_to_json<>(seq);
      else if (auto seq = dynamic_cast<const fhicl::Sequence<double, 3>* >(par))
      {
        return fhicl_sequence_to_json<>(seq);
      }
      else if (auto seq = dynamic_cast<const fhicl::Sequence<std::string>* >(par))
      {
        return fhicl_sequence_to_json<>(seq);
      }
      else if (auto seq = dynamic_cast<const fhicl::Sequence<int>* >(par))
      {
        return fhicl_sequence_to_json<>(seq);
      }
      else if (auto seq = dynamic_cast<const fhicl::Sequence<bool>* >(par))
      {
        return fhicl_sequence_to_json<>(seq);
      }
      else throw cet::exception("MarleyGen") << "Failed to deduce"
        << " the type of the FHiCL sequence " << par->key();

      return marley::JSON();
    }
    else if (fhicl::is_table(type)) {

      // Downcast the parameter pointer into a TableBase pointer so that we
      // can iterate over its members
      auto table = dynamic_cast<const fhicl::detail::TableBase*>(par);
      if (!table) throw cet::exception("MarleyGen") << "Failed dynamic_cast"
        << " of fhicl::detail::ParameterBase* to fhicl::detail::TableBase*"
        << " when fhicl_is_table() was true";

      // Create an empty JSON object
      marley::JSON object = marley::JSON::make(
        marley::JSON::DataType::Object);

      // Load it with the members of the table recursively
      for (const auto& m : table->members()) {
        if (m && m->should_use()) {
          marley::JSON temp = fhicl_parameter_to_json(m.get());
          // Skip members that are null (typically unused optional FHiCL
          // parameters)
          if (!temp.is_null()) object[ m->name() ] = temp;
        }
      }
      return object;
    }
  }

  // Return a null JSON object if something strange happened.
  return marley::JSON();
}
