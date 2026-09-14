//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYHelper.cxx
/// \brief LArSoft interface to the MARLEY (Model of Argon Reaction Low Energy
/// Yields) supernova neutrino event generator
///
/// \author Steven Gardiner <gardiner@fnal.gov>
//////////////////////////////////////////////////////////////////////////////

// framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larsim/EventGenerator/MARLEY/MARLEYHelper.h"
#include "larsim/EventGenerator/MARLEY/MarleyParameterSetWalker.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TFile.h"
#include "TInterpreter.h"
#include "TROOT.h"

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Data/GenRunInfoData.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "marley/JSONConfig.hh"
#include "marley/Reaction.hh"
#include "marley/hepmc3_utils.hh"

using ProcType = marley::Reaction::ProcessType;

namespace {
  // MARLEY's native energy units (MeV) may be converted to LArSoft's (GeV)
  // using this conversion factor
  constexpr double MeV_to_GeV = 1e-3;
}

//------------------------------------------------------------------------------
evgen::MARLEYHelper::MARLEYHelper(const fhicl::ParameterSet& pset,
                                  rndm::NuRandomService& rand_service,
                                  const std::string& helper_name,
                                  TTree* dump_tree)
  : fHelperName(helper_name), fDumpTree(dump_tree)
{
  // Configure MARLEY using the FHiCL parameters
  this->reconfigure(pset);

  // Register this MARLEY generator with the NuRandomService. For simplicity,
  // we use a lambda as the seeder function (see NuRandomService.h for
  // details). This allows the SeedService to automatically re-seed MARLEY
  // whenever necessary. The user can set an explicit seed for MARLEY in the
  // FHiCL configuration using the "seed" parameter. If you need to get the
  // seed for MARLEY from the SeedService, note that we're using use the value
  // of the input variable helper_name as its generator instance name.
  rndm::NuRandomService::seed_t marley_seed = rand_service.registerEngine(
    [this](rndm::NuRandomService::EngineId const& /* unused */,
           rndm::NuRandomService::seed_t lar_seed) -> void {
      if (fMarleyGenerator && fMarleyGenerator.get()) {
        auto seed = static_cast<uint_fast64_t>(lar_seed);
        fMarleyGenerator->reseed(seed);
      }
    },
    fHelperName,
    pset,
    {"seed"});

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
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::add_marley_particles(
  simb::MCTruth& truth,
  const std::vector<std::shared_ptr<HepMC3::GenParticle>>& particles,
  const TLorentzVector& vtx_pos,
  double conv_factor,
  bool track)
{
  // Loop over the vector of MARLEY particles and add simb::MCParticle
  // versions of each of them to the MCTruth object.
  for (const auto& p : particles) {
    // Treat all of these particles as primaries, which have negative
    // track IDs by convention
    int trackID = -1 * (truth.NParticles() + 1);

    int pdg = p->pid();
    double mass = p->generated_mass() * conv_factor;
    const auto& mom4 = p->momentum();
    double px = mom4.px() * conv_factor;
    double py = mom4.py() * conv_factor;
    double pz = mom4.pz() * conv_factor;
    double E = mom4.e() * conv_factor;
    TLorentzVector mom(px, py, pz, E);

    int status = 0; // don't track the particles in LArG4 by default
    if (track) status = 1;

    simb::MCParticle part(trackID /* trackID to use in Geant4 */,
                          pdg,
                          "MARLEY",
                          -1 /* primary particle */,
                          mass,
                          status);

    part.AddTrajectoryPoint(vtx_pos, mom);
    truth.Add(part);
  }
}

//------------------------------------------------------------------------------
simb::MCTruth evgen::MARLEYHelper::create_MCTruth(const TLorentzVector& vtx_pos)
{
  simb::MCTruth truth;

  truth.SetOrigin(simb::kSuperNovaNeutrino);

  auto event = fMarleyGenerator->create_event();

  // Double-check the energy units used in the event. If they're in MeV (the
  // usual case for MARLEY), then we will need to convert to the GeV units that
  // simb::MCTruth expects.
  double conv_factor = 1.;
  auto ev_energy_unit = event->momentum_unit();
  if (ev_energy_unit == HepMC3::Units::MEV) { conv_factor = MeV_to_GeV; }
  else if (ev_energy_unit != HepMC3::Units::GEV) {
    throw cet::exception("MARLEYHelper") << "Unrecognized energy unit"
                                         << " encountered in a MARLEY event";
  }

  // MARLEY v2 follows the NuHepMC standard
  // (https://doi.org/10.21468/SciPostPhysCodeb.57) in labeling initial
  // particles with either the "projectile" or "target" status codes. We
  // collect the projectile(s) first and then append the target(s).
  auto initial_particles =
    marley_hepmc3::get_particles_with_status(marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, *event);

  auto target_particles =
    marley_hepmc3::get_particles_with_status(marley_hepmc3::NUHEPMC_TARGET_STATUS, *event);

  for (const auto& t : target_particles)
    initial_particles.push_back(t);

  // Final-state particles have a single status code
  auto final_particles =
    marley_hepmc3::get_particles_with_status(marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, *event);

  // Add the initial and final state particles to the MCTruth object.
  add_marley_particles(truth, initial_particles, vtx_pos, conv_factor, false);
  add_marley_particles(truth, final_particles, vtx_pos, conv_factor, true);

  // calculate a few parameters for the call to SetNeutrino
  const auto& nu = marley_hepmc3::get_projectile(*event);
  const auto& p4_nu = nu->momentum();

  const auto& lep = marley_hepmc3::get_ejectile(*event);
  const auto& p4_lep = lep->momentum();

  double qt = (p4_nu.e() - p4_lep.e()) * conv_factor;
  double qx = (p4_nu.px() - p4_lep.px()) * conv_factor;
  double qy = (p4_nu.py() - p4_lep.py()) * conv_factor;
  double qz = (p4_nu.pz() - p4_lep.pz()) * conv_factor;

  double Q2 = qx * qx + qy * qy + qz * qz - qt * qt;

  // For definitions of Bjorken x, etc., a good reference is Mark Thomson's
  // set of slides on deep inelastic scattering (http://tinyurl.com/hcn5n6l)
  const auto& tgt = marley_hepmc3::get_target(*event);
  double m_tgt = tgt->generated_mass() * conv_factor;
  double bjorken_x = Q2 / (2. * m_tgt * qt);
  // Units cancel in the ratio, so no conv_factor is applied for y
  double inelasticity_y = 1. - p4_lep.e() / p4_nu.e();

  // Include the initial excitation energy of the final-state nucleus when
  // calculating W (the final-state invariant mass of the hadronic system)
  // since the other parameters (x, y) also take into account the 2-to-2
  // scattering reaction only.
  const auto& res = marley_hepmc3::get_residue(*event);
  double hadronic_mass_W = res->generated_mass() * conv_factor;

  // Retrieve the MARLEY process type code for the generated event
  ProcType proc_type = ProcType::Unknown;

  auto proc_attr = event->attribute<HepMC3::IntAttribute>("signal_process_id");
  if (proc_attr) { proc_type = marley_hepmc3::from_nuhepmc_proc_id(proc_attr->value()); }

  int cc_nc = simb::kCC;
  int mode = simb::kUnknownInteraction;
  int itype = simb::kNuanceOffset;
  int struck_nucleon_pdg = 0;

  constexpr int NUANCE_CCQE = 1;
  constexpr int NUANCE_NCEL = 2;
  constexpr int NUANCE_NuElectron = 98;
  if (proc_type == ProcType::NeutrinoCC_Discrete || proc_type == ProcType::NeutrinoCC_Continuum) {
    mode = simb::kQE;
    itype += NUANCE_CCQE; // CCQE in NUANCE labeling
    struck_nucleon_pdg = marley_utils::NEUTRON;
  }
  else if (proc_type == ProcType::AntiNeutrinoCC_Discrete ||
           proc_type == ProcType::AntiNeutrinoCC_Continuum) {
    mode = simb::kQE;
    itype += NUANCE_CCQE; // CCQE in NUANCE labeling
    struck_nucleon_pdg = marley_utils::PROTON;
  }
  else if (proc_type == ProcType::NC_Discrete || proc_type == ProcType::NC_Continuum) {
    cc_nc = simb::kNC;
    mode = simb::kQE;
    itype += NUANCE_NCEL; // NCEL in NUANCE labeling
    // Currently MARLEY doesn't label the struck nucleon for NC events (no
    // direct knockout). TODO: revisit if this changes
  }
  else if (proc_type == ProcType::NuElectronElastic) {
    mode = simb::kNuElectronElastic;
    itype += NUANCE_NuElectron; // NCEL in NUANCE labeling
  }

  truth.SetNeutrino(cc_nc,
                    mode,
                    itype,
                    tgt->pid(),
                    struck_nucleon_pdg,
                    0, // MARLEY handles low enough energies that we shouldn't need HitQuark
                    hadronic_mass_W,
                    bjorken_x,
                    inelasticity_y,
                    Q2);

  // Process the MARLEY logging messages (if any) captured by our
  // stringstream and forward them to the messagefacility logger
  std::string line;
  while (std::getline(fMarleyLogStream, line)) {
    MF_LOG_INFO(fHelperName) << line;
  }

  // Reset the MARLEY log stream
  fMarleyLogStream = std::stringstream();

  // If dumping has been enabled (indicated by a non-null fDumpTree), then ...
  if (fDumpTree) {

    // 1. Save the run information to the TFile associated with fDumpTree
    // (if it exists and the run information has not been saved previously).
    TFile* dump_file = fDumpTree->GetCurrentFile();
    if (dump_file) {
      auto ev_run_info = event->run_info();
      if (!fRunInfo && ev_run_info) {
        fRunInfo = ev_run_info;
        auto temp_run_info_data = std::make_unique<HepMC3::GenRunInfoData>();
        fRunInfo->write_data(*temp_run_info_data);

        dump_file->WriteObject(temp_run_info_data.get(), "MARLEY_run_info", "WriteDelete");
      }
      // Also a check for a drift in the run information, which should never
      // happen under correct code execution.
      else if (fRunInfo != ev_run_info) {
        throw cet::exception("MARLEYHelper") << "Unexpected change in MARLEY"
                                             << " run information";
      }
    }

    // 2. Create a branch to store the event data (if one does not already
    //    exist)
    if (!fEventData) {
      fEventData = std::make_unique<HepMC3::GenEventData>();
      // We use a bare pointer here so that the pointer-to-pointer branch
      // addressing machinery in ROOT is happy. We also use a std::unique_ptr
      // for convenient management of the associated memory.
      fEventDataPtr = fEventData.get();
      fDumpTree->Branch("event", &fEventDataPtr);
    }

    // 3. Update the event data associated with the branch (thus queueing
    //    it up for writing upon a call to TTree::Fill(), which is deferred
    //    to the caller rather than handled by MARLEYHelper itself).
    this->clear_event_data();
    event->write_data(*fEventData);
  }

  // Hand back the completed simb::MCTruth object
  return truth;
}

//------------------------------------------------------------------------------
std::string evgen::MARLEYHelper::find_file(const std::string& fileName, const std::string& fileType)
{
  cet::search_path searchPath("FW_SEARCH_PATH");

  std::string fullName;
  searchPath.find_file(fileName, fullName);

  if (fullName.empty())
    throw cet::exception("MARLEYHelper")
      << "Cannot find MARLEY " << fileType << " data file '" << fileName << '\'';

  return fullName;
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::load_full_paths_into_json(marley::JSON& json,
                                                    const std::string& key,
                                                    bool missing_ok)
{
  if (json.has_key(key)) {

    marley::JSON& value = json.at(key);

    if (value.is_array()) {
      // Replace each file name (which may appear in the FHiCL configuration
      // without a full path) with the full path found using cetlib
      for (auto& element : value.array_range()) {
        element = find_file(element.to_string(), key);
      }
    }

    else
      value = find_file(value.to_string(), key);
  }
  else if (!missing_ok)
    throw cet::exception("MARLEYHelper")
      << "Missing \"" << key << "\" key in the MARLEY parameters.";
}

//------------------------------------------------------------------------------
void evgen::MARLEYHelper::reconfigure(const fhicl::ParameterSet& pset)
{
  // Convert the FHiCL parameters into a JSON object that MARLEY can understand
  evgen::MarleyParameterSetWalker mpsw;
  pset.walk(mpsw);

  marley::JSON& json = mpsw.get_json();

  // Update the reaction and structure data file names to the full paths
  // using cetlib to search for them
  load_full_paths_into_json(json, "reactions", false);
  load_full_paths_into_json(json, "structure", true);

  // Also update the path for a neutrino source spectrum given in a ROOT
  // TFile
  if (json.has_key("source")) {
    marley::JSON& source_object = json.at("source");

    if (source_object.has_key("tfile")) { load_full_paths_into_json(source_object, "tfile"); }
  }

  // Create a new MARLEY configuration based on the JSON parameters
  MF_LOG_INFO("MARLEYHelper " + fHelperName) << "MARLEY will now use"
                                                " the JSON configuration\n"
                                             << json.dump_string() << '\n';
  marley::JSONConfig config(json);

  // Create a new marley::Generator object based on the current configuration
  fMarleyGenerator = std::make_unique<marley::Generator>(config.create_generator());
}

//------------------------------------------------------------------------------
// Removes any prior event information in the temporary storage used for dumping
void evgen::MARLEYHelper::clear_event_data()
{
  fEventData->particles.clear();
  fEventData->vertices.clear();
  fEventData->links1.clear();
  fEventData->links2.clear();
  fEventData->attribute_id.clear();
  fEventData->attribute_name.clear();
  fEventData->attribute_string.clear();
}
