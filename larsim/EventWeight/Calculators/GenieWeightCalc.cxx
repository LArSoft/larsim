// GenieWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated for merge into larsim develop branch on Feb 22 2021 by Steven Gardiner
// Heavily rewritten on Dec 9 2019
//   by Steven Gardiner <gardiner@fnal.gov>
// Updated by Marco Del Tutto on Feb 18 2017
// Ported from uboonecode to larsim on Feb 14 2017
//   by Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>

// Standard library includes
#include <map>
#include <memory>
#include <set>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/EventWeight/Base/WeightCalc.h"
#include "larsim/EventWeight/Base/WeightCalcCreator.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// GENIE includes
// TODO: add legacy support for GENIE v2 (mostly changes to header file
// locations would be needed)
#include "GENIE/Framework/Conventions/KineVar.h"
#include "GENIE/Framework/EventGen/EventRecord.h"
#include "GENIE/Framework/Interaction/Interaction.h"
#include "GENIE/Framework/Interaction/Kinematics.h"
#include "GENIE/Framework/Messenger/Messenger.h"
#include "GENIE/Framework/Utils/AppInit.h"

#include "GENIE/RwCalculators/GReWeightAGKY.h"
#include "GENIE/RwCalculators/GReWeightDISNuclMod.h"
#include "GENIE/RwCalculators/GReWeightFGM.h"
#include "GENIE/RwCalculators/GReWeightFZone.h"
#include "GENIE/RwCalculators/GReWeightINuke.h"
#include "GENIE/RwCalculators/GReWeightINukeParams.h"
#include "GENIE/RwCalculators/GReWeightNonResonanceBkg.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCQE.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCQEaxial.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCQEvec.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCRES.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCOH.h"
#include "GENIE/RwCalculators/GReWeightNuXSecDIS.h"
#include "GENIE/RwCalculators/GReWeightNuXSecNC.h"
#include "GENIE/RwCalculators/GReWeightNuXSecNCEL.h"
#include "GENIE/RwCalculators/GReWeightNuXSecNCRES.h"
#include "GENIE/RwCalculators/GReWeightResonanceDecay.h"
#include "GENIE/RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "GENIE/RwFramework/GReWeight.h"
#include "GENIE/RwFramework/GSyst.h"
#include "GENIE/RwFramework/GSystSet.h"

// MicroBooNE-specific reweighting tools go here. These require a special genie
// ups product.
#ifdef GENIE_UB_PATCH

// New weight calculator in GENIE v3.0.4 MicroBooNE patch 01
#include "GENIE/RwCalculators/GReWeightXSecMEC.h"

// New weight calculators in GENIE v3.0.4 MicroBooNE patch 02
#include "GENIE/RwCalculators/GReWeightDeltaradAngle.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCOHuB.h"
#include "GENIE/RwCalculators/GReWeightRESBugFix.h"

#endif

namespace {

  // These GENIE knobs are listed in the GSyst_t enum type but are not actually implemented.
  // They will be skipped and a warning message will be printed.
  // Last updated 9 Dec 2019 by S. Gardiner
  std::set<genie::rew::GSyst_t> UNIMPLEMENTED_GENIE_KNOBS = {
    kXSecTwkDial_RnubarnuCC,  // tweak the ratio of \sigma(\bar\nu CC) / \sigma(\nu CC)
    kXSecTwkDial_NormCCQEenu, // tweak CCQE normalization (maintains dependence on neutrino energy)
    kXSecTwkDial_NormDISCC,   // tweak the inclusive DIS CC normalization
    kXSecTwkDial_DISNuclMod   // unclear intent, does anyone else know? - S. Gardiner
  };

  // Some GENIE weight calculators can work with sets of knobs that should not be used simultaneously.
  // For instance, the CCQE weight calculator can vary parameters for the dipole axial form factor model
  // (the axial mass Ma) or for the z-expansion model, but using both together is invalid. The map below
  // is used to check that all of the requested knobs from the FHiCL input are compatible with each
  // other. Assuming they pass this check, the code will then configure the weight calculators to use
  // the correct "basis" of reweighting knobs as appropriate.

  // Outer keys are names of GENIE weight calculators that use sets of incompatible knobs. The names
  // need to match those used in GenieWeightCalc::SetupWeightCalculators(), specifically in the
  // calls to genie::rew::GReWeight::AdoptWeightCalc().
  // Inner keys are integer codes used to represent (and configure) each of the mutually exclusive
  // modes (i.e., to select one of the sets of incompatible knobs to use for the weight calculation).
  // Values are GSyst_t knob labels that belong to the given mode.
  std::map<std::string, std::map<int, std::set<genie::rew::GSyst_t>>> INCOMPATIBLE_GENIE_KNOBS = {
    // CCQE (genie::rew::GReWeightNuXSecCCQE)
    {"xsec_ccqe",
     {
       {genie::rew::GReWeightNuXSecCCQE::kModeNormAndMaShape,
        {
          // Norm + shape
          kXSecTwkDial_NormCCQE, // tweak CCQE normalization (energy independent)
          kXSecTwkDial_MaCCQEshape, // tweak Ma CCQE, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
          kXSecTwkDial_E0CCQEshape // tweak E0 CCQE RunningMA, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
        }},
       {genie::rew::GReWeightNuXSecCCQE::kModeMa,
        {
          // Ma
          kXSecTwkDial_MaCCQE, // tweak Ma CCQE, affects dsigma(CCQE)/dQ2 both in shape and normalization
          kXSecTwkDial_E0CCQE, // tweak E0 CCQE RunningMA, affects dsigma(CCQE)/dQ2 both in shape and normalization
        }},
       {genie::rew::GReWeightNuXSecCCQE::kModeZExp,
        {
          // Z-expansion
          kXSecTwkDial_ZNormCCQE, // tweak Z-expansion CCQE normalization (energy independent)
          kXSecTwkDial_ZExpA1CCQE, // tweak Z-expansion coefficient 1, affects dsigma(CCQE)/dQ2 both in shape and normalization
          kXSecTwkDial_ZExpA2CCQE, // tweak Z-expansion coefficient 2, affects dsigma(CCQE)/dQ2 both in shape and normalization
          kXSecTwkDial_ZExpA3CCQE, // tweak Z-expansion coefficient 3, affects dsigma(CCQE)/dQ2 both in shape and normalization
          kXSecTwkDial_ZExpA4CCQE // tweak Z-expansion coefficient 4, affects dsigma(CCQE)/dQ2 both in shape and normalization
        }},
     }},

    // CCRES (genie::rew::GReWeightNuXSecCCRES)
    {"xsec_ccres",
     {{genie::rew::GReWeightNuXSecCCRES::kModeNormAndMaMvShape,
       {
         // Norm + shape
         kXSecTwkDial_NormCCRES, /// tweak CCRES normalization
         kXSecTwkDial_MaCCRESshape, /// tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
         kXSecTwkDial_MvCCRESshape /// tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
       }},
      {genie::rew::GReWeightNuXSecCCRES::kModeMaMv,
       {
         // Ma + Mv
         kXSecTwkDial_MaCCRES, // tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
         kXSecTwkDial_MvCCRES // tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
       }}}},

    // NCRES (genie::rew::GReWeightNuXSecNCRES)
    {"xsec_ncres",
     {{genie::rew::GReWeightNuXSecNCRES::kModeNormAndMaMvShape,
       {
         // Norm + shape
         kXSecTwkDial_NormNCRES, /// tweak NCRES normalization
         kXSecTwkDial_MaNCRESshape, /// tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
         kXSecTwkDial_MvNCRESshape /// tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
       }},
      {genie::rew::GReWeightNuXSecNCRES::kModeMaMv,
       {
         // Ma + Mv
         kXSecTwkDial_MaNCRES, // tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
         kXSecTwkDial_MvNCRES // tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
       }}}},

    // DIS (genie::rew::GReWeightNuXSecDIS)
    {"xsec_dis",
     {{genie::rew::GReWeightNuXSecDIS::kModeABCV12u,
       {
         kXSecTwkDial_AhtBY, // tweak the Bodek-Yang model parameter A_{ht} - incl. both shape and normalization effect
         kXSecTwkDial_BhtBY, // tweak the Bodek-Yang model parameter B_{ht} - incl. both shape and normalization effect
         kXSecTwkDial_CV1uBY, // tweak the Bodek-Yang model parameter CV1u - incl. both shape and normalization effect
         kXSecTwkDial_CV2uBY // tweak the Bodek-Yang model parameter CV2u - incl. both shape and normalization effect
       }},
      {genie::rew::GReWeightNuXSecDIS::kModeABCV12uShape,
       {
         kXSecTwkDial_AhtBYshape, // tweak the Bodek-Yang model parameter A_{ht} - shape only effect to d2sigma(DIS)/dxdy
         kXSecTwkDial_BhtBYshape, // tweak the Bodek-Yang model parameter B_{ht} - shape only effect to d2sigma(DIS)/dxdy
         kXSecTwkDial_CV1uBYshape, // tweak the Bodek-Yang model parameter CV1u - shape only effect to d2sigma(DIS)/dxdy
         kXSecTwkDial_CV2uBYshape // tweak the Bodek-Yang model parameter CV2u - shape only effect to d2sigma(DIS)/dxdy
       }}}}};

  // Helper function that checks whether a string storing a GENIE knob name is
  // valid. Also stores the corresponding GSyst_t enum value for the knob.
  bool valid_knob_name(const std::string& knob_name, genie::rew::GSyst_t& knob)
  {
    knob = genie::rew::GSyst::FromString(knob_name);
    if (knob != kNullSystematic && knob != kNTwkDials) {
      if (UNIMPLEMENTED_GENIE_KNOBS.count(knob)) {
        MF_LOG_WARNING("GENIEWeightCalc") << "Ignoring unimplemented GENIE"
                                          << " knob " << knob_name;
        return false;
      }
    }
    else {
      throw cet::exception(__PRETTY_FUNCTION__) << "Encountered unrecognized"
                                                   "GENIE knob \""
                                                << knob_name << '\"';
      return false;
    }
    return true;
  }

  // Set of FHiCL weight calculator labels for which the tuned CV will be
  // ignored. If the name of the weight calculator doesn't appear in this set,
  // then variation weights will be thrown around the tuned CV.
  std::set<std::string> CALC_NAMES_THAT_IGNORE_TUNED_CV = {"RootinoFix"};

} // anonymous namespace

namespace evwgh {
  class GenieWeightCalc : public WeightCalc {
  public:
    GenieWeightCalc();
    void Configure(const fhicl::ParameterSet& pset, CLHEP::HepRandomEngine& engine) override;
    std::vector<std::vector<double>> GetWeight(art::Event& e) override;

  private:
    std::map<std::string, int> CheckForIncompatibleSystematics(
      const std::vector<genie::rew::GSyst_t>& knob_vec);

    void SetupWeightCalculators(genie::rew::GReWeight& rw,
                                const std::map<std::string, int>& modes_to_use);

    std::vector<genie::rew::GReWeight> reweightVector;

    std::string fGenieModuleLabel;

    bool fQuietMode;

    DECLARE_WEIGHTCALC(GenieWeightCalc)
  };

  GenieWeightCalc::GenieWeightCalc() {}

  void GenieWeightCalc::Configure(const fhicl::ParameterSet& p, CLHEP::HepRandomEngine& engine)
  {
    genie::Messenger* messenger = genie::Messenger::Instance();

    // By default, run GENIE reweighting in "quiet mode"
    // (use the logging settings in the "whisper" configuration
    // of the genie::Messenger class). The user can disable
    // quiet mode using the boolean FHiCL parameter "quiet_mode"
    fQuietMode = p.get<bool>("quiet_mode", true);
    if (fQuietMode) { genie::utils::app_init::MesgThresholds("Messenger_whisper.xml"); }
    else {
      // If quiet mode isn't enabled, then print detailed debugging
      // messages for GENIE reweighting
      messenger->SetPriorityLevel("ReW", log4cpp::Priority::DEBUG);

      MF_LOG_INFO("GENIEWeightCalc") << "Configuring GENIE weight"
                                     << " calculator " << this->GetName();
    }

    // Manually silence a couple of annoying GENIE logging messages
    // that appear a lot when running reweighting. This is done
    // whether or not "quiet mode" is enabled.
    messenger->SetPriorityLevel("TransverseEnhancementFFModel", log4cpp::Priority::WARN);
    messenger->SetPriorityLevel("Nieves", log4cpp::Priority::WARN);

    // Global Config
    fGenieModuleLabel = p.get<std::string>("genie_module_label");

    // Parameter central values from a table in the global config
    // Keys are GENIE knob names, values are knob settings that correspond to a tuned CV
    const fhicl::ParameterSet& param_CVs =
      p.get<fhicl::ParameterSet>("genie_central_values", fhicl::ParameterSet());
    std::vector<std::string> CV_knob_names = param_CVs.get_all_keys();

    // Map to store the CV knob settings
    std::map<genie::rew::GSyst_t, double> gsyst_to_cv_map;
    genie::rew::GSyst_t temp_knob;

    if (!CV_knob_names.empty() && !fQuietMode)
      MF_LOG_INFO("GENIEWeightCalc")
        << "Configuring non-default GENIE knob central values from input"
        << " FHiCL parameter set";

    for (const auto& knob_name : CV_knob_names) {
      double CV_knob_value = param_CVs.get<double>(knob_name);
      if (valid_knob_name(knob_name, temp_knob)) {
        if (gsyst_to_cv_map.count(temp_knob)) {
          throw cet::exception(__PRETTY_FUNCTION__)
            << "Duplicate central values"
            << " were configured for the " << knob_name << " GENIE knob.";
        }
        gsyst_to_cv_map[temp_knob] = CV_knob_value;
        if (!fQuietMode) {
          MF_LOG_INFO("GENIEWeightCalc")
            << "Central value for the " << knob_name << " knob was set to " << CV_knob_value;
        }
      }
    }

    // Calc Config
    const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());
    auto const pars =
      pset.get<std::vector<std::string>>("parameter_list", std::vector<std::string>());
    auto const par_sigmas = pset.get<std::vector<double>>("parameter_sigma", std::vector<double>());

    // Parameter limits (currently only used in minmax mode)
    // TODO: Revisit this. Could be useful for other modes.
    auto const par_mins = pset.get<std::vector<double>>("parameter_min", std::vector<double>());
    auto const par_maxes = pset.get<std::vector<double>>("parameter_max", std::vector<double>());

    auto const mode = pset.get<std::string>("mode");

    bool sigmas_ok = true;
    std::string array_name_for_exception;
    if (mode.find("central_value") == std::string::npos &&
        mode.find("minmax") == std::string::npos) {
      // For most reweighting modes, make sure that the number 1-sigma values
      // and the number of reweighting knobs match
      if (pars.size() != par_sigmas.size()) {
        sigmas_ok = false;
        array_name_for_exception = "parameter_sigma";
      }
    }
    else if (mode.find("minmax") != std::string::npos) {
      if (pars.size() != par_mins.size() || pars.size() != par_maxes.size()) {
        // For "minmax" mode, do the same for both the minimum and maximum
        // sigma values
        sigmas_ok = false;
        array_name_for_exception = "parameter_min and parameter_max";
      }
    }

    if (!sigmas_ok) {
      throw cet::exception(__PRETTY_FUNCTION__)
        << GetName() << "::Bad fcl configuration. parameter_list and " << array_name_for_exception
        << " need to have same number of parameters.";
    }

    if (!pars.empty() && !fQuietMode)
      MF_LOG_INFO("GENIEWeightCalc")
        << "Configuring"
        << " GENIE systematic knobs to be varied from the input FHiCL parameter set";

    // Convert the list of GENIE knob names from the input FHiCL configuration
    // into a vector of genie::rew::GSyst_t labels
    std::vector<genie::rew::GSyst_t> knobs_to_use;
    for (const auto& knob_name : pars) {
      if (valid_knob_name(knob_name, temp_knob)) knobs_to_use.push_back(temp_knob);
    }

    // We need to add all of the tuned CV knobs to the list when configuring
    // the weight calculators and checking for incompatibilities. Maybe all of
    // the systematic variation knobs are fine to use together, but they also
    // need to be compatible with the tuned CV. To perform the check, copy the
    // list of systematic variation knobs to use and add all the CV knobs that
    // aren't already present.
    std::vector<genie::rew::GSyst_t> all_knobs_vec = knobs_to_use;
    for (const auto& pair : gsyst_to_cv_map) {
      genie::rew::GSyst_t cv_knob = pair.first;
      auto begin = all_knobs_vec.cbegin();
      auto end = all_knobs_vec.cend();
      if (!std::count(begin, end, cv_knob)) all_knobs_vec.push_back(cv_knob);
    }

    // Check that the enabled knobs (both systematic variations and knobs used
    // for the CV tune) are all compatible with each other. The std::map
    // returned by this function call provides information needed to fine-tune
    // the configuration of the GENIE weight calculators.
    std::map<std::string, int> modes_to_use = this->CheckForIncompatibleSystematics(all_knobs_vec);

    // If we're working in "pm1sigma" or "minmax" mode, there should only be two
    // universes regardless of user input.
    size_t num_universes = 0u;
    if (mode.find("pm1sigma") != std::string::npos || mode.find("minmax") != std::string::npos) {
      num_universes = 2u;
    }
    else if (mode.find("multisim") != std::string::npos) {
      num_universes = pset.get<size_t>("number_of_multisims");

      // Since we're in multisim mode, force retrieval of the random
      // number seed. If it wasn't set, this will trigger an exception.
      // We want this check because otherwise the multisim universes
      // will not be easily reproducible.
      int dummy_seed = pset.get<int>("random_seed");
      MF_LOG_INFO("GENIEWeightCalc")
        << "GENIE weight calculator " << this->GetName() << " will generate " << num_universes
        << " multisim universes with random seed " << dummy_seed;
    }
    // If we're working in "central_value" or "default" mode, only a single
    // universe should be used
    else {
      num_universes = 1u;
    }

    // Create one default-constructed genie::rew::GReWeight object per universe
    reweightVector.resize(num_universes);

    // Set up the weight calculators for each universe
    for (auto& rwght : reweightVector) {
      this->SetupWeightCalculators(rwght, modes_to_use);
    }

    // Prepare sigmas
    size_t num_usable_knobs = knobs_to_use.size();
    std::vector<std::vector<double>> reweightingSigmas(num_usable_knobs);

    for (size_t k = 0u; k < num_usable_knobs; ++k) {
      reweightingSigmas[k].resize(num_universes);

      genie::rew::GSyst_t current_knob = knobs_to_use.at(k);

      for (size_t u = 0u; u < num_universes; ++u) {
        if (mode.find("multisim") != std::string::npos) {
          reweightingSigmas[k][u] = par_sigmas[k] * CLHEP::RandGaussQ::shoot(&engine, 0., 1.);
        }
        else if (mode.find("pm1sigma") != std::string::npos) {
          // u == 0 => +1*sigma; u == 1 => -1*sigma if pm1sigma is specified
          reweightingSigmas[k][u] = (u == 0 ? 1. : -1.) * par_sigmas.at(k);
        }
        else if (mode.find("minmax") != std::string::npos) {
          // u == 0 => max; u == 1 => min if minmax is specified
          reweightingSigmas[k][u] = (u == 0 ? par_maxes.at(k) : par_mins.at(k));
        }
        else if (mode.find("central_value") != std::string::npos) {
          // We'll correct for a modified CV below if needed
          reweightingSigmas[k][u] = 0.;
        }
        else {
          // By default, use the exact sigma value given for each knob
          reweightingSigmas[k][u] = par_sigmas[k];
        }

        if (!fQuietMode)
          MF_LOG_INFO("GENIEWeightCalc")
            << "Set sigma for the " << genie::rew::GSyst::AsString(current_knob)
            << " knob in universe #" << u << ". sigma = " << reweightingSigmas[k][u];

        // Add an offset if the central value for the current knob has been
        // configured (and is thus probably nonzero). Ignore this for minmax
        // mode (the limits should be chosen to respect a modified central
        // value)
        if (mode.find("minmax") == std::string::npos) {
          auto iter = gsyst_to_cv_map.find(current_knob);
          if (iter != gsyst_to_cv_map.end()) {
            reweightingSigmas[k][u] += iter->second;
            if (!fQuietMode)
              MF_LOG_INFO("GENIEWeightCalc")
                << "CV offset added to the " << genie::rew::GSyst::AsString(current_knob)
                << " knob. New sigma for universe #" << u << " is " << reweightingSigmas[k][u];
          }
        }
      }
    }

    // Don't adjust knobs to reflect the tuned CV if this weight calculator
    // should ignore those (as determined by whether it has one of the special
    // FHiCL names)
    if (!CALC_NAMES_THAT_IGNORE_TUNED_CV.count(this->GetName())) {

      // Add tuned CV knobs which have not been tweaked, and set them to their
      // modified central values. This ensures that weights are always thrown
      // around the modified CV.
      for (const auto& pair : gsyst_to_cv_map) {
        genie::rew::GSyst_t cv_knob = pair.first;
        double cv_value = pair.second;

        // If the current tuned CV knob is not present in the list of tweaked
        // knobs, then add it to the list with its tuned central value
        if (!std::count(knobs_to_use.cbegin(), knobs_to_use.cend(), cv_knob)) {
          ++num_usable_knobs;
          knobs_to_use.push_back(cv_knob);

          // The tuned CV knob will take the same value in every universe
          reweightingSigmas.emplace_back(std::vector<double>(num_universes, cv_value));

          if (!fQuietMode)
            MF_LOG_INFO("GENIEWeightCalc") << "Tuned knob " << genie::rew::GSyst::AsString(cv_knob)
                                           << " was not configured. Setting it to " << cv_value
                                           << " in all " << num_universes << " universes.";
        }
      }
    }

    // TODO: deal with parameters that have a priori bounds (e.g., FFCCQEVec,
    // which can vary on the interval [0,1])

    // Set up the knob values for each universe
    for (size_t u = 0; u < reweightVector.size(); ++u) {

      auto& rwght = reweightVector.at(u);
      genie::rew::GSystSet& syst = rwght.Systematics();

      for (unsigned int k = 0; k < knobs_to_use.size(); ++k) {
        genie::rew::GSyst_t knob = knobs_to_use.at(k);

        double twk_dial_value = reweightingSigmas.at(k).at(u);
        syst.Set(knob, twk_dial_value);

        if (!fQuietMode) {
          MF_LOG_INFO("GENIEWeightCalc") << "In universe #" << u << ", knob #" << k << " ("
                                         << genie::rew::GSyst::AsString(knob) << ") was set to"
                                         << " the value " << twk_dial_value;
        }
      } // loop over tweaked knobs

      rwght.Reconfigure();
      rwght.Print();
    } // loop over universes
  }

  // Returns a vector of weights for each neutrino interaction in the event
  std::vector<std::vector<double>> GenieWeightCalc::GetWeight(art::Event& e)
  {
    // Get the MC generator information out of the event
    // These are both handles to MC information.
    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    art::Handle<std::vector<simb::GTruth>> gTruthHandle;

    // Actually go and get the stuff
    e.getByLabel(fGenieModuleLabel, mcTruthHandle);
    e.getByLabel(fGenieModuleLabel, gTruthHandle);

    std::vector<art::Ptr<simb::MCTruth>> mclist;
    art::fill_ptr_vector(mclist, mcTruthHandle);

    std::vector<art::Ptr<simb::GTruth>> glist;
    art::fill_ptr_vector(glist, gTruthHandle);

    size_t num_neutrinos = mclist.size();
    size_t num_knobs = reweightVector.size();

    // Calculate weight(s) here
    std::vector<std::vector<double>> weights(num_neutrinos);
    for (size_t v = 0u; v < num_neutrinos; ++v) {

      // Convert the MCTruth and GTruth objects from the event
      // back into the original genie::EventRecord needed to
      // compute the weights
      std::unique_ptr<genie::EventRecord> genie_event(evgb::RetrieveGHEP(*mclist[v], *glist[v]));

      // Set the final lepton kinetic energy and scattering cosine
      // in the owned GENIE kinematics object. This is done during
      // event generation but is not reproduced by evgb::RetrieveGHEP().
      // Several new CCMEC weight calculators developed for MicroBooNE
      // expect the variables to be set in this way (so that differential
      // cross sections can be recomputed). Failing to set them results
      // in inf and NaN weights.
      // TODO: maybe update evgb::RetrieveGHEP to handle this instead.
      genie::Interaction* interaction = genie_event->Summary();
      genie::Kinematics* kine_ptr = interaction->KinePtr();

      // Final lepton mass
      double ml = interaction->FSPrimLepton()->Mass();
      // Final lepton 4-momentum
      const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
      // Final lepton kinetic energy
      double Tl = p4l.E() - ml;
      // Final lepton scattering cosine
      double ctl = p4l.CosTheta();

      kine_ptr->SetKV(kKVTl, Tl);
      kine_ptr->SetKV(kKVctl, ctl);

      // All right, the event record is fully ready. Now ask the GReWeight
      // objects to compute the weights.
      weights[v].resize(num_knobs);
      for (size_t k = 0u; k < num_knobs; ++k) {
        weights[v][k] = reweightVector.at(k).CalcWeight(*genie_event);
      }
    }
    return weights;
  }

  std::map<std::string, int> GenieWeightCalc::CheckForIncompatibleSystematics(
    const std::vector<genie::rew::GSyst_t>& knob_vec)
  {
    std::map<std::string, int> modes_to_use;

    for (const auto& knob : knob_vec) {
      for (const auto& pair1 : INCOMPATIBLE_GENIE_KNOBS) {
        std::string calc_name = pair1.first;
        const auto& mode_map = pair1.second;
        for (const auto& pair2 : mode_map) {
          int mode = pair2.first;
          std::set<genie::rew::GSyst_t> knob_set = pair2.second;

          if (knob_set.count(knob)) {
            auto search = modes_to_use.find(calc_name);
            if (search != modes_to_use.end()) {
              if (search->second != mode) {
                auto knob_str = genie::rew::GSyst::AsString(knob);
                throw cet::exception(__PRETTY_FUNCTION__)
                  << this->GetName() << ": the GENIE knob " << knob_str << " is incompatible"
                  << " with others that are already configured";
              }
            }
            else
              modes_to_use[calc_name] = mode;
          }
        }
      }
    }

    return modes_to_use;
  }

  void GenieWeightCalc::SetupWeightCalculators(genie::rew::GReWeight& rw,
                                               const std::map<std::string, int>& modes_to_use)
  {
    // Based on the list from the GENIE command-line tool grwght1p
    rw.AdoptWghtCalc("xsec_ncel", new GReWeightNuXSecNCEL);
    rw.AdoptWghtCalc("xsec_ccqe", new GReWeightNuXSecCCQE);
    rw.AdoptWghtCalc("xsec_ccqe_axial", new GReWeightNuXSecCCQEaxial);
    rw.AdoptWghtCalc("xsec_ccqe_vec", new GReWeightNuXSecCCQEvec);
    rw.AdoptWghtCalc("xsec_ccres", new GReWeightNuXSecCCRES);
    rw.AdoptWghtCalc("xsec_ncres", new GReWeightNuXSecNCRES);
    rw.AdoptWghtCalc("xsec_nonresbkg", new GReWeightNonResonanceBkg);
    rw.AdoptWghtCalc("xsec_coh", new GReWeightNuXSecCOH);
    rw.AdoptWghtCalc("xsec_dis", new GReWeightNuXSecDIS);
    rw.AdoptWghtCalc("nuclear_qe", new GReWeightFGM);
    rw.AdoptWghtCalc("hadro_res_decay", new GReWeightResonanceDecay);
    rw.AdoptWghtCalc("hadro_fzone", new GReWeightFZone);
    rw.AdoptWghtCalc("hadro_intranuke", new GReWeightINuke);
    rw.AdoptWghtCalc("hadro_agky", new GReWeightAGKY);
    rw.AdoptWghtCalc("xsec_nc", new GReWeightNuXSecNC);
    rw.AdoptWghtCalc("res_dk", new GReWeightResonanceDecay);
    rw.AdoptWghtCalc("xsec_empmec", new GReWeightXSecEmpiricalMEC);
    // GReWeightDISNuclMod::CalcWeight() is not implemented, so we won't
    // bother to use it here. - S. Gardiner, 9 Dec 2019
    //rw.AdoptWghtCalc( "nuclear_dis",   new GReWeightDISNuclMod );

#ifdef GENIE_UB_PATCH

    // New weight calculator in GENIE v3.0.4 MicroBooNE patch 01
    rw.AdoptWghtCalc("xsec_mec", new GReWeightXSecMEC);

    // New weight calculators in GENIE v3.0.4 MicroBooNE patch 02
    rw.AdoptWghtCalc("deltarad_angle", new GReWeightDeltaradAngle);
    rw.AdoptWghtCalc("xsec_coh_ub", new GReWeightNuXSecCOHuB);
    rw.AdoptWghtCalc("res_bug_fix", new GReWeightRESBugFix);

#endif

    // Set the modes for the weight calculators that need them to be specified
    for (const auto& pair : modes_to_use) {
      std::string calc_name = pair.first;
      int mode = pair.second;

      genie::rew::GReWeightI* calc = rw.WghtCalc(calc_name);
      if (!calc)
        throw cet::exception(__PRETTY_FUNCTION__)
          << "Failed to retrieve the GENIE weight calculator labeled \"" << calc_name << '\"';

      // The GReWeightI base class doesn't have a SetMode(int) function,
      // so we'll just try dynamic casting until we get the right one.
      // If none work, then throw an exception.
      // TODO: Add a virtual function GReWeightI::SetMode( int ) in GENIE's
      // Reweight framework. Then we can avoid the hacky dynamic casts here.
      auto* calc_ccqe = dynamic_cast<genie::rew::GReWeightNuXSecCCQE*>(calc);
      auto* calc_ccres = dynamic_cast<genie::rew::GReWeightNuXSecCCRES*>(calc);
      auto* calc_ncres = dynamic_cast<genie::rew::GReWeightNuXSecNCRES*>(calc);
      auto* calc_dis = dynamic_cast<genie::rew::GReWeightNuXSecDIS*>(calc);
      if (calc_ccqe)
        calc_ccqe->SetMode(mode);
      else if (calc_ccres)
        calc_ccres->SetMode(mode);
      else if (calc_ncres)
        calc_ncres->SetMode(mode);
      else if (calc_dis)
        calc_dis->SetMode(mode);
      else
        throw cet::exception(__PRETTY_FUNCTION__)
          << "Request to set the mode of an unrecognized GENIE weight calculator \"" << calc_name
          << '\"';
    }
  }

  REGISTER_WEIGHTCALC(GenieWeightCalc)
}
