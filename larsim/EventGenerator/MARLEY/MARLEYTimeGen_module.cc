//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYTimeGen_module.cc
/// \brief Module that allows for sampling neutrino energies and times from
/// time-dependent supernova spectra. This module uses MARLEY to help generate
/// events.
///
/// \author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

// standard library includes
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Simulation/SupernovaTruth.h"
#include "larsim/EventGenerator/MARLEY/MARLEYGenerator.h"
#include "larsim/EventGenerator/MARLEY/ActiveVolumeVertexSampler.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TVectorD.h"

// MARLEY includes
#include "marley/marley_root.hh"
#include "marley/marley_utils.hh"

// Anonymous namespace for definitions local to this source file
namespace {
  constexpr int MAX_UNIFORM_ENERGY_ITERATIONS = 1000;
}

namespace evgen {
  class MarleyTimeGen;
}

class evgen::MarleyTimeGen : public art::EDProducer {

  public:
  
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
  
    /// Collection of configuration parameters for the module
    struct Config {
  
      fhicl::Table<evgen::ActiveVolumeVertexSampler::Config> vertex_ {
        Name("vertex"),
        Comment("Configuration for selecting the vertex location(s)")
      };
  
      fhicl::Table<evgen::MARLEYGenerator::Config> marley_parameters_ {
        Name("marley_parameters"),
        Comment("Configuration for the MARLEY generator")
      };
  
      fhicl::Atom<std::string> module_type_ {
        Name("module_type"),
        Comment(""),
        "MARLEYTimeGen" // default value
      };

      fhicl::Atom<std::string> sampling_mode_ {
        Name("sampling_mode"),
        Comment("Technique to use when sampling times and energies. Valid"
          " options are \"histogram\", \"uniform time\","
          " and \"uniform energy\""),
        std::string("histogram") // default value
      };

      fhicl::Atom<unsigned int> nu_per_event_ {
        Name("nu_per_event"),
        Comment("The number of neutrino vertices to generate in"
          " each art::Event"),
        1u // default value
      };

      fhicl::Atom<std::string> spectrum_file_ {
        Name("spectrum_file"),
        Comment("Name of a ROOT file containing a TH2D object to use"
          " to represent the incident (not cross section weighted)"
          " neutrino spectrum as a function of time (X axis, seconds) and"
          " energy (Y axis, MeV).")
      };

      fhicl::Atom<std::string> namecycle_ {
        Name("namecycle"),
        Comment("Name of the TH2D object to use to represent the"
          " incident neutrino spectrum. This value should match the"
          " name of the TH2D as given in the ROOT file specified"
          " in the \"spectrum_file\" parameter.")
      };
  
    }; // struct Config
  
    // Type to enable FHiCL parameter validation by art
    using Parameters = art::EDProducer::Table<Config>;
  
    // Configuration-checking constructors
    explicit MarleyTimeGen(const Parameters& p);
  
    virtual ~MarleyTimeGen();
  
    virtual void produce(art::Event& e) override;
    virtual void beginRun(art::Run& run) override;
  
    virtual void reconfigure(const Parameters& p);
  
  protected:
  
    // Object that provides an interface to the MARLEY event generator
    std::unique_ptr<evgen::MARLEYGenerator> fMarleyGenerator;
  
    // Algorithm that allows us to sample vertex locations within the active
    // volume(s) of the detector
    std::unique_ptr<evgen::ActiveVolumeVertexSampler> fVertexSampler;
  
    // unique_ptr to the current event created by MARLEY
    std::unique_ptr<marley::Event> fEvent;

    // ROOT TH2D that contains the time-dependent spectrum to use when sampling
    // neutrino times and energies
    std::unique_ptr<TH2D> fSpectrumHist;

    // Sampling mode to use when selecting neutrino times and energies
    enum class TimeGenSamplingMode { HISTOGRAM, UNIFORM_TIME, UNIFORM_ENERGY };
    TimeGenSamplingMode fSamplingMode;
  
    // the MARLEY event TTree
    TTree* fEventTree;
  
    // Run, subrun, and event numbers from the art::Event being processed
    uint_fast32_t fRunNumber;
    uint_fast32_t fSubRunNumber;
    uint_fast32_t fEventNumber;

    // Time since start of supernova for the current event
    double fTNu;

    // Weight for the current event
    double fWeight;

    // Flux-averaged total cross section (over all energies and times for
    // all defined reactions). This will be used to assign weights to each
    // neutrino vertex (simb::MCTruth object).
    double fFluxAveragedCrossSection; // fm^2

    // The number of neutrino vertices to generate in each art::Event
    unsigned int fNeutrinosPerEvent;
};

//------------------------------------------------------------------------------
evgen::MarleyTimeGen::MarleyTimeGen(const Parameters& p)
  : fEvent(new marley::Event), fRunNumber(0), fSubRunNumber(0),
  fEventNumber(0), fTNu(0.), fFluxAveragedCrossSection(0.)
{
  // Configure the module (including MARLEY itself) using the FHiCL parameters
  this->reconfigure(p);

  // Create a ROOT TTree using the TFileService that will store the MARLEY
  // event objects (useful for debugging purposes)
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("MARLEY Event Tree",
    "A tree of marley::Event objects");
  fEventTree->Branch("events", "marley::Event", fEvent.get());

  // Add branches that give the art::Event run, subrun, and event numbers for
  // easy match-ups between the MARLEY and art TTrees. All three are recorded
  // as 32-bit unsigned integers.
  fEventTree->Branch("run_number", &fRunNumber, "run_number/i");
  fEventTree->Branch("subrun_number", &fSubRunNumber, "subrun_number/i");
  fEventTree->Branch("event_number", &fEventNumber, "event_number/i");
  fEventTree->Branch("tSN", &fTNu, "tSN/D");
  fEventTree->Branch("weight", &fWeight, "weight/D");

  produces< std::vector<simb::MCTruth> >();
  produces< std::vector<sim::SupernovaTruth> >();
  produces< art::Assns<simb::MCTruth, sim::SupernovaTruth> >();
  produces< sumdata::RunData, art::InRun >();
}

//------------------------------------------------------------------------------
evgen::MarleyTimeGen::~MarleyTimeGen()
{
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::beginRun(art::Run& run) {
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;
  std::unique_ptr<sumdata::RunData>
    runcol(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runcol));
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::produce(art::Event& e)
{
  art::ServiceHandle<geo::Geometry> geo;

  // Get the run, subrun, and event numbers from the current art::Event
  fRunNumber = e.run();
  fSubRunNumber = e.subRun();
  fEventNumber = e.event();

  // Prepare associations and vectors of truth objects that will be produced
  // and loaded into the current art::Event
  std::unique_ptr< std::vector<simb::MCTruth> >
    truthcol(new std::vector<simb::MCTruth>);

  std::unique_ptr< std::vector<sim::SupernovaTruth> >
    sn_truthcol(new std::vector<sim::SupernovaTruth>);

  std::unique_ptr< art::Assns<simb::MCTruth, sim::SupernovaTruth> >
    truth_assns(new art::Assns<simb::MCTruth, sim::SupernovaTruth>);

  // Create temporary truth objects that we will use to load the event
  simb::MCTruth truth;
  sim::SupernovaTruth sn_truth;

  // Get a reference to the generator object created by MARLEY (we'll need
  // to do a few fancy things with it other than just creating events)
  marley::Generator& gen = fMarleyGenerator->get_generator();

  for (unsigned int n = 0; n < fNeutrinosPerEvent; ++n) {

    // Sample a primary vertex location for this event
    TLorentzVector vertex_pos = fVertexSampler->sample_vertex_pos(*geo);

    // Reset the neutrino's time-since-supernova to zero (for now)
    fTNu = 0.;

    if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM)
    {
      // Generate a MARLEY event using the time-integrated spectrum
      // (the generator was already configured to use it by reconfigure())
      truth = fMarleyGenerator->create_MCTruth(vertex_pos,
        fEvent.get());

      // Find the time distribution corresponding to the selected energy bin
      double E_nu = fEvent->projectile().total_energy();
      int E_bin_index = fSpectrumHist->GetYaxis()->FindBin(E_nu);
      TH1D* t_hist = fSpectrumHist->ProjectionX("dummy_time_hist", E_bin_index,
        E_bin_index);
      double* time_bin_weights = t_hist->GetArray();

      // Sample a time bin from the distribution
      std::discrete_distribution<int> time_dist;
      std::discrete_distribution<int>::param_type time_params(time_bin_weights,
        &(time_bin_weights[t_hist->GetNbinsX() + 1]));
      int time_bin_index = gen.discrete_sample(time_dist, time_params);

      // Sample a time uniformly from within the selected time bin
      double t_min = t_hist->GetBinLowEdge(time_bin_index);
      double t_max = t_min + t_hist->GetBinWidth(time_bin_index);
      // sample a time on [ t_min, t_max )
      fTNu = gen.uniform_random_double(t_min, t_max, false);
      fWeight = fFluxAveragedCrossSection;

      sn_truth = sim::SupernovaTruth(fTNu, fWeight, sim::kUnbiased);
    }

    else if (fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
    {
      // Generate a MARLEY event using the time-integrated spectrum
      // (the generator was already configured to use it by reconfigure())
      truth = fMarleyGenerator->create_MCTruth(vertex_pos,
        fEvent.get());

      // Sample a time uniformly
      TAxis* time_axis = fSpectrumHist->GetXaxis();
      // underflow bin has index zero
      double t_min = time_axis->GetBinLowEdge(1);
      double t_max = time_axis->GetBinLowEdge(time_axis->GetNbins() + 1);
      // sample a time on [ t_min, t_max )
      fTNu = gen.uniform_random_double(t_min, t_max, false);

      // Get the value of the true dependent probability density (probability
      // of the sampled time given the sampled energy) to use as a biasing
      // correction in the neutrino vertex weight.
      double E_nu = fEvent->projectile().total_energy();
      int E_bin_index = fSpectrumHist->GetYaxis()->FindBin(E_nu);
      TH1D* t_hist = fSpectrumHist->ProjectionX("dummy_time_hist", E_bin_index,
        E_bin_index);
      int t_bin_index = t_hist->FindBin(fTNu);
      double weight_bias = t_hist->GetBinContent(t_bin_index)
        / ( t_hist->Integral() * t_hist->GetBinWidth(t_bin_index) 
        * (t_max - t_min) );

      fWeight = fFluxAveragedCrossSection * weight_bias;

      sn_truth = sim::SupernovaTruth(fTNu, fWeight, sim::kUniformTime);
    }

    else if (fSamplingMode == TimeGenSamplingMode::UNIFORM_ENERGY)
    {
      // Select a time bin using the energy-integrated spectrum
      TH1D* t_hist = fSpectrumHist->ProjectionX("dummy_time_hist");
      double* time_bin_weights = t_hist->GetArray();

      // Sample a time bin from the distribution
      std::discrete_distribution<int> time_dist;
      std::discrete_distribution<int>::param_type time_params(time_bin_weights,
        &(time_bin_weights[t_hist->GetNbinsX() + 1]));
      int time_bin_index = gen.discrete_sample(time_dist, time_params);

      // Sample a time uniformly from within the selected time bin
      double t_min = t_hist->GetBinLowEdge(time_bin_index);
      double t_max = t_min + t_hist->GetBinWidth(time_bin_index);
      // sample a time on [ t_min, t_max )
      fTNu = gen.uniform_random_double(t_min, t_max, false);

      // Sample an energy uniformly over the entire allowed range
      // underflow bin has index zero
      TAxis* energy_axis = fSpectrumHist->GetYaxis();
      double E_min = energy_axis->GetBinLowEdge(1);
      double E_max = energy_axis->GetBinLowEdge(energy_axis->GetNbins() + 1);
      double E_nu, total_xs;
      int j = 0;
      do {
        // We have to check that the cross section is nonzero for the sampled
        // energy (otherwise we'll generate an unphysical event). However, if the
        // user has given us a histogram that is below threshold, the
        // program could get stuck here endlessly, sampling rejected energy
        // after rejected energy. Just in case, we cap the total number of tries
        // and quit if things don't work out.
        if (j >= MAX_UNIFORM_ENERGY_ITERATIONS) {
          throw cet::exception("MARLEYTimeGen") << "Exceeded the maximum of "
            << MAX_UNIFORM_ENERGY_ITERATIONS << " while attempting to sample"
            << " a neutrino energy uniformly.";
        }
        // Sample an energy uniformly on [ E_min, E_max )
        E_nu = gen.uniform_random_double(E_min, E_max, false);
        total_xs = 0.;
        // Check that at least one defined reaction has a non-zero total
        // cross section at the sampled energy. If this is not the case, try
        // again.
        for (const auto& react : gen.get_reactions()) {
          total_xs += react->total_xs(marley_utils::ELECTRON_NEUTRINO, E_nu);
        }

        ++j;
      } while (total_xs <= 0.);

      // Replace the existing neutrino source with a monoenergetic one at the
      // neutrino energy that was just sampled above
      std::unique_ptr<marley::NeutrinoSource> nu_source
        = std::make_unique<marley::MonoNeutrinoSource>(
        marley_utils::ELECTRON_NEUTRINO, E_nu);
      gen.set_source(std::move(nu_source));

      // Generate a MARLEY event using the new monoenergetic source
      truth = fMarleyGenerator->create_MCTruth(vertex_pos, fEvent.get());

      // Get the value of the true dependent probability density (probability
      // of the sampled energy given the sampled time) to use as a biasing
      // correction in the neutrino vertex weight.
      //
      // Get a 1D projection of the energy spectrum for the sampled time bin
      TH1D* energy_spect = fSpectrumHist->ProjectionY("energy_spect",
        time_bin_index, time_bin_index);

      // Create a new MARLEY neutrino source object using this projection (this
      // will create a normalized probability density that we can use) and load
      // it into the generator.
      nu_source = marley_root::make_root_neutrino_source(
        marley_utils::ELECTRON_NEUTRINO, energy_spect);
      double new_source_E_min = nu_source->get_Emin();
      double new_source_E_max = nu_source->get_Emax();
      gen.set_source(std::move(nu_source));
      // NOTE: The marley::Generator object normalizes the E_pdf to unity
      // automatically, but just in case, we redo it here.
      double E_pdf_integ = marley_utils::num_integrate([&gen](double E_nu)
        -> double { return gen.E_pdf(E_nu); }, new_source_E_min,
        new_source_E_max, 1e4);

      // Compute the likelihood ratio that we need to bias the neutrino vertex
      // weight
      double weight_bias = (gen.E_pdf(E_nu) / E_pdf_integ) / (E_max - E_min);

      fWeight = fFluxAveragedCrossSection * weight_bias;
      
      sn_truth = sim::SupernovaTruth(fTNu, fWeight, sim::kUniformEnergy);
    }

    else {
      throw cet::exception("MARLEYTimeGen") << "Unrecognized sampling mode"
        << " encountered in evgen::MarleyTimeGen::produce()";
    }

    // Write the marley::Event object to the event tree
    fEventTree->Fill();

    // Add the truth objects to the appropriate vectors
    truthcol->push_back(truth);

    sn_truthcol->push_back(sn_truth);

    // Associate the last entries in each of the truth object vectors (the
    // truth objects that were just created for the current neutrino vertex)
    // with each other
    util::CreateAssn(*this, e, *truthcol, *sn_truthcol, *truth_assns,
      truthcol->size() - 1, truthcol->size()/*, sn_truthcol->size() - 1*/);
  }

  // Load the completed truth object vectors and associations into the event
  e.put(std::move(truthcol));

  e.put(std::move(sn_truthcol));

  e.put(std::move(truth_assns));
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::reconfigure(const Parameters& p)
{
  const auto& seed_service = art::ServiceHandle<rndm::NuRandomService>();
  const auto& geom_service = art::ServiceHandle<geo::Geometry>();

  // Create a new evgen::ActiveVolumeVertexSampler object based on the current
  // configuration
  fVertexSampler = std::make_unique<evgen::ActiveVolumeVertexSampler>(
    p().vertex_, *seed_service, *geom_service, "MARLEY_Vertex_Sampler");

  // Create a new marley::Generator object based on the current configuration
  fMarleyGenerator = std::make_unique<MARLEYGenerator>(p().marley_parameters_,
    *seed_service, "MARLEY");

  // Retrieve the time-dependent neutrino flux from a ROOT file
  std::string full_spectrum_file_name
    = fMarleyGenerator->find_file(p().spectrum_file_(), "spectrum");

  std::unique_ptr<TFile> spectrum_file
    = std::make_unique<TFile>(full_spectrum_file_name.c_str(), "read");
  TH2D* temp_h2 = nullptr;
  spectrum_file->GetObject(p().namecycle_().c_str(), temp_h2);
  fSpectrumHist.reset(temp_h2);

  // Disassociate the TH2D from its parent TFile. If we fail to do this,
  // then ROOT will auto-delete the TH2D when the TFile goes out of scope.
  fSpectrumHist->SetDirectory(nullptr);

  // Get the number of neutrino vertices per event from the FHiCL parameters
  fNeutrinosPerEvent = p().nu_per_event_();

  // Determine the current sampling mode from the FHiCL parameters
  const std::string& samp_mode_str = p().sampling_mode_();
  if (samp_mode_str == "histogram")
    fSamplingMode = TimeGenSamplingMode::HISTOGRAM;
  else if (samp_mode_str == "uniform time")
    fSamplingMode = TimeGenSamplingMode::UNIFORM_TIME;
  else if (samp_mode_str == "uniform energy")
    fSamplingMode = TimeGenSamplingMode::UNIFORM_ENERGY;
  else throw cet::exception("MARLEYTimeGen")
      << "Invalid sampling mode \"" << samp_mode_str << "\""
      << " specified for the MARLEYTimeGen module.";

  LOG_INFO("MARLEYTimeGen") << fNeutrinosPerEvent << " neutrino vertices"
    << " will be generated for each art::Event using the \"" << samp_mode_str
    << "\" sampling mode.";
 
  // Compute the flux-averaged total cross section using MARLEY. This will be
  // used to compute neutrino vertex weights for the sim::SupernovaTruth
  // objects.
  marley::Generator& gen = fMarleyGenerator->get_generator();

  // Get a 1D projection of the energy spectrum (integrated over time)
  TH1D* energy_spect = fSpectrumHist->ProjectionY("energy_spect");

  // Create a new MARLEY neutrino source object using this projection
  // TODO: replace the hard-coded electron neutrino PDG code here (and in
  // several other places in this source file) when you're ready to use
  // MARLEY with multiple neutrino flavors
  std::unique_ptr<marley::NeutrinoSource> nu_source
    = marley_root::make_root_neutrino_source(marley_utils::ELECTRON_NEUTRINO,
    energy_spect);

  // Get an integral of the source PDF (in case it isn't normalized to 1)
  double source_integ = marley_utils::num_integrate([&nu_source](double E_nu)
    -> double { return nu_source->pdf(E_nu); }, nu_source->get_Emin(),
    nu_source->get_Emax(), 1e4);

  // Factor of hbar_c^2 converts from MeV^(-2) to fm^2
  fFluxAveragedCrossSection = marley_utils::hbar_c2 / source_integ
    * marley_utils::num_integrate([&nu_source, &gen](double E_nu) -> double
    {
      double xs = 0.;
      for (const auto& react : gen.get_reactions()) {
        xs += react->total_xs(marley_utils::ELECTRON_NEUTRINO, E_nu);
      }
      return xs * nu_source->pdf(E_nu);
    }, nu_source->get_Emin(), nu_source->get_Emax(), 1e4);

  LOG_INFO("MARLEYTimeGen") << "The flux-averaged total cross section"
    << " predicted by MARLEY for the current supernova spectrum is "
    << fFluxAveragedCrossSection << " fm^2";

  // For speed, sample energies first whenever possible (and then sample from
  // an energy-dependent timing distribution). This avoids unnecessary calls to
  // MARLEY to change the energy spectrum.
  if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM ||
    fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
  {
    gen.set_source(std::move(nu_source));
  }
}

DEFINE_ART_MODULE(evgen::MarleyTimeGen)
