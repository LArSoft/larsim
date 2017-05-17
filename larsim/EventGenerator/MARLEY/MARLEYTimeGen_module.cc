//////////////////////////////////////////////////////////////////////////////
/// \file MARLEYTimeGen_module.cc
/// \brief Module that allows for sampling neutrino energies and times from
/// time-dependent supernova spectra. This module uses MARLEY to help generate
/// events.
///
/// \author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

// standard library includes
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
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/EventGenerator/MARLEY/MARLEYGenerator.h"
#include "larsim/EventGenerator/MARLEY/ActiveVolumeVertexSampler.h"

// ROOT includes
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TVectorD.h"

// MARLEY includes
#include "marley/marley_root.hh"
#include "marley/marley_utils.hh"

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
};

//------------------------------------------------------------------------------
evgen::MarleyTimeGen::MarleyTimeGen(const Parameters& p)
  : fEvent(new marley::Event), fRunNumber(0), fSubRunNumber(0), fEventNumber(0)
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

  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();
}

//------------------------------------------------------------------------------
evgen::MarleyTimeGen::~MarleyTimeGen()
{
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;
  std::unique_ptr<sumdata::RunData>
    runcol(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runcol));
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::produce(art::Event& e)
{
  // Get the run, subrun, and event numbers from the current art::Event
  fRunNumber = e.run();
  fSubRunNumber = e.subRun();
  fEventNumber = e.event();

  std::unique_ptr< std::vector<simb::MCTruth> >
    truthcol(new std::vector<simb::MCTruth>);

  // Get the primary vertex location for this event
  art::ServiceHandle<geo::Geometry> geo;
  TLorentzVector vertex_pos = fVertexSampler->sample_vertex_pos(*geo);

  fTNu = 0.;

  marley::Generator& gen = fMarleyGenerator->get_generator();

  simb::MCTruth truth;

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
  }

  else if (fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
  {
    // Generate a MARLEY event using the time-integrated spectrum
    // (the generator was already configured to use it by reconfigure())
    truth = fMarleyGenerator->create_MCTruth(vertex_pos,
      fEvent.get());

    // Sample a time uniformly
    TAxis* time_axis = fSpectrumHist->GetXaxis();
    double t_min = time_axis->GetBinLowEdge(1); // underflow bin has index zero
    double t_max = time_axis->GetBinLowEdge(time_axis->GetNbins() + 1);
    // sample a time on [ t_min, t_max )
    fTNu = gen.uniform_random_double(t_min, t_max, false);
  }

  if (fSamplingMode == TimeGenSamplingMode::UNIFORM_ENERGY)
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
    // sample an energy uniformly on [ E_min, E_max )
    double E_nu = gen.uniform_random_double(E_min, E_max, false);

    // TODO: implement event creation using uniformly sampled energies
    throw cet::exception("MARLEYTimeGen")
      << "The \"uniform energy\" sampling mode is not fully implemented yet.";
  }

  // Write the marley::Event object to the event tree
  fEventTree->Fill();

  truthcol->push_back(truth);

  e.put(std::move(truthcol));
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
    = fMarleyGenerator->find_file("nusperbin2d.root", "spectrum");

  std::unique_ptr<TFile> spectrum_file
    = std::make_unique<TFile>(full_spectrum_file_name.c_str(), "read");
  TH2D* temp_h2 = nullptr;
  spectrum_file->GetObject("nusperbin2d_nue", temp_h2);
  fSpectrumHist.reset(temp_h2);

  // Disassociate the TH2D from its parent TFile. If we fail to do this,
  // then ROOT will auto-delete the TH2D when the TFile goes out of scope.
  fSpectrumHist->SetDirectory(nullptr);

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

  // For speed, sample energies first whenever possible (and then sample from
  // an energy-dependent timing distribution). This avoids unnecessary calls to
  // MARLEY to change the energy spectrum.
  if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM ||
    fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
  {
    marley::Generator& gen = fMarleyGenerator->get_generator();

    // Get a 1D projection of the energy spectrum (integrated over time)
    TH1D* energy_spect = fSpectrumHist->ProjectionY("energy_spect");

    // Replace the existing neutrino source used by MARLEY with this
    // projection
    // TODO: replace hard-coded electron neutrino PDG code here when
    // you're ready to use MARLEY with multiple neutrino flavors
    std::unique_ptr<marley::NeutrinoSource> nu_source
      = marley_root::make_root_neutrino_source(marley_utils::ELECTRON_NEUTRINO,
      energy_spect);

    gen.set_source(std::move(nu_source));
  }
}

DEFINE_ART_MODULE(evgen::MarleyTimeGen)
