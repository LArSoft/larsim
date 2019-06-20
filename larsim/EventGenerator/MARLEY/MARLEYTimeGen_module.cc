//////////////////////////////////////////////////////////////////////////////
/// @file MARLEYTimeGen_module.cc
/// @brief Module that allows for sampling neutrino energies and times from
/// time-dependent supernova spectra. This module uses MARLEY to help generate
/// events.
///
/// @author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

// standard library includes
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <regex>
#include <string>

// framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Simulation/SupernovaTruth.h"
#include "larsim/EventGenerator/MARLEY/MARLEYHelper.h"
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
#include "marley/Integrator.hh"

// Anonymous namespace for definitions local to this source file
namespace {

  // The number of times to attempt to sample an energy uniformly
  // before giving up
  constexpr int MAX_UNIFORM_ENERGY_ITERATIONS = 1000;

  // Neutrino vertices generated using unbiased sampling are assigned
  // unit weight
  constexpr double ONE = 1.;

  // Number of sampling points to use for numerical integration
  // (via the Clenshaw-Curtis method)
  constexpr int NUM_INTEGRATION_SAMPLING_POINTS = 100;

  // Multiply an energy in MeV by this factor to convert to erg
  // (based on 2017 PDG "Physical Constants" article)
  constexpr double MeV_to_erg = 1.6021766208e-6; // erg / MeV

  // The PDG codes that represent one of the varieties of "nu_x"
  constexpr std::array<int, 4> NU_X_PDG_CODES
  {
    marley_utils::MUON_NEUTRINO, marley_utils::MUON_ANTINEUTRINO,
    marley_utils::TAU_NEUTRINO, marley_utils::TAU_ANTINEUTRINO,
  };

  // Helper function that tests whether a given PDG code represents
  // a "nu_x"
  bool is_nux(int pdg_code) {
    return std::find(std::begin(NU_X_PDG_CODES),
      std::end(NU_X_PDG_CODES), pdg_code) != std::end(NU_X_PDG_CODES);
  }

  // Matches comment lines and empty lines in a "fit"-format spectrum file
  const std::regex rx_comment_or_empty = std::regex("\\s*(#.*)?");

  // Compute the flux-averaged total cross section (MeV^[-2]) for all
  // defined reactions and for a particular neutrino source
  double flux_averaged_total_xs(marley::NeutrinoSource& nu_source,
    marley::Generator& gen);

  // Overloaded version that allows access to the flux integral
  // (which will be loaded into source_integ) and the xs * flux integral
  // (which will be loaded into tot_xs_integ)
  double flux_averaged_total_xs(marley::NeutrinoSource& nu_source,
    marley::Generator& gen, double& source_integ, double& tot_xs_integ);

  // Returns a numerical integral of the function f on the interval
  // [x_min, x_max]
  double integrate(const std::function<double(double)>& f, double x_min,
    double x_max);
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

      fhicl::Table<evgen::MARLEYHelper::Config> marley_parameters_ {
        Name("marley_parameters"),
        Comment("Configuration for the MARLEY generator. Note that for"
          " MARLEYTimeGen, the source configuration given here is ignored.")
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

      fhicl::Atom<std::string> spectrum_file_format_ {
        Name("spectrum_file_format"),
        Comment("Format to assume for the neutrino spectrum file."
          " Valid options are \"th2d\" (a ROOT file containing a "
          " TH2D object) and \"fit\" (an ASCII text file containing"
          " fit parameters for each time bin). This parameter is not"
          " case-sensitive."),
        "th2d" // default value
      };

      fhicl::Atom<std::string> spectrum_file_ {
        Name("spectrum_file"),
        Comment("Name of a file that contains a representation"
          " of the incident (not cross section weighted) neutrino spectrum"
          " as a function of time and energy.")
      };

      fhicl::OptionalAtom<std::string> pinching_parameter_type_ {
        Name("pinching_parameter_type"),
        Comment("Type of pinching parameter to assume when parsing"
          " the time-dependent fit parameters for the incident neutrino"
          " spectrum. Valid options are \"alpha\" and \"beta\". This"
          " parameter is not case-sensitive."),
        [this]() -> bool {
          auto spectrum_file_format
            = marley_utils::to_lowercase(spectrum_file_format_());
          return (spectrum_file_format == "fit");
        }
      };

      fhicl::OptionalAtom<std::string> namecycle_ {
        Name("namecycle"),
        Comment("Name of the TH2D object to use to represent the"
          " incident neutrino spectrum. This value should match the"
          " name of the TH2D as given in the ROOT file specified"
          " in the \"spectrum_file\" parameter. The TH2D should use "
          " time bins on the X axis (seconds) and energy bins on the "
          " Y axis (MeV)."),
        [this]() -> bool {
          auto spectrum_file_format
            = marley_utils::to_lowercase(spectrum_file_format_());
          return (spectrum_file_format == "th2d");
        }
      };

      fhicl::OptionalAtom<double> fit_Emin_ {
        Name("fit_Emin"),
        Comment("Minimum allowed neutrino energy (MeV) for a \"fit\" format"
          " spectrum file"),
        [this]() -> bool {
          auto spectrum_file_format
            = marley_utils::to_lowercase(spectrum_file_format_());
          return (spectrum_file_format == "fit");
        }
      };

      fhicl::OptionalAtom<double> fit_Emax_ {
        Name("fit_Emax"),
        Comment("Maximum allowed neutrino energy (MeV) for a \"fit\" format"
          " spectrum file"),
        [this]() -> bool {
          auto spectrum_file_format
            = marley_utils::to_lowercase(spectrum_file_format_());
          return (spectrum_file_format == "fit");
        }
      };

    }; // struct Config

    // Type to enable FHiCL parameter validation by art
    using Parameters = art::EDProducer::Table<Config>;

    // @brief Configuration-checking constructor
    explicit MarleyTimeGen(const Parameters& p);

    virtual void produce(art::Event& e) override;
    virtual void beginRun(art::Run& run) override;

    virtual void reconfigure(const Parameters& p);

  protected:

    /// @brief Stores parsed fit parameters from a single time bin and neutrino
    /// type in a "fit"-format spectrum file
    class FitParameters {
      public:

        FitParameters(double Emean, double alpha, double luminosity)
          : fEmean(Emean), fAlpha(alpha), fLuminosity(luminosity) {}

        /// @brief Mean neutrino energy (MeV)
        double Emean() const { return fEmean; }
        /// @brief Pinching parameter (dimensionless)
        double Alpha() const { return fAlpha; }
        /// @brief Luminosity (erg / s)
        double Luminosity() const { return fLuminosity; }

        /// @brief Set the mean neutrino energy (MeV)
        void set_Emean(double Emean) { fEmean = Emean; }
        /// @brief Set the pinching parameter
        void set_Alpha(double alpha) { fAlpha = alpha; }
        /// @brief Set the luminosity (erg / s)
        void set_Luminosity(double lum) { fLuminosity = lum; }

        /// @brief Converts an iterator that points to a FitParameters
        /// object into an iterator that points to that object's
        /// fLuminosity member.
        /// @details This function helps us to be able to sample time bins
        /// with a std::discrete_distribution using the bin luminosities
        /// without redundnant storage.
        template<typename It> static marley::IteratorToMember<It,
          FitParameters, double> make_luminosity_iterator(It it)
        {
          return marley::IteratorToMember<It, FitParameters,
            double>(it, &FitParameters::fLuminosity);
        }

      protected:
        double fEmean; ///< Mean neutrino energy (MeV)
        double fAlpha; ///< Pinching parameter
        double fLuminosity; ///< Luminosity (erg / s)
    };

    /// @brief Stores fitting parameters for all neutrino types from a single
    /// time bin in a "fit"-format spectrum file
    class TimeFit {
      public:
        TimeFit(double time, double Emean_nue, double alpha_nue,
          double lum_nue, double Emean_nuebar, double alpha_nuebar,
          double lum_nuebar, double Emean_nux, double alpha_nux,
          double lum_nux) : fTime(time),
          fNueFitParams(Emean_nue, alpha_nue, lum_nue),
          fNuebarFitParams(Emean_nuebar, alpha_nuebar, lum_nuebar),
          fNuxFitParams(Emean_nux, alpha_nux, lum_nux) {}

        /// @brief Returns the low edge (in s) of the time bin that uses
        /// these fit parameters
        /// @details Time zero is defined to be the core bounce
        double Time() const { return fTime; }

      protected:
        /// @brief Helper function that returns a pointer-to-member
        /// for the FitParameters object appropriate for a given
        /// neutrino type
        /// @param pdg_code PDG code for the neutrino type of interest
        static FitParameters TimeFit::* GetFitParametersMemberPointer(
          int pdg_code)
        {
          if (pdg_code == marley_utils::ELECTRON_NEUTRINO) {
            return &TimeFit::fNueFitParams;
          }
          else if (pdg_code == marley_utils::ELECTRON_ANTINEUTRINO) {
            return &TimeFit::fNuebarFitParams;
          }
          else if ( is_nux(pdg_code) )
          {
            // The PDG code represents one of the varieties of nu_x
            return &TimeFit::fNuxFitParams;
          }
          else throw cet::exception("MARLEYTimeGen") << "Invalid neutrino"
            << " PDG code " << pdg_code << " encountered in MARLEYTimeGen"
            << "::TimeFit::GetFitParametersMemberPointer()";
          return nullptr;
        }

      public:

        /// @brief Retrieves fit parameters for a specific neutrino type
        /// for this time bin
        /// @param pdg_code The PDG code for the desired neutrino type
        const FitParameters& GetFitParameters(int pdg_code) const {
          return this->*GetFitParametersMemberPointer(pdg_code);
        }

        /// @brief Converts an iterator that points to a TimeFit
        /// object into an iterator that points to the appropriate
        /// (based on the neutrino PDG code) FitParameters member
        /// owned by the TimeFit object.
        /// @details This function helps us to be able to sample time bins
        /// with a std::discrete_distribution using the bin luminosities
        /// without redundnant storage.
        /// @param pdg_code PDG code for the neutrino type of interest
        /// @param iterator An iterator to a TimeFit object that will be
        /// converted
        template<typename It> static marley::IteratorToMember<It,
          TimeFit, FitParameters> make_FitParameters_iterator(int pdg_code,
          It iterator)
        {
          return marley::IteratorToMember<It, TimeFit, FitParameters>(
            iterator, GetFitParametersMemberPointer(pdg_code) );
        }

      protected:
        /// @brief Time bin low edge (s)
        double fTime;

        /// @brief Fitting parameters for electron neutrinos in this time bin
        FitParameters fNueFitParams;

        /// @brief Fitting parameters for electron antineutrinos in this time
        /// bin
        FitParameters fNuebarFitParams;

        /// @brief Fitting parameters for non-electron-flavor neutrinos in this
        /// time bin
        FitParameters fNuxFitParams;
    };

    /// @brief Creates a simb::MCTruth object using a uniformly-sampled
    /// neutrino energy
    /// @details This function samples a reacting neutrino energy uniformly
    /// from the full range of energies allowed by the incident spectrum and
    /// the currently defined reactions.
    simb::MCTruth make_uniform_energy_mctruth(double E_min, double E_max,
      double& E_nu, const TLorentzVector& vertex_pos);

    /// @brief Create a MARLEY neutrino source object using a set of fit
    /// parameters for a particular time bin
    std::unique_ptr<marley::NeutrinoSource> source_from_time_fit(
      const TimeFit& fit);

    /// @brief Create simb::MCTruth and sim::SupernovaTruth objects using
    /// spectrum information from a ROOT TH2D
    void create_truths_th2d(simb::MCTruth& mc_truth,
      sim::SupernovaTruth& sn_truth, const TLorentzVector& vertex_pos);

    /// @brief Create simb::MCTruth and sim::SupernovaTruth objects using a
    /// neutrino spectrum described by a previously-parsed "fit"-format file
    void create_truths_time_fit(simb::MCTruth& mc_truth,
      sim::SupernovaTruth& sn_truth, const TLorentzVector& vertex_pos);

    /// @brief Helper function that makes a final dummy TimeFit object so that
    /// the final real time bin can have a right edge
    void make_final_timefit(double time);

    /// @brief Makes ROOT histograms showing the emitted neutrinos in each time
    /// bin when using a "fit"-format spectrum file
    void make_nu_emission_histograms() const;

    /// @brief Object that provides an interface to the MARLEY event generator
    std::unique_ptr<evgen::MARLEYHelper> fMarleyHelper;

    /// @brief Algorithm that allows us to sample vertex locations within the
    /// active volume(s) of the detector
    std::unique_ptr<evgen::ActiveVolumeVertexSampler> fVertexSampler;

    /// @brief unique_ptr to the current event created by MARLEY
    std::unique_ptr<marley::Event> fEvent;

    /// @brief ROOT TH2D that contains the time-dependent spectrum to use when
    /// sampling neutrino times and energies.
    /// @details This member is only used when reading the spectrum from a ROOT
    /// file.
    std::unique_ptr<TH2D> fSpectrumHist;

    /// @brief Vector that contains the fit parameter information for each time
    /// bin when using a "fit"-format spectrum file.
    /// @details This member is unused when the spectrum is read from a ROOT
    /// file.
    std::vector<TimeFit> fTimeFits;

    /// @enum TimeGenSamplingMode
    /// @brief Enumerated type that defines the allowed ways that a neutrino's
    /// energy and arrival time may be sampled
    ///
    /// @var TimeGenSamplingMode::HISTOGRAM
    /// @brief Sample from the input spectrum histogram directly (after cross
    /// section weighting), without any biasing
    ///
    /// @var TimeGenSamplingMode::UNIFORM_TIME
    /// @brief Sample energies without biasing, sample times uniformly
    /// over the entire allowed range
    ///
    /// @var TimeGenSamplingMode::UNIFORM_ENERGY
    /// @brief Sample energies uniformly over the entire allowed range,
    /// sample times without biasing
    enum class TimeGenSamplingMode { HISTOGRAM, UNIFORM_TIME, UNIFORM_ENERGY };

    /// @brief Represents the sampling mode to use when selecting neutrino
    /// times and energies
    TimeGenSamplingMode fSamplingMode;

    /// @enum PinchParamType
    /// @brief Enumerated type that defines the pinching parameter conventions
    /// that are understood by this module
    ///
    /// @var PinchParamType::ALPHA
    /// @brief Use the fitting convention in which the pinching parameter
    /// is called @f$\alpha@f$
    /// @details In the "alpha" convention, the probability density
    /// function @f$p(E_\nu)@f$ describing the distribution of neutrino
    /// energies @f$E_\nu@f$ is given by
    /// @f$p(E_\nu) = @f$\left(\frac{\alpha + 1}{\left<E_\nu\right>}
    /// \right)^{\alpha + 1}\frac{E_\nu^\alpha}{\Gamma(\alpha + 1)}
    /// \exp\left(\frac{-[\alpha + 1]E_\nu}{\left<E_\nu\right>}\right)@f$
    ///
    /// @var PinchParamType::BETA
    /// @brief Use the fitting convention in which the pinching parameter
    /// is called @f$\beta@f$
    /// @details In the "beta" convention, the probability density
    /// function @f$p(E_\nu)@f$ describing the distribution of neutrino
    /// energies @f$E_\nu@f$ is given by
    /// @f$p(E_\nu) = @f$\left(\frac{\beta}{\left<E_\nu\right>}
    /// \right)^{\beta}\frac{E_\nu^(\beta - 1)}{\Gamma(\beta)}
    /// \exp\left(\frac{-\beta E_\nu}{\left<E_\nu\right>}\right)@f$
    /// Note that this differs from the "alpha" convention via
    /// @f$\beta = \alpha + 1@f$.
    enum class PinchParamType { ALPHA, BETA };

    /// @brief The pinching parameter convention to use when interpreting the
    /// time-dependent fits
    PinchParamType fPinchType;

    /// @enum SpectrumFileFormat
    /// @brief Enumerated type that defines the allowed neutrino spectrum input
    /// file formats
    ///
    /// @var SpectrumFileFormat::RootTH2D
    /// @brief The incident neutrino spectrum is described by a TH2D object
    /// stored in a ROOT file
    /// @details The X axis of the TH2D should be the time axis (s) and the
    /// Y axis should be the energy axis (MeV).
    // TODO: add information about the bin value units
    ///
    /// @var SpectrumFileFormat::FIT
    /// @brief The incident neutrino spectrum is described by a set of fitting
    /// parameters for each time bin
    // TODO: add a "fit"-format description here
    enum class SpectrumFileFormat { RootTH2D, FIT };

    /// @brief Format to assume for the neutrino spectrum input file
    SpectrumFileFormat fSpectrumFileFormat;

    /// @brief The event TTree created by MARLEY
    /// @details This tree will be saved to the "hist" output file for
    /// validation purposes. The tree contains the same information as the
    /// generated simb::MCTruth objects, but in MARLEY's internal format
    TTree* fEventTree;

    /// @brief Run number from the art::Event being processed
    uint_fast32_t fRunNumber;
    /// @brief Subrun number from the art::Event being processed
    uint_fast32_t fSubRunNumber;
    /// @brief Event number from the art::Event being processed
    uint_fast32_t fEventNumber;

    /// @brief Time since the supernova core bounce for the current MARLEY
    /// neutrino vertex
    double fTNu;

    /// @brief Statistical weight for the current MARLEY neutrino vertex
    double fWeight;

    /// @brief Flux-averaged total cross section (fm<sup>2</sup>, average is
    /// taken over all energies and times for all defined reactions) used by
    /// MARLEY to generate neutrino vertices
    double fFluxAveragedCrossSection; // fm^2

    /// @brief The number of MARLEY neutrino vertices to generate in each
    /// art::Event
    unsigned int fNeutrinosPerEvent;

    /// @brief Minimum neutrino energy to consider when using a "fit"-format
    /// spectrum file
    double fFitEmin;

    /// @brief Maximum neutrino energy to consider when using a "fit"-format
    /// spectrum file
    double fFitEmax;
};

//------------------------------------------------------------------------------
evgen::MarleyTimeGen::MarleyTimeGen(const Parameters& p)
  : EDProducer{p}, fEvent(new marley::Event), fRunNumber(0), fSubRunNumber(0),
  fEventNumber(0), fTNu(0.), fFluxAveragedCrossSection(0.)
{
  // Configure the module (including MARLEY itself) using the FHiCL parameters
  this->reconfigure(p);

  // Create a ROOT TTree using the TFileService that will store the MARLEY
  // event objects (useful for debugging purposes)
  art::ServiceHandle<art::TFileService const> tfs;
  fEventTree = tfs->make<TTree>("MARLEY_event_tree",
    "Neutrino events generated by MARLEY");
  fEventTree->Branch("event", "marley::Event", fEvent.get());

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
void evgen::MarleyTimeGen::beginRun(art::Run& run) {
  art::ServiceHandle<geo::Geometry const> geo;
  run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::create_truths_th2d(simb::MCTruth& mc_truth,
  sim::SupernovaTruth& sn_truth, const TLorentzVector& vertex_pos)
{
  // Get a reference to the generator object created by MARLEY (we'll need
  // to do a few fancy things with it other than just creating events)
  marley::Generator& gen = fMarleyHelper->get_generator();

  if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM)
  {
    // Generate a MARLEY event using the time-integrated spectrum
    // (the generator was already configured to use it by reconfigure())
    mc_truth = fMarleyHelper->create_MCTruth(vertex_pos,
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
    int time_bin_index = gen.sample_from_distribution(time_dist, time_params);

    // Sample a time uniformly from within the selected time bin
    double t_min = t_hist->GetBinLowEdge(time_bin_index);
    double t_max = t_min + t_hist->GetBinWidth(time_bin_index);
    // sample a time on [ t_min, t_max )
    fTNu = gen.uniform_random_double(t_min, t_max, false);
    // Unbiased sampling was used, so assign this neutrino vertex a
    // unit statistical weight
    fWeight = ONE;

    sn_truth = sim::SupernovaTruth(fTNu, fWeight, fFluxAveragedCrossSection,
      sim::kUnbiased);
  }

  else if (fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
  {
    // Generate a MARLEY event using the time-integrated spectrum
    // (the generator was already configured to use it by reconfigure())
    mc_truth = fMarleyHelper->create_MCTruth(vertex_pos,
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
    double weight_bias = t_hist->GetBinContent(t_bin_index) * (t_max - t_min)
      / ( t_hist->Integral() * t_hist->GetBinWidth(t_bin_index) );

    fWeight = weight_bias;

    sn_truth = sim::SupernovaTruth(fTNu, fWeight, fFluxAveragedCrossSection,
      sim::kUniformTime);
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
    int time_bin_index = gen.sample_from_distribution(time_dist, time_params);

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

    double E_nu = std::numeric_limits<double>::lowest();

    // Generate a MARLEY event using a uniformly sampled energy
    mc_truth = make_uniform_energy_mctruth(E_min, E_max, E_nu, vertex_pos);

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
    auto nu_source = marley_root::make_root_neutrino_source(
      marley_utils::ELECTRON_NEUTRINO, energy_spect);
    double new_source_E_min = nu_source->get_Emin();
    double new_source_E_max = nu_source->get_Emax();
    gen.set_source(std::move(nu_source));
    // NOTE: The marley::Generator object normalizes the E_pdf to unity
    // automatically, but just in case, we redo it here.
    double E_pdf_integ = integrate([&gen](double E_nu)
      -> double { return gen.E_pdf(E_nu); }, new_source_E_min,
      new_source_E_max);

    // Compute the likelihood ratio that we need to bias the neutrino vertex
    // weight
    double weight_bias = (gen.E_pdf(E_nu) / E_pdf_integ) * (E_max - E_min);

    fWeight = weight_bias;

    sn_truth = sim::SupernovaTruth(fTNu, fWeight, fFluxAveragedCrossSection,
      sim::kUniformEnergy);
  }

  else {
    throw cet::exception("MARLEYTimeGen") << "Unrecognized sampling mode"
      << " encountered in evgen::MarleyTimeGen::produce()";
  }

}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::produce(art::Event& e)
{
  art::ServiceHandle<geo::Geometry const> geo;

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

  for (unsigned int n = 0; n < fNeutrinosPerEvent; ++n) {

    // Sample a primary vertex location for this event
    TLorentzVector vertex_pos = fVertexSampler->sample_vertex_pos(*geo);

    // Reset the neutrino's time-since-supernova to a bogus value (for now)
    fTNu = std::numeric_limits<double>::lowest();

    if (fSpectrumFileFormat == SpectrumFileFormat::RootTH2D) {
      create_truths_th2d(truth, sn_truth, vertex_pos);
    }
    else if (fSpectrumFileFormat == SpectrumFileFormat::FIT) {
      create_truths_time_fit(truth, sn_truth, vertex_pos);
    }
    else {
      throw cet::exception("MARLEYTimeGen") << "Invalid spectrum file"
        << " format encountered in evgen::MarleyTimeGen::produce()";
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
  const auto& geom_service = art::ServiceHandle<geo::Geometry const>();

  // Create a new evgen::ActiveVolumeVertexSampler object based on the current
  // configuration
  fVertexSampler = std::make_unique<evgen::ActiveVolumeVertexSampler>(
    p().vertex_, *seed_service, *geom_service, "MARLEY_Vertex_Sampler");

  // Create a new marley::Generator object based on the current configuration
  fMarleyHelper = std::make_unique<MARLEYHelper>(p().marley_parameters_,
    *seed_service, "MARLEY");

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

  MF_LOG_INFO("MARLEYTimeGen") << fNeutrinosPerEvent << " neutrino vertices"
    << " will be generated for each art::Event using the \"" << samp_mode_str
    << "\" sampling mode.";

  // Retrieve the time-dependent neutrino spectrum from the spectrum file.
  // Use different methods depending on the file's format.
  std::string spectrum_file_format = marley_utils::to_lowercase(
    p().spectrum_file_format_());

  if (spectrum_file_format == "th2d")
    fSpectrumFileFormat = SpectrumFileFormat::RootTH2D;
  else if (spectrum_file_format == "fit") {
    fSpectrumFileFormat = SpectrumFileFormat::FIT;

    std::string pinch_type;
    if ( !p().pinching_parameter_type_(pinch_type) ) {
      throw cet::exception("MARLEYTimeGen") << "Missing pinching parameter"
        << " type for a \"fit\" format spectrum file";
    }

    marley_utils::to_lowercase_inplace(pinch_type);
    if (pinch_type == "alpha") fPinchType = PinchParamType::ALPHA;
    else if (pinch_type == "beta") fPinchType = PinchParamType::BETA;
    else throw cet::exception("MARLEYTimeGen")
      << "Invalid pinching parameter type \"" << pinch_type
      << "\" specified for the MARLEYTimeGen module.";

    if ( !p().fit_Emin_(fFitEmin) ) throw cet::exception("MARLEYTimeGen")
      << "Missing minimum energy for a \"fit\" format spectrum"
      << " used by the MARLEYTimeGen module.";

    if ( !p().fit_Emax_(fFitEmax) ) throw cet::exception("MARLEYTimeGen")
      << "Missing maximum energy for a \"fit\" format spectrum"
      << " used by the MARLEYTimeGen module.";

    if (fFitEmax < fFitEmin) throw cet::exception("MARLEYTimeGen")
      << "Maximum energy is less than the minimum energy for"
      << " a \"fit\" format spectrum used by the MARLEYTimeGen module.";
  }
  else throw cet::exception("MARLEYTimeGen")
      << "Invalid spectrum file format \"" << p().spectrum_file_format_()
      << "\" specified for the MARLEYTimeGen module.";

  // Determine the full file name (including path) of the spectrum file
  std::string full_spectrum_file_name
    = fMarleyHelper->find_file(p().spectrum_file_(), "spectrum");

  marley::Generator& gen = fMarleyHelper->get_generator();

  if (fSpectrumFileFormat == SpectrumFileFormat::RootTH2D) {

    // Retrieve the time-dependent neutrino flux from a ROOT file
    std::unique_ptr<TFile> spectrum_file
      = std::make_unique<TFile>(full_spectrum_file_name.c_str(), "read");
    TH2D* temp_h2 = nullptr;
    std::string namecycle;
    if ( !p().namecycle_(namecycle) ) {
      throw cet::exception("MARLEYTimeGen") << "Missing namecycle for"
        << " a TH2D spectrum file";
    }

    spectrum_file->GetObject(namecycle.c_str(), temp_h2);
    fSpectrumHist.reset(temp_h2);

    // Disassociate the TH2D from its parent TFile. If we fail to do this,
    // then ROOT will auto-delete the TH2D when the TFile goes out of scope.
    fSpectrumHist->SetDirectory(nullptr);

    // Compute the flux-averaged total cross section using MARLEY. This will be
    // used to compute neutrino vertex weights for the sim::SupernovaTruth
    // objects.

    // Get a 1D projection of the energy spectrum (integrated over time)
    TH1D* energy_spect = fSpectrumHist->ProjectionY("energy_spect");

    // Create a new MARLEY neutrino source object using this projection
    // TODO: replace the hard-coded electron neutrino PDG code here (and in
    // several other places in this source file) when you're ready to use
    // MARLEY with multiple neutrino flavors
    std::unique_ptr<marley::NeutrinoSource> nu_source
      = marley_root::make_root_neutrino_source(marley_utils::ELECTRON_NEUTRINO,
      energy_spect);

    // Factor of hbar_c^2 converts from MeV^(-2) to fm^2
    fFluxAveragedCrossSection = marley_utils::hbar_c2
      * flux_averaged_total_xs(*nu_source, gen);

    // For speed, sample energies first whenever possible (and then sample from
    // an energy-dependent timing distribution). This avoids unnecessary calls
    // to MARLEY to change the energy spectrum.
    if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM ||
      fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
    {
      gen.set_source(std::move(nu_source));
    }

  } // spectrum_file_format == "th2d"

  else if (fSpectrumFileFormat == SpectrumFileFormat::FIT) {

    // Clear out the old parameterized spectrum, if one exists
    fTimeFits.clear();

    std::ifstream fit_file(full_spectrum_file_name);
    std::string line;

    bool found_end = false;

    // current line number
    int line_num = 0;
    // number of lines checked in last call to marley_utils::get_next_line()
    int lines_checked = 0;

    double old_time = std::numeric_limits<double>::lowest();

    while ( line = marley_utils::get_next_line(fit_file, rx_comment_or_empty,
      false, lines_checked), line_num += lines_checked, !line.empty() )
    {
      if (found_end) {
        MF_LOG_WARNING("MARLEYTimeGen") << "Trailing content after last time"
          << " bin found on line " << line_num << " of the spectrum file "
          << full_spectrum_file_name;
      }

      double time;
      double Emean_nue, alpha_nue, luminosity_nue;
      double Emean_nuebar, alpha_nuebar, luminosity_nuebar;
      double Emean_nux, alpha_nux, luminosity_nux;
      std::istringstream iss(line);
      bool ok_first = static_cast<bool>( iss >> time );

      if (time <= old_time) throw cet::exception("MARLEYTimeGen")
        << "Time bin left edges given in the spectrum file must be"
        << " strictly increasing. Invalid time bin value found on line "
        << line_num << " of the spectrum file " << full_spectrum_file_name;
      else old_time = time;

      bool ok_rest = static_cast<bool>( iss >> Emean_nue >> alpha_nue
        >> luminosity_nue >> Emean_nuebar >> alpha_nuebar >> luminosity_nuebar
        >> Emean_nux >> alpha_nux >> luminosity_nux
      );

      if (ok_first) {
        // We haven't reached the final bin, so add another time bin
        // in the typical way.
        if (ok_rest) {
          fTimeFits.emplace_back(time, Emean_nue, alpha_nue, luminosity_nue,
            Emean_nuebar, alpha_nuebar, luminosity_nuebar, Emean_nux,
            alpha_nux, luminosity_nux);
        }
        else {
          make_final_timefit(time);
          found_end = true;
        }
      }
      else throw cet::exception("MARLEYTimeGen") << "Parse error on line "
          << line_num << " of the spectrum file " << full_spectrum_file_name;
    }

    if (!found_end) {

      size_t num_time_bins = fTimeFits.size();
      if (num_time_bins < 2) throw cet::exception("MARLEYTimeGen")
        << "Missing right edge for the final time bin in the spectrum file "
        << full_spectrum_file_name << ". Unable to guess a bin width for the "
        << " final bin.";

      // Guess that the width of the penultimate time bin and the width of
      // the final time bin are the same
      double delta_t_bin = fTimeFits.back().Time()
        - fTimeFits.at(num_time_bins - 2).Time();

      double last_bin_right_edge = fTimeFits.back().Time() + delta_t_bin;

      make_final_timefit(last_bin_right_edge);

      MF_LOG_WARNING("MARLEYTimeGen") << "Missing right"
      << " edge for the final time bin in the spectrum file "
      << full_spectrum_file_name << ". Assuming a width of "
      << delta_t_bin << " s for the final bin.";
    }


    // Compute the flux-averaged total cross section for the fitted spectrum.
    // We will need this to compute neutrino vertex weights.
    std::vector< std::unique_ptr< marley::NeutrinoSource > > fit_sources;
    for (const auto& fit : fTimeFits) {
      fit_sources.emplace_back( source_from_time_fit(fit) );
    }

    // TODO: add handling for non-nu_e neutrino types when suitable data become
    // available in MARLEY
    auto temp_source = std::make_unique<marley::FunctionNeutrinoSource>(
      marley_utils::ELECTRON_NEUTRINO, fFitEmin, fFitEmax,
      [&fit_sources, this](double E_nu) -> double {
        double flux = 0.;
        for (size_t s = 0; s < fit_sources.size(); ++s) {

          const TimeFit& current_fit = this->fTimeFits.at(s);
          const auto& current_source = fit_sources.at(s);
          const FitParameters& params = current_fit.GetFitParameters(
            current_source->get_pid() );

          double lum = params.Luminosity();

          // Skip entries with zero luminosity, since they won't contribute
          // anything to the overall integral. Skip negative luminosity ones as
          // well, just in case.
          if (lum <= 0.) continue;

          flux += lum * current_source->pdf(E_nu);
        }
        return flux;
      }
    );

    double flux_integ = 0.;
    double tot_xs_integ = 0.;
    flux_averaged_total_xs(*temp_source, gen, flux_integ, tot_xs_integ);

    // Factor of hbar_c^2 converts from MeV^(-2) to fm^2
    fFluxAveragedCrossSection =  marley_utils::hbar_c2
      * tot_xs_integ / flux_integ;

    make_nu_emission_histograms();

  } // spectrum_file_format == "fit"

  else {
    throw cet::exception("MARLEYTimeGen") << "Unrecognized neutrino spectrum"
      << " file format \"" << p().spectrum_file_format_() << "\" encountered"
      << " in evgen::MarleyTimeGen::reconfigure()";
  }

  MF_LOG_INFO("MARLEYTimeGen") << "The flux-averaged total cross section"
    << " predicted by MARLEY for the current supernova spectrum is "
    << fFluxAveragedCrossSection << " fm^2";

}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::create_truths_time_fit(simb::MCTruth& mc_truth,
  sim::SupernovaTruth& sn_truth, const TLorentzVector& vertex_pos)
{
  // Get a reference to the generator object created by MARLEY (we'll need
  // to do a few fancy things with it other than just creating events)
  marley::Generator& gen = fMarleyHelper->get_generator();

  // Initialize the time bin index to something absurdly large. This will help
  // us detect strange bugs that arise when it is sampled incorrectly.
  size_t time_bin_index = std::numeric_limits<size_t>::max();

  // Create an object to represent the discrete time distribution given in
  // the spectrum file. Use the luminosity for each bin as its sampling weight.
  // This distribution will actually be used to sample a neutrino arrival time
  // unless we're using the uniform time sampling mode, in which case it will
  // be used to help calculate the neutrino vertex weight.
  int source_pdg_code = gen.get_source().get_pid();

  const auto fit_params_begin = TimeFit::make_FitParameters_iterator(
    source_pdg_code, fTimeFits.begin() );
  const auto fit_params_end = TimeFit::make_FitParameters_iterator(
    source_pdg_code, fTimeFits.end() );

  const auto lum_begin = FitParameters::make_luminosity_iterator(
    fit_params_begin);
  const auto lum_end = FitParameters::make_luminosity_iterator(
    fit_params_end);

  std::discrete_distribution<size_t> time_dist(lum_begin, lum_end);

  if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM
    || fSamplingMode == TimeGenSamplingMode::UNIFORM_ENERGY)
  {
    time_bin_index = gen.sample_from_distribution(time_dist);
  }

  else if (fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME) {
    int last_time_index = 0;
    if (fTimeFits.size() > 0) last_time_index = fTimeFits.size() - 1;
    std::uniform_int_distribution<size_t> uid(0, last_time_index);
    // for c2: time_bin_index has already been declared
    time_bin_index = gen.sample_from_distribution(uid);
  }

  else {
    throw cet::exception("MARLEYTimeGen") << "Unrecognized sampling mode"
      << " encountered in evgen::MarleyTimeGen::produce()";
  }

  // Sample a time uniformly from within the selected time bin. Note that
  // the entries in fTimeFits use the time_ member to store the bin left
  // edges. The module creates fTimeFits in such a way that its last element
  // will always have luminosity_ == 0. (zero sampling weight), so we may
  // always add one to the sampled bin index without worrying about going
  // off the edge.
  double t_min = fTimeFits.at(time_bin_index).Time();
  double t_max = fTimeFits.at(time_bin_index + 1).Time();
  // sample a time on [ t_min, t_max )
  fTNu = gen.uniform_random_double(t_min, t_max, false);

  // Create a "beta-fit" neutrino source using the correct parameters for the
  // sampled time bin. This will be used to sample a neutrino energy unless
  // we're using the uniform time sampling mode. For uniform time sampling,
  // it will be used to determine the neutrino event weight.
  const auto& fit = fTimeFits.at(time_bin_index);
  std::unique_ptr<marley::NeutrinoSource> nu_source = source_from_time_fit(fit);

  if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM
    || fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME)
  {
    // Replace the generator's old source with the new one for the current
    // time bin
    gen.set_source(std::move(nu_source));

    // Generate a MARLEY event using the updated source
    mc_truth = fMarleyHelper->create_MCTruth(vertex_pos, fEvent.get());

    if (fSamplingMode == TimeGenSamplingMode::HISTOGRAM) {
      // Unbiased sampling creates neutrino vertices with unit weight
      fWeight = ONE;
      sn_truth = sim::SupernovaTruth(fTNu, fWeight, fFluxAveragedCrossSection,
        sim::kUnbiased);
    }
    else {
      // fSamplingMode == TimeGenSamplingMode::UNIFORM_TIME

      // Multiply by the likelihood ratio in order to correct for uniform
      // time sampling if we're using that biasing method
      double weight_bias = time_dist.probabilities().at(time_bin_index)
        / (t_max - t_min) * ( fTimeFits.back().Time()
        - fTimeFits.front().Time() );

      fWeight = weight_bias;
      sn_truth = sim::SupernovaTruth(fTNu, fWeight, fFluxAveragedCrossSection,
        sim::kUniformTime);
    }
  }

  else if (fSamplingMode == TimeGenSamplingMode::UNIFORM_ENERGY)
  {
    double E_nu = std::numeric_limits<double>::lowest();
    mc_truth = make_uniform_energy_mctruth(fFitEmin, fFitEmax, E_nu,
      vertex_pos);

    // Get the value of the true dependent probability density (probability
    // of the sampled energy given the sampled time) to use as a biasing
    // correction in the neutrino vertex weight.

    // Load the generator with the neutrino source that represents the
    // true (i.e., unbiased) energy probability distribution. This will
    // create a normalized probability density that we can use to determine
    // the neutrino vertex weight.
    double nu_source_E_min = nu_source->get_Emin();
    double nu_source_E_max = nu_source->get_Emax();
    gen.set_source(std::move(nu_source));

    // NOTE: The marley::Generator object normalizes the E_pdf to unity
    // automatically, but just in case, we redo it here.
    double E_pdf_integ = integrate([&gen](double E_nu)
      -> double { return gen.E_pdf(E_nu); }, nu_source_E_min,
      nu_source_E_max);

    // Compute the likelihood ratio that we need to bias the neutrino vertex
    // weight
    double weight_bias = (gen.E_pdf(E_nu) / E_pdf_integ)
      * (fFitEmax - fFitEmin);

    fWeight = weight_bias;
    sn_truth = sim::SupernovaTruth(fTNu, fWeight, fFluxAveragedCrossSection,
      sim::kUniformEnergy);
  }

  else {
    throw cet::exception("MARLEYTimeGen") << "Unrecognized sampling mode"
      << " encountered in evgen::MarleyTimeGen::produce()";
  }

}

//------------------------------------------------------------------------------
std::unique_ptr<marley::NeutrinoSource>
  evgen::MarleyTimeGen::source_from_time_fit(const TimeFit& fit)
{
  // Create a "beta-fit" neutrino source using the given fit parameters.

  // The two common fitting schemes (alpha and beta) differ in their
  // definitions by \beta = \alpha + 1.
  // TODO: add handling for non-nu_e neutrino types
  const FitParameters& nue_parameters
    = fit.GetFitParameters(marley_utils::ELECTRON_NEUTRINO);

  double beta = nue_parameters.Alpha();
  if (fPinchType == PinchParamType::ALPHA) beta += 1.;
  else if (fPinchType != PinchParamType::BETA) {
    throw cet::exception("MARLEYTimeGen") << "Unreognized pinching parameter"
      << " type encountered in evgen::MarleyTimeGen::source_from_time_fit()";
  }

  // Create the new source
  std::unique_ptr<marley::NeutrinoSource> nu_source
    = std::make_unique<marley::BetaFitNeutrinoSource>(
    marley_utils::ELECTRON_NEUTRINO, fFitEmin, fFitEmax,
    nue_parameters.Emean(), beta);

  return nu_source;
}

//------------------------------------------------------------------------------
simb::MCTruth evgen::MarleyTimeGen::make_uniform_energy_mctruth(double E_min,
  double E_max, double& E_nu, const TLorentzVector& vertex_pos)
{
  marley::Generator& gen = fMarleyHelper->get_generator();

  // Sample an energy uniformly over the entire allowed range
  double total_xs;
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
  auto mc_truth = fMarleyHelper->create_MCTruth(vertex_pos, fEvent.get());

  return mc_truth;
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::make_final_timefit(double time)
{
  // The final time bin has zero luminosity, and therefore zero sampling
  // weight. We need it to be present so that the last nonzero weight bin
  // has a right edge.
  fTimeFits.emplace_back(time, 0., 0., 0., 0., 0., 0., 0., 0., 0.);
}

//------------------------------------------------------------------------------
void evgen::MarleyTimeGen::make_nu_emission_histograms() const
{
  // To make a histogram with variable size bins, ROOT needs the bin
  // low edges to be contiguous in memory. This is not true of the
  // stored times in fTimeFits, so we'll need to make a temporary copy
  // of them.
  std::vector<double> temp_times;
  for (const auto& fit : fTimeFits) temp_times.push_back(fit.Time());

  // If, for some reason, there are fewer than two time bins, just return
  // without making the histograms.
  // TODO: consider throwing an exception here
  if (temp_times.size() < 2) return;

  // Get the number of time bins
  int num_bins = temp_times.size() - 1;

  // Create some ROOT TH1D objects using the TFileService. These will store
  // the number of emitted neutrinos of each type in each time bin
  art::ServiceHandle<art::TFileService const> tfs;
  TH1D* nue_hist = tfs->make<TH1D>("NueEmission",
    "Number of emitted #nu_{e}; time (s); number of neutrinos in time bin",
    num_bins, temp_times.data());
  TH1D* nuebar_hist = tfs->make<TH1D>("NuebarEmission",
    "Number of emitted #bar{#nu}_{e}; time (s); number of neutrinos in"
    " time bin", num_bins, temp_times.data());
  TH1D* nux_hist = tfs->make<TH1D>("NuxEmission",
    "Number of emitted #nu_{x}; time (s); number of neutrinos in time bin",
    num_bins, temp_times.data());

  // Load the histograms with the emitted neutrino counts
  for (size_t b = 1; b < temp_times.size(); ++b) {
    const TimeFit& current_fit = fTimeFits.at(b - 1);
    const TimeFit& next_fit = fTimeFits.at(b);
    double bin_deltaT = next_fit.Time() - current_fit.Time();

    const auto& nue_params = current_fit.GetFitParameters(
      marley_utils::ELECTRON_NEUTRINO);
    const auto& nuebar_params = current_fit.GetFitParameters(
      marley_utils::ELECTRON_ANTINEUTRINO);
    const auto& nux_params = current_fit.GetFitParameters(
      marley_utils::MUON_NEUTRINO);

    // Convert from bin luminosity to number of neutrinos by
    // multiplying by the bin time interval and dividing by the
    // mean neutrino energy
    double num_nue = 0.;
    double num_nuebar = 0.;
    double num_nux = 0.;
    if (nue_params.Emean() != 0.) num_nue = nue_params.Luminosity()
      * bin_deltaT / (nue_params.Emean() * MeV_to_erg);
    if (nuebar_params.Emean() != 0.) num_nuebar = nuebar_params.Luminosity()
      * bin_deltaT / (nuebar_params.Emean() * MeV_to_erg);
    if (nux_params.Emean() != 0.) num_nux = nux_params.Luminosity()
      * bin_deltaT / (nux_params.Emean() * MeV_to_erg);

    nue_hist->SetBinContent(b, num_nue);
    nuebar_hist->SetBinContent(b, num_nuebar);
    nux_hist->SetBinContent(b, num_nux);
  }

}


//------------------------------------------------------------------------------
// Anonymous namespace function definitions
namespace {

  double flux_averaged_total_xs(marley::NeutrinoSource& nu_source,
    marley::Generator& gen, double& source_integ, double& tot_xs_integ)
  {
    // Get an integral of the source PDF (in case it isn't normalized to 1)
    source_integ = integrate(
      [&nu_source](double E_nu) -> double { return nu_source.pdf(E_nu); },
      nu_source.get_Emin(), nu_source.get_Emax()
    );

    tot_xs_integ = integrate(
      [&nu_source, &gen](double E_nu) -> double
      {
        double xs = 0.;
        for (const auto& react : gen.get_reactions()) {
          xs += react->total_xs(marley_utils::ELECTRON_NEUTRINO, E_nu);
        }
        return xs * nu_source.pdf(E_nu);
      }, nu_source.get_Emin(), nu_source.get_Emax()
    );

    return tot_xs_integ / source_integ;
  }

  double flux_averaged_total_xs(marley::NeutrinoSource& nu_source,
    marley::Generator& gen)
  {
    double dummy1, dummy2;
    return flux_averaged_total_xs(nu_source, gen, dummy1, dummy2);
  }

  double integrate(const std::function<double(double)>& f, double x_min,
    double x_max)
  {
    static marley::Integrator integrator(NUM_INTEGRATION_SAMPLING_POINTS);
    return integrator.num_integrate(f, x_min, x_max);
  }

}

DEFINE_ART_MODULE(evgen::MarleyTimeGen)
