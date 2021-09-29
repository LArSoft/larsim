////////////////////////////////////////////////////////////////////////
// Class:       PDFastSimPAR
// Plugin Type: producer
// File:        PDFastSimPAR_module.cc
// Description:
// - acts on sim::SimEnergyDeposit from LArG4Main,
// - simulate (fast, photon visibility service) the OpDet response to optical
// photons Input: 'sim::SimEnergyDeposit' Output: 'sim::OpDetBacktrackerRecord'
// Fast simulation of propagating the photons created from SimEnergyDeposits.

// This module does a fast simulation of propagating the photons created from
// SimEnergyDeposits, This simulation is done using the PhotonLibrary, which
// stores the visibilities of each optical channel with respect to each optical
// voxel in the TPC volume, to avoid propagating single photons using Geant4. At
// the end of this module a collection of the propagated photons either as
//'sim::OpDetBacktrackerRecord' are placed into the art event.

// The steps this module takes are:
//  - to take number of photon and the vertex information from
//  'sim::SimEnergyDeposits',
//  - use the PhotonLibrary (visibilities) to determine the amount of visible
//  photons at each optical channel,
//  - visible photons: the number of photons times the visibility at the middle
//  of the Geant4 step for a given optical channel.
//  - other photon information is got from 'sim::SimEnergyDeposits'
//  - add 'sim::OpDetBacktrackerRecord' to event
// Aug. 19 by Mu Wei
////////////////////////////////////////////////////////////////////////

// Art libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"         // geo::vect::fillCoords()
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h" // geo::vect::toTVector3()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/PhotonPropagation/PhotonVisibilityTypes.h" // phot::MappedT0s_t
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"

#include "larsim/IonizationScintillation/ISTPC.h"

// Random numbers
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
//#include "CLHEP/Random/RandGauss.h"

#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4EmProcessSubType.hh"
#include "Geant4/G4EmSaturation.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4OpticalPhoton.hh"
#include "Geant4/G4ParticleMomentum.hh"
#include "Geant4/G4PhysicsOrderedFreeVector.hh"
#include "Geant4/G4PhysicsTable.hh"
#include "Geant4/G4Poisson.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4VRestDiscreteProcess.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/globals.hh"
#include "Geant4/templates.hh"

// support libraries
#include "cetlib_except/exception.h"

#include "Math/SpecFuncMathMore.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"
#include "TF1.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

// Define a new policy *not* internally promoting RealType to double:
typedef boost::math::policies::policy<
  // boost::math::policies::digits10<8>,
  boost::math::policies::promote_double<false>>
  noLDoublePromote;

namespace {
  //......................................................................
  double
  finter_d(double* x, double* par)
  {
    double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
    double y2 = TMath::Exp(par[3] + x[0] * par[4]);

    return TMath::Abs(y1 - y2);
  }
  //  double LandauPlusExpoFinal(double *x, double *par)
  //  {
  //        // par0 = joining point
  //        // par1 = Landau MPV
  //        // par2 = Landau widt
  //        // par3 = normalization
  //        // par4 = Expo cte
  //        // par5 = Expo tau
  //        double y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
  //        double y2 = TMath::Exp(par[4]+x[0]*par[5]);
  //        if(x[0] > par[0]) y1 = 0.;
  //        if(x[0] < par[0]) y2 = 0.;
  //
  //        return (y1 + y2);
  //  }

  //......................................................................
  //  double finter_r(double *x, double *par)
  //  {
  //        double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
  //        double y2 = par[5]*TMath::Landau(x[0],par[3],par[4]);
  //
  //        return TMath::Abs(y1 - y2);
  //  }

  double
  model_close(double* x, double* par)
  {
    // par0 = joining point
    // par1 = Landau MPV
    // par2 = Landau width
    // par3 = normalization
    // par4 = Expo cte
    // par5 = Expo tau
    // par6 = t_min

    double y1 = par[3] * TMath::Landau(x[0], par[1], par[2]);
    double y2 = TMath::Exp(par[4] + x[0] * par[5]);
    if (x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
    if (x[0] < par[0]) y2 = 0.;

    return (y1 + y2);
  }

  //......................................................................
  double
  model_far(double* x, double* par)
  {
    // par1 = Landau MPV
    // par2 = Landau width
    // par3 = normalization
    // par0 = t_min

    double y = par[3] * TMath::Landau(x[0], par[1], par[2]);
    if (x[0] <= par[0]) y = 0.;

    return y;
  }

  //......................................................................
  double
  Gaisser_Hillas(const double x, const double* par)
  {
    double X_mu_0 = par[3];
    double Normalization = par[0];
    double Diff = par[1] - X_mu_0;
    double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
    double Exponential = std::exp((par[1] - x) / par[2]);

    return (Normalization * Term * Exponential);
  }

  //......................................................................
  double
  fast_acos(double x)
  {
    double negate = double(x < 0);
    x = std::abs(x);
    x -= double(x > 1.0) * (x - 1.0); // <- equivalent to min(1.0,x), but faster
    double ret = -0.0187293;
    ret = ret * x;
    ret = ret + 0.0742610;
    ret = ret * x;
    ret = ret - 0.2121144;
    ret = ret * x;
    ret = ret + 1.5707288;
    ret = ret * std::sqrt(1.0 - x);
    ret = ret - 2. * negate * ret;
    return negate * 3.14159265358979 + ret;
  }

  // implements relative method - do not use for comparing with zero
  // use this most of the time, tolerance needs to be meaningful in your context
  template <typename TReal>
  inline constexpr static bool
  isApproximatelyEqual(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    TReal diff = std::fabs(a - b);
    if (diff <= tolerance) return true;
    if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
    return false;
  }

  // supply tolerance that is meaningful in your context
  // for example, default tolerance may not work if you are comparing double with
  // float
  template <typename TReal>
  inline constexpr static bool
  isApproximatelyZero(TReal a, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    if (std::fabs(a) <= tolerance) return true;
    return false;
  }

  // use this when you want to be on safe side
  // for example, don't start rover unless signal is above 1
  template <typename TReal>
  inline constexpr static bool
  isDefinitelyLessThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    TReal diff = a - b;
    if (diff < tolerance) return true;
    if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
    return false;
  }

  template <typename TReal>
  inline constexpr static bool
  isDefinitelyGreaterThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    TReal diff = a - b;
    if (diff > tolerance) return true;
    if (diff > std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
    return false;
  }

} // namespace

namespace phot {
  class PDFastSimPAR : public art::EDProducer {
  public:

  // Define the fhicl configuration 
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      using DP  = fhicl::DelegatedParameter;
      using ODP = fhicl::OptionalDelegatedParameter;

      fhicl::Atom<art::InputTag> SimulationLabel  { Name("SimulationLabel"),  Comment("SimEnergyDeposit label.") };
      fhicl::Atom<bool>          DoFastComponent  { Name("DoFastComponent"),  Comment("Simulate slow scintillation light, default true"), true };
      fhicl::Atom<bool>          DoSlowComponent  { Name("DoSlowComponent"),  Comment("Simulate slow scintillation light") };
      fhicl::Atom<bool>          DoReflectedLight { Name("DoReflectedLight"), Comment("Simulate reflected visible light") };
      fhicl::Atom<bool>          IncludeAnodeReflections { Name("IncludeAnodeReflections"), Comment("Simulate anode reflections, default false"), false };
      fhicl::Atom<bool>          IncludePropTime  { Name("IncludePropTime"),  Comment("Simulate light propagation time") };
      fhicl::Atom<bool>          GeoPropTimeOnly  { Name("GeoPropTimeOnly"),  Comment("Simulate light propagation time geometric approximation, default false"), false };
      fhicl::Atom<bool>          UseLitePhotons   { Name("UseLitePhotons"),   Comment("Store SimPhotonsLite/OpDetBTRs instead of SimPhotons") };
      fhicl::Atom<bool>          OpaqueCathode    { Name("OpaqueCathode"),    Comment("Photons cannot cross the cathode") };
      fhicl::Atom<bool>          OnlyActiveVolume { Name("OnlyActiveVolume"), Comment("PAR fast sim usually only for active volume, default true"), true };
      fhicl::Atom<bool>          OnlyOneCryostat  { Name("OnlyOneCryostat"),  Comment("Set to true if light is only supported in C:1") };
      DP                         ScintTimeTool    { Name("ScintTimeTool"),    Comment("Tool describing scintillation time structure")}; 
      ODP                        VUVTiming        { Name("VUVTiming"),        Comment("Configuration for UV timing parameterization")}; 
      ODP                        VISTiming        { Name("VISTiming"),        Comment("Configuration for visible timing parameterization")}; 
      DP                         VUVHits          { Name("VUVHits"),          Comment("Configuration for UV visibility parameterization")}; 
      ODP                        VISHits          { Name("VISHits"),          Comment("Configuration for visibile visibility parameterization")}; 

      

    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit PDFastSimPAR(Parameters const & config);
    void produce(art::Event&) override;

  private:
    larg4::ISTPC fISTPC;
    // structure definition for solid angle of rectangle function
    struct Dims {
      double h, w; // height, width
    };

    struct OpticalDetector {
      double h; // height
      double w; // width
      geo::Point_t OpDetPoint;
      int type;
      int orientation;
    };

    void Initialization();

    void getVUVTimes(std::vector<double>& arrivalTimes, const double distance_in_cm, const size_t angle_bin);
    void getVUVTimesGeo(std::vector<double>& arrivalTimes, const double distance_in_cm);
    void getVISTimes(std::vector<double>& arrivalTimes, const TVector3 &ScintPoint, const TVector3 &OpDetPoint);

    void generateParam(const size_t index, const size_t angle_bin);

    void AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                     std::map<size_t, int>& ChannelMap,
                     sim::OpDetBacktrackerRecord btr);

    void detectedDirectHits(std::map<size_t, int>& DetectedNumFast,
                            std::map<size_t, int>& DetectedNumSlow,
                            const double NumFast,
                            const double NumSlow,
                            geo::Point_t const& ScintPoint);
    void detectedReflecHits(std::map<size_t, int>& ReflDetectedNumFast,
                            std::map<size_t, int>& ReflDetectedNumSlow,
                            const double NumFast,
                            const double NumSlow,
                            geo::Point_t const& ScintPoint,
                            bool AnodeMode = false);
    
    void VUVHits(const double NumFast,
                const double NumSlow,
                geo::Point_t const& ScintPoint,
                OpticalDetector const& opDet,
                std::vector<int> &DetThis);

    void VISHits(geo::Point_t const& ScintPoint,
                OpticalDetector const& opDet,
                const double cathode_hits_rec_fast,
                const double cathode_hits_rec_slow,
                geo::Point_t const& hotspot,
                std::vector<int> &ReflDetThis,
                bool AnodeMode = false);

    void propagationTime(std::vector<double>& arrival_time_dist,
                         geo::Point_t const& x0,
                         const size_t OpChannel,
                         bool Reflected = false); // const;

    double interpolate(const std::vector<double>& xData,
                       const std::vector<double>& yData,
                       double x,
                       bool extrapolate,
                       size_t i = 0);
    void interpolate3(std::array<double, 3>& inter,
                      const std::vector<double>& xData,
                      const std::vector<double>& yData1,
                      const std::vector<double>& yData2,
                      const std::vector<double>& yData3,
                      double x,
                      bool extrapolate);

    // solid angle of rectangular aperture calculation functions
    double Rectangle_SolidAngle(const double a, const double b, const double d);
    double Rectangle_SolidAngle(Dims const& o, geo::Vector_t const& v, double OpDetOrientation);
    // solid angle of circular aperture calculation functions
    double Disk_SolidAngle(const double d, const double h, const double b);
     // solid angle of a dome aperture calculation functions
    double Omega_Dome_Model(const double distance, const double theta) const;


    CLHEP::HepRandomEngine& fPhotonEngine;
    std::unique_ptr<CLHEP::RandPoissonQ> fRandPoissPhot;
    CLHEP::HepRandomEngine& fScintTimeEngine;

    size_t nOpDets; // Pulled from geom during Initialization()

    std::map<size_t, int> PDChannelToSOCMapDirect; // Where each OpChan is.
    std::map<size_t, int> PDChannelToSOCMapReflect; // Where each OpChan is.


    // geometry properties
    double fplane_depth, fcathode_zdimension, fcathode_ydimension;
    double fanode_plane_depth, fanode_ydimension, fanode_zdimension;
    TVector3 fcathode_centre;
    TVector3 fanode_centre;
    std::vector<geo::BoxBoundedGeo> fActiveVolumes;
    int fNTPC;
    // Optical detector properties for semi-analytic hits
    double fradius;
    Dims fcathode_plane;
    Dims fanode_plane;
    int fL_abs_vuv;
    std::vector<geo::Point_t> fOpDetCenter;
    std::vector<int> fOpDetType;
    std::vector<int> fOpDetOrientation;
    std::vector<double> fOpDetLength;
    std::vector<double> fOpDetHeight;


    bool isOpDetInSameTPC(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) const;
    

    //////////////////////
    // Input Parameters //
    //////////////////////

    // Module behavior
    art::InputTag simTag;
    bool fDoFastComponent;
    bool fDoSlowComponent;
    bool fDoReflectedLight;
    bool fIncludeAnodeReflections;
    bool fIncludePropTime;
    bool fGeoPropTimeOnly;
    bool fUseLitePhotons;
    bool fOpaqueCathode;
    bool fOnlyActiveVolume;
    bool fOnlyOneCryostat;
    std::unique_ptr<ScintTime> fScintTime; // Tool to retrive timinig of scintillation

    // Parameterized Simulation
    fhicl::ParameterSet fVUVTimingParams;
    fhicl::ParameterSet fVISTimingParams;
    fhicl::ParameterSet fVUVHitsParams;
    fhicl::ParameterSet fVISHitsParams;

    // For VUV transport time parametrization
    double fstep_size, fmax_d, fmin_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance, fangle_bin_timing_vuv;
    std::vector<std::vector<double>> fparameters[7];
    // vector containing generated VUV timing parameterisations
    std::vector<std::vector<TF1>> VUV_timing;
    // vector containing min and max range VUV timing parameterisations are sampled to
    std::vector<std::vector<double>> VUV_max;
    std::vector<std::vector<double>> VUV_min;

    // For VIS transport time parameterisation
    double fvis_vmean, fangle_bin_timing_vis;
    std::vector<double> fdistances_refl;
    std::vector<double> fradial_distances_refl;
    std::vector<std::vector<std::vector<double>>> fcut_off_pars;
    std::vector<std::vector<std::vector<double>>> ftau_pars;

    // For VUV semi-analytic hits
    double fdelta_angulo_vuv;
    // flat PDs
    bool fIsFlatPDCorr;
    std::vector<std::vector<double>> fGHvuvpars_flat;
    std::vector<double> fborder_corr_angulo_flat;
    std::vector<std::vector<double>> fborder_corr_flat;
    // lateral PDs
    bool fIsFlatPDCorrLat;
    std::vector<double> fGH_distances_anode;
    std::vector<std::vector<std::vector<double>>> fGHvuvpars_flat_lateral;    
    // dome PDs
    bool fIsDomePDCorr;
    std::vector<std::vector<double>> fGHvuvpars_dome;
    std::vector<double> fborder_corr_angulo_dome;
    std::vector<std::vector<double>> fborder_corr_dome;
    // Field cage scaling
    bool fApplyFieldCageTransparency;
    double fFieldCageTransparencyLateral;
    double fFieldCageTransparencyCathode;

    // For VIS semi-analytic hits
    // correction parameters for VIS Nhits estimation
    double fdelta_angulo_vis;
    double fAnodeReflectivity;
    // flat PDs
    std::vector<double> fvis_distances_x_flat;
    std::vector<double> fvis_distances_r_flat;
    std::vector<std::vector<std::vector<double>>> fvispars_flat;
    // lateral PDs
    std::vector<double> fvis_distances_x_flat_lateral;
    std::vector<double> fvis_distances_r_flat_lateral;
    std::vector<std::vector<std::vector<double>>> fvispars_flat_lateral;
    // dome PDs
    std::vector<double> fvis_distances_x_dome;
    std::vector<double> fvis_distances_r_dome;
    std::vector<std::vector<std::vector<double>>> fvispars_dome;
    

  };

  //......................................................................
  PDFastSimPAR::PDFastSimPAR(Parameters const & config)
    : art::EDProducer{config}
    , fISTPC{*(lar::providerFrom<geo::Geometry>())}
    , fPhotonEngine(   art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, 
                                                                                 "HepJamesRandom",
                                                                                 "photon", 
                                                                                 config.get_PSet(), 
                                                                                 "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,
                                                                                 "HepJamesRandom",
                                                                                 "scinttime",
                                                                                 config.get_PSet(),
                                                                                 "SeedScintTime"))
    , simTag(config().SimulationLabel())
    , fDoFastComponent(config().DoFastComponent())
    , fDoSlowComponent(config().DoSlowComponent())
    , fDoReflectedLight(config().DoReflectedLight())
    , fIncludeAnodeReflections(config().IncludeAnodeReflections())
    , fIncludePropTime(config().IncludePropTime())
    , fGeoPropTimeOnly(config().GeoPropTimeOnly())
    , fUseLitePhotons(config().UseLitePhotons())
    , fOpaqueCathode(config().OpaqueCathode())
    , fOnlyActiveVolume(config().OnlyActiveVolume())
    , fOnlyOneCryostat(config().OnlyOneCryostat())
    , fScintTime{art::make_tool<ScintTime>(config().ScintTimeTool.get<fhicl::ParameterSet>())}
    , fVUVHitsParams(config().VUVHits.get<fhicl::ParameterSet>())
  {

    // Validate configuration options
    if(fIncludePropTime && !config().VUVTiming.get_if_present<fhicl::ParameterSet>(fVUVTimingParams)) {
      throw art::Exception(art::errors::Configuration)
          << "Propagation time simulation requested, but VUVTiming not specified." << "\n";
    }
    
    if(fDoReflectedLight && !config().VISHits.get_if_present<fhicl::ParameterSet>(fVISHitsParams)) {
      throw art::Exception(art::errors::Configuration)
          << "Reflected light simulation requested, but VisHits not specified." << "\n";
    }

    if (fDoReflectedLight && fIncludePropTime && !config().VISTiming.get_if_present<fhicl::ParameterSet>(fVISTimingParams)) {
      throw art::Exception(art::errors::Configuration)
          << "Reflected light propagation time simulation requested, but VISTiming not specified." << "\n";
    }

    if(fIncludeAnodeReflections && !config().VISHits.get_if_present<fhicl::ParameterSet>(fVISHitsParams)) {
      throw art::Exception(art::errors::Configuration)
          << "Anode reflections light simulation requested, but VisHits not specified." << "\n";
    }   


    Initialization();
    if (fUseLitePhotons)
    {
        mf::LogInfo("PDFastSimPAR") << "Using Lite Photons";
        produces< std::vector<sim::SimPhotonsLite> >();
        produces< std::vector<sim::OpDetBacktrackerRecord> >();

        if(fDoReflectedLight)
        {
            mf::LogInfo("PDFastSimPAR") << "Storing Reflected Photons";
            produces< std::vector<sim::SimPhotonsLite> >("Reflected");
            produces< std::vector<sim::OpDetBacktrackerRecord> >("Reflected");
        }
    }
    else
    {
        mf::LogInfo("PDFastSimPAR") << "Using Sim Photons";
        produces< std::vector<sim::SimPhotons> >();
        if(fDoReflectedLight)
        {
            mf::LogInfo("PDFastSimPAR") << "Storing Reflected Photons";
            produces< std::vector<sim::SimPhotons> >("Reflected");
        }
    }
  }

  //......................................................................
  void
  PDFastSimPAR::produce(art::Event& event)
  {
    mf::LogTrace("PDFastSimPAR") << "PDFastSimPAR Module Producer"
                                 << "EventID: " << event.event();

    auto phot = std::make_unique<std::vector<sim::SimPhotons>>();
    auto phlit = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto opbtr = std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();

    auto phot_ref = std::make_unique<std::vector<sim::SimPhotons>>();
    auto phlit_ref = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto opbtr_ref = std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();

    auto& dir_photcol(*phot);
    auto& ref_photcol(*phot_ref);
    auto& dir_phlitcol(*phlit);
    auto& ref_phlitcol(*phlit_ref);
    dir_photcol.resize(nOpDets);
    ref_photcol.resize(nOpDets);
    dir_phlitcol.resize(nOpDets);
    ref_phlitcol.resize(nOpDets);
    for (unsigned int i = 0; i < nOpDets; i ++)
    {
        dir_photcol[i].fOpChannel  = i;
        ref_photcol[i].fOpChannel  = i;
        dir_phlitcol[i].OpChannel  = i;
        ref_phlitcol[i].OpChannel  = i;
    }

    art::Handle<std::vector<sim::SimEnergyDeposit>> edepHandle;
    if (!event.getByLabel(simTag, edepHandle)) {
      mf::LogError("PDFastSimPAR") << "PDFastSimPAR Module Cannot getByLabel: " << simTag;
      return;
    }


    auto const& edeps = edepHandle;

    int num_points = 0;
    int num_fastph = 0;
    int num_slowph = 0;
    int num_fastdp = 0;
    int num_slowdp = 0;

    for (auto const& edepi : *edeps) {
      num_points++;

      int trackID = edepi.TrackID();
      double nphot = edepi.NumPhotons();
      double edeposit = edepi.Energy() / nphot;
      double pos[3] = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
      geo::Point_t const ScintPoint = {pos[0], pos[1], pos[2]};

      if (fOnlyActiveVolume && !fISTPC.isScintInActiveVolume(ScintPoint)) continue;

      double nphot_fast = edepi.NumFPhotons();
      double nphot_slow = edepi.NumSPhotons();

      num_fastph += nphot_fast;
      num_slowph += nphot_slow;

      // direct light
      std::map<size_t, int> DetectedNumFast;
      std::map<size_t, int> DetectedNumSlow;

      bool needHits = (nphot_fast > 0 && fDoFastComponent) || (nphot_slow > 0 && fDoSlowComponent);
      if ( needHits ) {
        detectedDirectHits(DetectedNumFast, DetectedNumSlow, nphot_fast, nphot_slow, ScintPoint);
        if ( fIncludeAnodeReflections ) {
          std::map<size_t, int> AnodeDetectedNumFast;
          std::map<size_t, int> AnodeDetectedNumSlow;
          detectedReflecHits(AnodeDetectedNumFast, AnodeDetectedNumSlow, nphot_fast, nphot_slow, ScintPoint, true);
          // add to exiting count
          for (size_t const OpDet : util::counter(nOpDets)){
            DetectedNumFast[OpDet] += AnodeDetectedNumFast[OpDet];
            DetectedNumSlow[OpDet] += AnodeDetectedNumSlow[OpDet];
          }
        }
      }

      // reflected light, if enabled
      std::map<size_t, int> ReflDetectedNumFast;
      std::map<size_t, int> ReflDetectedNumSlow;
      if (fDoReflectedLight && needHits)
        detectedReflecHits(ReflDetectedNumFast, ReflDetectedNumSlow, nphot_fast, nphot_slow, ScintPoint);

      // propagation time
      std::vector<double> transport_time;

      // loop through direct photons then reflected photons cases
      for (size_t Reflected = 0; Reflected <= 1; ++Reflected) {

        // only do the reflected loop if including reflected light
        if (Reflected && !fDoReflectedLight) continue;

        for (size_t channel = 0; channel < nOpDets; channel++) {

          if (fOpaqueCathode && !isOpDetInSameTPC(ScintPoint, fOpDetCenter[channel])) continue;

          int ndetected_fast = DetectedNumFast[channel];
          int ndetected_slow = DetectedNumSlow[channel];
          if (Reflected) {
            ndetected_fast = ReflDetectedNumFast[channel];
            ndetected_slow = ReflDetectedNumSlow[channel];
          }

          // calculate propagation time, does not matter whether fast or slow photon
          transport_time.resize(ndetected_fast + ndetected_slow);
          if (fIncludePropTime && needHits)
            propagationTime(transport_time, ScintPoint, channel, Reflected);

          // SimPhotonsLite case
          if (fUseLitePhotons) {

            sim::OpDetBacktrackerRecord tmpbtr(channel);

            if (ndetected_fast > 0 && fDoFastComponent) {
              int n = ndetected_fast;
              num_fastdp += n;
              for (long i = 0; i < n; ++i) {
                // calculates the time at which the photon was produced
                fScintTime->GenScintTime(true, fScintTimeEngine);
                auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[i]);
                if (Reflected) ++ref_phlitcol[channel].DetectedPhotons[time];
                else ++dir_phlitcol[channel].DetectedPhotons[time];
                tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
              }
            }

            if (ndetected_slow > 0 && fDoSlowComponent) {
              int n = ndetected_slow;
              num_slowdp += n;
              for (long i = 0; i < n; ++i) {
                fScintTime->GenScintTime(false, fScintTimeEngine);
                auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[ndetected_fast + i]);
                if (Reflected) ++ref_phlitcol[channel].DetectedPhotons[time];
                else ++dir_phlitcol[channel].DetectedPhotons[time];
                tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
              }
            }

            if (Reflected) AddOpDetBTR(*opbtr_ref, PDChannelToSOCMapReflect, tmpbtr);
            else AddOpDetBTR(*opbtr, PDChannelToSOCMapDirect, tmpbtr);
          }
          // SimPhotons case
          else {

            sim::OnePhoton photon;
            photon.SetInSD         = false;
            photon.InitialPosition = edepi.End();
            if (Reflected) photon.Energy = 2.9 * CLHEP::eV; // 430 nm
            else photon.Energy = 9.7 * CLHEP::eV; // 128 nm

            if (ndetected_fast > 0 && fDoFastComponent) {
              int n = ndetected_fast;
              num_fastdp += n;
              for (long i = 0; i < n; ++i) {
                // calculates the time at which the photon was produced
                fScintTime->GenScintTime(true, fScintTimeEngine);
                auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[i]);
                photon.Time = time;
                if(Reflected) ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                else dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
              }
            }

            if (ndetected_slow > 0 && fDoSlowComponent) {
              int n = ndetected_slow;
              num_slowdp += n;
              for (long i = 0; i < n; ++i) {
                fScintTime->GenScintTime(false, fScintTimeEngine);
                auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + transport_time[ndetected_fast + i]);
                photon.Time = time;
                if(Reflected) ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                else dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
              }
            }
          }
        }
      }
    }

    mf::LogTrace("PDFastSimPAR") << "Total points: " << num_points
                                 << ", total fast photons: " << num_fastph
                                 << ", total slow photons: " << num_slowph
                                 << "\ndetected fast photons: " << num_fastdp
                                 << ", detected slow photons: " << num_slowdp;

    PDChannelToSOCMapDirect.clear();
    PDChannelToSOCMapReflect.clear();

    if (fUseLitePhotons) {
        event.put(move(phlit));
        event.put(move(opbtr));
        if (fDoReflectedLight) {
            event.put(move(phlit_ref), "Reflected");
            event.put(move(opbtr_ref), "Reflected");
        }
    }
    else {
        event.put(move(phot));
        if (fDoReflectedLight) {
            event.put(move(phot_ref), "Reflected");
        }
    }

    return;
  }

  //......................................................................
  void
  PDFastSimPAR::AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                            std::map<size_t, int>& ChannelMap,
                            sim::OpDetBacktrackerRecord btr
                            )
  {
    size_t iChan = btr.OpDetNum();
    auto channelPosition = ChannelMap.find(iChan);

    if (channelPosition == ChannelMap.end()) {
      ChannelMap[iChan] = opbtr.size();
      opbtr.emplace_back(std::move(btr));
    }
    else {
      unsigned int idtest = channelPosition->second;
      auto const& timePDclockSDPsMap = btr.timePDclockSDPsMap();

      for (auto const& timePDclockSDP : timePDclockSDPsMap) {
        for (auto const& sdp : timePDclockSDP.second) {
          double xyz[3] = {sdp.x, sdp.y, sdp.z};
          opbtr.at(idtest).AddScintillationPhotons(
            sdp.trackID, timePDclockSDP.first, sdp.numPhotons, xyz, sdp.energy);
        }
      }
    }
  }

  //......................................................................
  void
  PDFastSimPAR::Initialization()
  {
    std::cout << "PDFastSimPAR Initialization" << std::endl;
    std::cout << "Initializing the geometry of the detector." << std::endl;
    std::cout << "Simulate using semi-analytic model for number of hits." << std::endl;

    fRandPoissPhot = std::make_unique<CLHEP::RandPoissonQ>(fPhotonEngine);
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

    // Store info from the Geometry service
    nOpDets = geom.NOpDets();
    fActiveVolumes = fISTPC.extractActiveLArVolume(geom);
    fNTPC = geom.NTPC();

    {
      auto log = mf::LogTrace("PDFastSimPAR") << "PDFastSimPAR: active volume boundaries from "
                                              << fActiveVolumes.size() << " volumes:";
      for (auto const& [iCryo, box] : util::enumerate(fActiveVolumes)) {
        log << "\n - C:" << iCryo << ": " << box.Min() << " -- " << box.Max() << " cm";
      }
    } // local scope

    if (geom.Ncryostats() > 1U) {
      if (fOnlyOneCryostat) {
        mf::LogWarning("PDFastSimPAR")
          << std::string(80, '=') << "\nA detector with " << geom.Ncryostats()
          << " cryostats is configured"
          << " , and semi-analytic model is requested for scintillation photon propagation."
          << " THIS CONFIGURATION IS NOT SUPPORTED and it is open to bugs"
          << " (e.g. scintillation may be detected only in cryostat #0)."
          << "\nThis would be normally a fatal error, but it has been forcibly overridden."
          << "\n"
          << std::string(80, '=');
      }
      else {
        throw art::Exception(art::errors::Configuration)
          << "Photon propagation via semi-analytic model is not supported yet"
          << " on detectors with more than one cryostat.";
      }
    }

    fcathode_centre = {geom.TPC(0, 0).GetCathodeCenter().X(),
                       fActiveVolumes[0].CenterY(),
                       fActiveVolumes[0].CenterZ()};

    fanode_centre = {geom.TPC(0, 0).FirstPlane().GetCenter().X(),
                      fActiveVolumes[0].CenterY(),
                      fActiveVolumes[0].CenterZ() };

    for (size_t const i : util::counter(nOpDets)) {
      geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpDet(i);
      fOpDetCenter.push_back(opDet.GetCenter());

      if (opDet.isSphere()) {  // dome PMTs
        fOpDetType.push_back(1); // dome
        fOpDetOrientation.push_back(0); // anode/cathode (default)
        fOpDetLength.push_back(-1);
        fOpDetHeight.push_back(-1);
      }
      else if (opDet.isBar()) {
        fOpDetType.push_back(0); // (X)Arapucas/Bars
        // determine orientation to get correction OpDet dimensions
        fOpDetLength.push_back(opDet.Length());
        if (opDet.Width() > opDet.Height()) { // laterals, Y dimension smallest
          fOpDetOrientation.push_back(1);          
          fOpDetHeight.push_back(opDet.Width());  
        }
        else {  // anode/cathode (default), X dimension smallest
          fOpDetOrientation.push_back(0);
          fOpDetHeight.push_back(opDet.Height());
        }        
      }
      else {
        fOpDetType.push_back(2); // disk PMTs
        fOpDetOrientation.push_back(0); // anode/cathode (default)
        fOpDetLength.push_back(-1);
        fOpDetHeight.push_back(-1);
      }

    }

    if (fIncludePropTime && !fGeoPropTimeOnly) {
      mf::LogInfo("PDFastSimPAR") << "Using VUV timing parameterization";

      fparameters[0] = std::vector(1, fVUVTimingParams.get<std::vector<double>>("Distances_landau"));
      fparameters[1] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Norm_over_entries");
      fparameters[2] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Mpv");
      fparameters[3] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Width");
      fparameters[4] = std::vector(1, fVUVTimingParams.get<std::vector<double>>("Distances_exp"));
      fparameters[5] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Slope");
      fparameters[6] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Expo_over_Landau_norm");

      fstep_size                = fVUVTimingParams.get<double>("step_size");
      fmax_d                    = fVUVTimingParams.get<double>("max_d");
      fmin_d                    = fVUVTimingParams.get<double>("min_d");
      fvuv_vgroup_mean          = fVUVTimingParams.get<double>("vuv_vgroup_mean");
      fvuv_vgroup_max           = fVUVTimingParams.get<double>("vuv_vgroup_max");
      finflexion_point_distance = fVUVTimingParams.get<double>("inflexion_point_distance");
      fangle_bin_timing_vuv     = fVUVTimingParams.get<double>("angle_bin_timing_vuv");

      // create vector of empty TF1s that will be replaces with the parameterisations
      // that are generated as they are required
      // default TF1() constructor gives function with 0 dimensions, can then check
      // numDim to qucikly see if a parameterisation has been generated
      const size_t num_params = (fmax_d - fmin_d) / fstep_size; // for d < fmin_d, no parameterisaton, a delta function is used instead
      size_t num_angles = std::round(90/fangle_bin_timing_vuv);
      VUV_timing = std::vector(num_angles, std::vector(num_params, TF1()));

      // initialise vectors to contain range parameterisations sampled to in each case
      // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling
      // is regenerated, this is the slow part!
      VUV_max = std::vector(num_angles, std::vector(num_params, 0.0));
      VUV_min = std::vector(num_angles, std::vector(num_params, 0.0));

      // VIS time parameterisation
      if (fDoReflectedLight) {
        mf::LogInfo("PDFastSimPAR") << "Using VIS (reflected) timing parameterization";

        // load parameters
        fdistances_refl        = fVISTimingParams.get<std::vector<double>>("Distances_refl");
        fradial_distances_refl = fVISTimingParams.get<std::vector<double>>("Distances_radial_refl");
        fcut_off_pars          = fVISTimingParams.get<std::vector<std::vector<std::vector<double>>>>("Cut_off");
        ftau_pars              = fVISTimingParams.get<std::vector<std::vector<std::vector<double>>>>("Tau");
        fvis_vmean             = fVISTimingParams.get<double>("vis_vmean");
        fangle_bin_timing_vis  = fVISTimingParams.get<double>("angle_bin_timing_vis");
      }
    }
    if (fIncludePropTime && fGeoPropTimeOnly) {
      mf::LogInfo("PDFastSimPAR") << "Using geometric VUV time propagation";
      fvuv_vgroup_mean          = fVUVTimingParams.get<double>("vuv_vgroup_mean");
    }

    // LAr absorption length in cm
    std::map<double, double> abs_length_spectrum = lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
    std::vector<double> x_v, y_v;
    for (auto elem : abs_length_spectrum) {
      x_v.push_back(elem.first);
      y_v.push_back(elem.second);
    }
    fL_abs_vuv = std::round(interpolate(x_v, y_v, 9.7, false)); // 9.7 eV: peak of VUV emission spectrum

    // Load Gaisser-Hillas corrections for VUV semi-analytic hits
    mf::LogInfo("PDFastSimPAR") << "Using VUV visibility parameterization";
    
    fIsFlatPDCorr     = fVUVHitsParams.get<bool>("FlatPDCorr", false);
    fIsFlatPDCorrLat  = fVUVHitsParams.get<bool>("FlatPDCorrLat", false);
    fIsDomePDCorr     = fVUVHitsParams.get<bool>("DomePDCorr", false);
    fdelta_angulo_vuv = fVUVHitsParams.get<double>("delta_angulo_vuv", 10);
    fradius           = fVUVHitsParams.get<double>("PMT_radius", 10.16);
    fApplyFieldCageTransparency = fVUVHitsParams.get<bool>("ApplyFieldCageTransparency", false);
    fFieldCageTransparencyLateral = fVUVHitsParams.get<double>("FieldCageTransparencyLateral", 1.0);
    fFieldCageTransparencyCathode = fVUVHitsParams.get<double>("FieldCageTransparencyCathode", 1.0);


    if (!fIsFlatPDCorr && !fIsDomePDCorr && !fIsFlatPDCorrLat) {
      throw cet::exception("PDFastSimPAR")
          << "Both isFlatPDCorr/isFlatPDCorrLat and isDomePDCorr parameters are false, at least one type of parameterisation is required for the semi-analytic light simulation." << "\n";
    }
    if (fIsFlatPDCorr) {
        fGHvuvpars_flat          = fVUVHitsParams.get<std::vector<std::vector<double>>>("GH_PARS_flat");
        fborder_corr_angulo_flat = fVUVHitsParams.get<std::vector<double>>("GH_border_angulo_flat");
        fborder_corr_flat        = fVUVHitsParams.get<std::vector<std::vector<double>>>("GH_border_flat");
    }
    if (fIsFlatPDCorrLat) {
        fGHvuvpars_flat_lateral  = fVUVHitsParams.get<std::vector<std::vector<std::vector<double>>>>("GH_PARS_flat_lateral");
        fGH_distances_anode      = fVUVHitsParams.get<std::vector<double>>("GH_distances_anode");
    }
    if (fIsDomePDCorr) {
        fGHvuvpars_dome          = fVUVHitsParams.get<std::vector<std::vector<double>>>("GH_PARS_dome");
        fborder_corr_angulo_dome = fVUVHitsParams.get<std::vector<double>>("GH_border_angulo_dome");
        fborder_corr_dome        = fVUVHitsParams.get<std::vector<std::vector<double>>>("GH_border_dome");
    }

    // Load corrections for VIS semi-analytic hits
    if (fDoReflectedLight) {
      mf::LogInfo("PDFastSimPAR") << "Using VIS (reflected) visibility parameterization";
      fdelta_angulo_vis = fVISHitsParams.get<double>("delta_angulo_vis");

      if (fIsFlatPDCorr) {
        fvis_distances_x_flat = fVISHitsParams.get<std::vector<double>>("VIS_distances_x_flat");
        fvis_distances_r_flat = fVISHitsParams.get<std::vector<double>>("VIS_distances_r_flat");
        fvispars_flat         = fVISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_flat");
      }
      if (fIsDomePDCorr) {
        fvis_distances_x_dome = fVISHitsParams.get<std::vector<double>>("VIS_distances_x_dome");
        fvis_distances_r_dome = fVISHitsParams.get<std::vector<double>>("VIS_distances_r_dome");
        fvispars_dome         = fVISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_dome");
      }

      // cathode dimensions
      fcathode_ydimension = fActiveVolumes[0].SizeY();
      fcathode_zdimension = fActiveVolumes[0].SizeZ();

      // set cathode plane struct for solid angle function
      fcathode_plane.h = fcathode_ydimension;
      fcathode_plane.w = fcathode_zdimension;
      fplane_depth = std::abs(fcathode_centre[0]);
    }

    // Load corrections for Anode reflections configuration
    if (fIncludeAnodeReflections) {
      mf::LogInfo("PDFastSimPAR") << "Using anode reflections parameterization";
      fdelta_angulo_vis = fVISHitsParams.get<double>("delta_angulo_vis");
      fAnodeReflectivity = fVISHitsParams.get<double>("AnodeReflectivity");

      if (fIsFlatPDCorr) {
        fvis_distances_x_flat = fVISHitsParams.get<std::vector<double>>("VIS_distances_x_flat");
        fvis_distances_r_flat = fVISHitsParams.get<std::vector<double>>("VIS_distances_r_flat");
        fvispars_flat         = fVISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_flat");
      }

      if (fIsFlatPDCorrLat) {
        fvis_distances_x_flat_lateral = fVISHitsParams.get<std::vector<double>>("VIS_distances_x_flat_lateral");
        fvis_distances_r_flat_lateral = fVISHitsParams.get<std::vector<double>>("VIS_distances_r_flat_lateral");
        fvispars_flat_lateral         = fVISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_flat_lateral");
      }

      // anode dimensions
      fanode_ydimension = fActiveVolumes[0].SizeY();
      fanode_zdimension = fActiveVolumes[0].SizeZ();

      // set anode plane struct for solid angle function
      fanode_plane.h = fanode_ydimension;
      fanode_plane.w = fanode_zdimension;
      fanode_plane_depth = fanode_centre[0];
    }
  }

  //......................................................................
  // VUV semi-analytic hits calculation
  void
  PDFastSimPAR::detectedDirectHits(std::map<size_t, int>& DetectedNumFast,
                                   std::map<size_t, int>& DetectedNumSlow,
                                   const double NumFast,
                                   const double NumSlow,
                                   geo::Point_t const& ScintPoint)
  {
    for (size_t const OpDet : util::counter(nOpDets)) {
      if (!isOpDetInSameTPC(ScintPoint, fOpDetCenter[OpDet])) continue;

      // set detector struct for solid angle function
      const PDFastSimPAR::OpticalDetector op{
        fOpDetHeight[OpDet], fOpDetLength[OpDet],
        fOpDetCenter[OpDet], fOpDetType[OpDet], fOpDetOrientation[OpDet]};

      std::vector<int> DetThis(2, 0);
      VUVHits(NumFast, NumSlow, ScintPoint, op, DetThis);

      DetectedNumFast[OpDet] = DetThis[0];
      DetectedNumSlow[OpDet] = DetThis[1];

      //   mf::LogInfo("PDFastSimPAR") << "FastScint: " <<
      //   //   it->second<<" " << Num << " " << DetThisPMT;
      //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
    }
  }

  void
  PDFastSimPAR::VUVHits(const double NumFast,
                        const double NumSlow,
                        geo::Point_t const& ScintPoint,
                        OpticalDetector const& opDet,
                        std::vector<int>& DetThis)
  {
    // distance and angle between ScintPoint and OpDetPoint
    geo::Vector_t const relative = ScintPoint - opDet.OpDetPoint;
    const double distance = relative.R();
    double cosine;
    if (opDet.orientation == 1) cosine = std::abs(relative.Y()) / distance;
    else cosine = std::abs(relative.X()) / distance;
    // const double theta = std::acos(cosine) * 180. / CLHEP::pi;
    const double theta = fast_acos(cosine) * 180. / CLHEP::pi;

    double solid_angle = 0.;
    // ARAPUCAS/Bars (rectangle)
    if (opDet.type == 0) {
      // get scintillation point coordinates relative to arapuca window centre
      geo::Vector_t const abs_relative{
        std::abs(relative.X()), std::abs(relative.Y()), std::abs(relative.Z())};
      solid_angle = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_relative, opDet.orientation);
    }
    // PMTs (dome)
    else if (opDet.type == 1) {
      solid_angle = Omega_Dome_Model(distance, theta);
    }
    // PMTs (disk)
    else if (opDet.type == 2) {
      const double zy_offset = std::sqrt(relative.Y() * relative.Y() + relative.Z() * relative.Z());
      const double x_distance = std::abs(relative.X());
      solid_angle = Disk_SolidAngle(zy_offset, x_distance, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk" << std::endl;
    }

    // calculate number of photons hits by geometric acceptance for fast and slow components
    // accounting for solid angle and LAr absorbtion length
    double hits_geo_fast =
      std::exp(-1. * distance / fL_abs_vuv) * (solid_angle / (4 * CLHEP::pi)) * NumFast;
    double hits_geo_slow =
      std::exp(-1. * distance / fL_abs_vuv) * (solid_angle / (4 * CLHEP::pi)) * NumSlow;

    // apply Gaisser-Hillas correction for Rayleigh scattering distance
    // and angular dependence offset angle bin
    const size_t j = (theta / fdelta_angulo_vuv);

    // determine GH parameters, accounting for border effects
    // radial distance from centre of detector (Y-Z)
    double r = std::hypot(ScintPoint.Y() - fcathode_centre[1], ScintPoint.Z() - fcathode_centre[2]);

    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0; double s2 = 0; double s3 = 0;
    // flat PDs
    if ((opDet.type == 0 || opDet.type == 2) && (fIsFlatPDCorr || fIsFlatPDCorrLat)){
      if (opDet.orientation == 1 && fIsFlatPDCorrLat) { // laterals, alternate parameterisation method
        // distance to anode plane
        double d_anode = std::abs(fanode_centre[0] - ScintPoint.X()); 
        
        // build arrays for interpolation
        int n_distances = fGH_distances_anode.size();
        std::vector<double> p1, p2, p3, p4;
        p1.reserve(n_distances); p2.reserve(n_distances); p3.reserve(n_distances); p4.reserve(n_distances);
        for (int i = 0; i < n_distances; i++) {
          p1.push_back(fGHvuvpars_flat_lateral[0][i][j]);
          p2.push_back(fGHvuvpars_flat_lateral[1][i][j]);
          p3.push_back(fGHvuvpars_flat_lateral[2][i][j]);
          p4.push_back(fGHvuvpars_flat_lateral[3][i][j]);
        }

        // interpolate in distance to anode
        pars_ini[0] = interpolate(fGH_distances_anode, p1, d_anode, false);
        pars_ini[1] = interpolate(fGH_distances_anode, p2, d_anode, false);
        pars_ini[2] = interpolate(fGH_distances_anode, p3, d_anode, false);
        pars_ini[3] = interpolate(fGH_distances_anode, p4, d_anode, false);
      
      }
      else if (opDet.orientation == 0 && fIsFlatPDCorr) { // cathode/anode, default parameterisation method
        pars_ini[0] = fGHvuvpars_flat[0][j];
        pars_ini[1] = fGHvuvpars_flat[1][j];
        pars_ini[2] = fGHvuvpars_flat[2][j];
        pars_ini[3] = fGHvuvpars_flat[3][j];
        s1 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[0], theta, true);
        s2 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[1], theta, true);
        s3 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[2], theta, true);
      }
      else std::cout << "Error: corrections for chosen optical detector type missing." << std::endl;
    }
    // dome PDs
    else if (opDet.type == 1 && fIsDomePDCorr) {
      pars_ini[0] = fGHvuvpars_dome[0][j];
      pars_ini[1] = fGHvuvpars_dome[1][j];
      pars_ini[2] = fGHvuvpars_dome[2][j];
      pars_ini[3] = fGHvuvpars_dome[3][j];
      s1 = interpolate( fborder_corr_angulo_dome, fborder_corr_dome[0], theta, true);
      s2 = interpolate( fborder_corr_angulo_dome, fborder_corr_dome[1], theta, true);
      s3 = interpolate( fborder_corr_angulo_dome, fborder_corr_dome[2], theta, true);
    }
    else std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk. Or corrections for chosen optical detector type missing." << std::endl;

    // add border correction to parameters
    pars_ini[0] = pars_ini[0] + s1 * r;
    pars_ini[1] = pars_ini[1] + s2 * r;
    pars_ini[2] = pars_ini[2] + s3 * r;
    pars_ini[3] = pars_ini[3];

    // calculate correction
    double GH_correction = Gaisser_Hillas(distance, pars_ini);

    // apply field cage transparency factor 
    if (fApplyFieldCageTransparency) {
      if (opDet.orientation == 1) GH_correction = GH_correction * fFieldCageTransparencyLateral;
      else if (opDet.orientation == 0) GH_correction = GH_correction * fFieldCageTransparencyCathode;
    }

    // calculate number photons for fast and slow componenets
    DetThis[0] = fRandPoissPhot->fire(GH_correction * hits_geo_fast / cosine);
    DetThis[1] = fRandPoissPhot->fire(GH_correction * hits_geo_slow / cosine);

  }

  //......................................................................
  // VIS hits semi-analytic model calculation
  void
  PDFastSimPAR::detectedReflecHits(std::map<size_t, int>& ReflDetectedNumFast,
                                   std::map<size_t, int>& ReflDetectedNumSlow,
                                   const double NumFast,
                                   const double NumSlow,
                                   geo::Point_t const& ScintPoint,
                                   bool AnodeMode)
  {
    // 1). calculate total number of hits of VUV photons on
    // reflective foils via solid angle + Gaisser-Hillas
    // corrections:  

    // get scintpoint coords relative to centre of cathode plane and set plane dimensions
    geo::Vector_t ScintPoint_relative;
    Dims plane_dimensions;
    double plane_depth;
    if (AnodeMode) {
      plane_dimensions = fanode_plane;
      plane_depth = fanode_plane_depth;
      ScintPoint_relative.SetCoordinates(std::abs(ScintPoint.X() - fanode_plane_depth), std::abs(ScintPoint.Y() - fanode_centre[1]), std::abs(ScintPoint.Z() - fanode_centre[2]));
    }
    else {
      plane_dimensions = fcathode_plane;
      plane_depth = ScintPoint.X() < 0. ? -fplane_depth : fplane_depth;
      ScintPoint_relative.SetCoordinates(std::abs(ScintPoint.X() - plane_depth), std::abs(ScintPoint.Y() - fcathode_centre[1]), std::abs(ScintPoint.Z() - fcathode_centre[2]));
    }

    // calculate solid angle of cathode from the scintillation point, orientation always = 0 (cathode)
    double solid_angle_cathode = Rectangle_SolidAngle(plane_dimensions, ScintPoint_relative, 0);

    // calculate distance and angle between ScintPoint and hotspot
    // vast majority of hits in hotspot region directly infront of scintpoint,
    // therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
    double distance_cathode = std::abs(plane_depth - ScintPoint.X());
    // calculate hits on cathode plane via geometric acceptance
    double cathode_hits_geo_fast = std::exp(-1. * distance_cathode / fL_abs_vuv) *
                              (solid_angle_cathode / (4. * CLHEP::pi)) * NumFast;
    double cathode_hits_geo_slow = std::exp(-1. * distance_cathode / fL_abs_vuv) *
                              (solid_angle_cathode / (4. * CLHEP::pi)) * NumSlow;

    // determine Gaisser-Hillas correction including border effects
    // use flat correction
    double r = std::hypot(ScintPoint.Y() - fcathode_centre[1], ScintPoint.Z() - fcathode_centre[2]);
    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0; double s2 = 0; double s3 = 0;
    if(fIsFlatPDCorr) {
      pars_ini[0] = fGHvuvpars_flat[0][0];
      pars_ini[1] = fGHvuvpars_flat[1][0];
      pars_ini[2] = fGHvuvpars_flat[2][0];
      pars_ini[3] = fGHvuvpars_flat[3][0];
      s1 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[0], 0, true);
      s2 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[1], 0, true);
      s3 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[2], 0, true);
    }
    else std::cout << "Error: flat optical detector VUV correction required for reflected semi-analytic hits." << std::endl;

    // add border correction
    pars_ini[0] = pars_ini[0] + s1 * r;
    pars_ini[1] = pars_ini[1] + s2 * r;
    pars_ini[2] = pars_ini[2] + s3 * r;
    pars_ini[3] = pars_ini[3];


    // calculate corrected number of hits
    double GH_correction = Gaisser_Hillas(distance_cathode, pars_ini);
    const double cathode_hits_rec_fast = GH_correction * cathode_hits_geo_fast;
    const double cathode_hits_rec_slow = GH_correction * cathode_hits_geo_slow;

    // detemine hits on each PD
    const geo::Point_t hotspot = {plane_depth, ScintPoint.Y(), ScintPoint.Z()};
    for (size_t const OpDet : util::counter(nOpDets)) {
      if (!isOpDetInSameTPC(ScintPoint, fOpDetCenter[OpDet])) continue;

      // set detector struct for solid angle function
      const  PDFastSimPAR::OpticalDetector op{
        fOpDetHeight[OpDet], fOpDetLength[OpDet],
        fOpDetCenter[OpDet], fOpDetType[OpDet], fOpDetOrientation[OpDet]};

      std::vector<int> ReflDetThis(2, 0);
      VISHits(ScintPoint, op, cathode_hits_rec_fast, cathode_hits_rec_slow, hotspot, ReflDetThis, AnodeMode);

      ReflDetectedNumFast[OpDet] = ReflDetThis[0];
      ReflDetectedNumSlow[OpDet] = ReflDetThis[1];
    }
  }

  void
  PDFastSimPAR::VISHits(geo::Point_t const& ScintPoint,
                        OpticalDetector const& opDet,
                        const double cathode_hits_rec_fast,
                        const double cathode_hits_rec_slow,
                        geo::Point_t const& hotspot,
                        std::vector<int> &ReflDetThis,
                        bool AnodeMode)
  {

    // set correct plane_depth
    double plane_depth;
    if (AnodeMode) plane_depth = fanode_plane_depth;
    else plane_depth = ScintPoint.X() < 0. ? -fplane_depth : fplane_depth;

    // calculate number of these hits which reach the optical
    // detector from the hotspot using solid angle:

    geo::Vector_t const emission_relative = hotspot - opDet.OpDetPoint;

    // calculate distances and angles for application of corrections
    // distance from hotspot to optical detector
    const double distance_vis = emission_relative.R();
    //  angle between hotspot and optical detector
    double cosine_vis;
    if (opDet.orientation == 1) { // lateral
      cosine_vis = std::abs(emission_relative.Y()) / distance_vis;
    }
    else { // anode/cathode (default)
      cosine_vis = std::abs(emission_relative.X()) / distance_vis;
    }
    // const double theta_vis = std::acos(cosine_vis) * 180. / CLHEP::pi;
    const double theta_vis = fast_acos(cosine_vis) * 180. / CLHEP::pi;

    // calculate solid angle of optical channel
    double solid_angle_detector = 0.;
    // ARAPUCAS/Bars (rectangle)
    if (opDet.type == 0) {
      // get hotspot coordinates relative to opDet
      geo::Vector_t const abs_emission_relative{std::abs(emission_relative.X()),
                                                std::abs(emission_relative.Y()),
                                                std::abs(emission_relative.Z())};
      solid_angle_detector = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_emission_relative, opDet.orientation);
    }
    // PMTS (dome)
    else if (opDet.type == 1) {
      solid_angle_detector = Omega_Dome_Model(distance_vis, theta_vis);
    }
    // PMTs (disk)
    else if (opDet.type == 2) {
      const double zy_offset = std::sqrt(emission_relative.Y() * emission_relative.Y() +
                                         emission_relative.Z() * emission_relative.Z());
      const double x_distance = std::abs(emission_relative.X());
      solid_angle_detector = Disk_SolidAngle(zy_offset, x_distance, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk" << std::endl;
    }

    // calculate number of hits via geometeric acceptance
    double hits_geo_fast = (solid_angle_detector / (2. * CLHEP::pi)) *
                      cathode_hits_rec_fast; // 2*pi due to presence of reflective foils
    double hits_geo_slow = (solid_angle_detector / (2. * CLHEP::pi)) *
                      cathode_hits_rec_slow; // 2*pi due to presence of reflective foils

    // determine correction factor, depending on PD type
    const size_t k = (theta_vis / fdelta_angulo_vis);         // off-set angle bin
    double r = std::hypot(ScintPoint.Y() - fcathode_centre[1], ScintPoint.Z() - fcathode_centre[2]);
    double d_c = std::abs(ScintPoint.X() - plane_depth);       // distance to cathode
    double border_correction = 0;
    // flat PDs
    if ((opDet.type == 0 || opDet.type == 2) && (fIsFlatPDCorr || fIsFlatPDCorrLat)){

      // select correct parameter set depending on PD orientation
      std::vector<double> vis_distances_x;
      std::vector<double> vis_distances_r;
      std::vector<std::vector<std::vector<double>>> vispars;
      if (opDet.orientation == 1 && fIsFlatPDCorrLat) { // laterals
        vis_distances_x = fvis_distances_x_flat_lateral;
        vis_distances_r = fvis_distances_r_flat_lateral;
        vispars = fvispars_flat_lateral;
      }
      else if (opDet.orientation == 0 && fIsFlatPDCorr) { // cathode/anode
        vis_distances_x = fvis_distances_x_flat;
        vis_distances_r = fvis_distances_r_flat;
        vispars = fvispars_flat;
      }
      else std::cout << "Error: corrections for chosen optical detector type missing." << std::endl;
      
      // interpolate in d_c for each r bin
      const size_t nbins_r = vispars[k].size();
      std::vector<double> interp_vals(nbins_r, 0.0);
      {
        size_t idx = 0;
        size_t size = vis_distances_x.size();
        if (d_c >= vis_distances_x[size - 2])
          idx = size - 2;
        else {
          while (d_c > vis_distances_x[idx + 1])
            idx++;
        }
        for (size_t i = 0; i < nbins_r; ++i) {
          interp_vals[i] = interpolate(vis_distances_x,
                                       vispars[k][i],
                                       d_c,
                                       false,
                                       idx);
        }
      }
      // interpolate in r
      border_correction = interpolate(vis_distances_r, interp_vals, r, false); 
    }
    // dome PDs
    else if (opDet.type == 1 && fIsDomePDCorr) {
      // interpolate in d_c for each r bin
      const size_t nbins_r = fvispars_dome[k].size();
      std::vector<double> interp_vals(nbins_r, 0.0);
      {
        size_t idx = 0;
        size_t size = fvis_distances_x_dome.size();
        if (d_c >= fvis_distances_x_dome[size - 2])
          idx = size - 2;
        else {
          while (d_c > fvis_distances_x_dome[idx + 1])
            idx++;
        }
        for (size_t i = 0; i < nbins_r; ++i) {
          interp_vals[i] = interpolate(fvis_distances_x_dome,
                                       fvispars_dome[k][i],
                                       d_c,
                                       false,
                                       idx);
        }
      }
      // interpolate in r
      border_correction = interpolate(fvis_distances_r_dome, interp_vals, r, false);
    }
    else {
     std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk. Or corrections for chosen optical detector type missing." << std::endl;
    }

    // apply anode reflectivity factor
    if (AnodeMode) border_correction = border_correction * fAnodeReflectivity;

    // apply field cage transparency factor 
    if (fApplyFieldCageTransparency) {
      if (opDet.orientation == 1) border_correction = border_correction * fFieldCageTransparencyLateral;
      else if (opDet.orientation == 0) border_correction = border_correction * fFieldCageTransparencyCathode;
    }

    ReflDetThis[0] = fRandPoissPhot->fire(border_correction * hits_geo_fast / cosine_vis);
    ReflDetThis[1] = fRandPoissPhot->fire(border_correction * hits_geo_slow / cosine_vis);
  }


  bool
  PDFastSimPAR::isOpDetInSameTPC(geo::Point_t const& ScintPoint,
                                 geo::Point_t const& OpDetPoint) const
  {
    // check optical channel is in same TPC as scintillation light, if not return 0 hits
    // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero, otherwise return 0 hits
    if (((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) &&
        std::abs(OpDetPoint.X()) > 10. && fNTPC == 2) { // TODO: unhardcode
      return false;
    }
    return true;
  }

  //......................................................................
  void
  PDFastSimPAR::propagationTime(std::vector<double>& arrival_time_dist,
                                geo::Point_t const& x0,
                                const size_t OpChannel,
                                bool Reflected)
  {
    if (fIncludePropTime && !fGeoPropTimeOnly) {
      // Get VUV photons arrival time distribution from the parametrization
      geo::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
      if (!Reflected) {
        double distance = std::hypot(x0.X() - opDetCenter.X(), x0.Y() - opDetCenter.Y(), x0.Z() - opDetCenter.Z());
        double cosine; 
        if (fOpDetOrientation[OpChannel] == 1) cosine = std::abs(x0.Y() - opDetCenter.Y()) / distance;
        else cosine = std::abs(x0.X() - opDetCenter.X()) / distance;
      
        double theta = fast_acos(cosine)*180./CLHEP::pi;
        int angle_bin = theta/fangle_bin_timing_vuv;
        getVUVTimes(arrival_time_dist, distance, angle_bin); // in ns
      }
      else {
        getVISTimes(arrival_time_dist, geo::vect::toTVector3(x0),
                    geo::vect::toTVector3(opDetCenter)); // in ns
      }
    }
    else if (fIncludePropTime && fGeoPropTimeOnly && !Reflected) {
      // Get VUV photons arrival time geometrically
      geo::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
      double distance = std::hypot(x0.X() - opDetCenter.X(), x0.Y() - opDetCenter.Y(), x0.Z() - opDetCenter.Z());
      getVUVTimesGeo(arrival_time_dist, distance); // in ns
    }
    else {
      throw cet::exception("PDFastSimPAR")
        << "Propagation time model not found.";
    }
  }

  //......................................................................
  // VUV arrival times calculation function
  void
  PDFastSimPAR::getVUVTimes(std::vector<double>& arrivalTimes, const double distance, const size_t angle_bin)
  {
    if (distance < fmin_d) {
      // times are fixed shift i.e. direct path only
      double t_prop_correction = distance / fvuv_vgroup_mean;
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = t_prop_correction;
      }
    }
    else { // distance >= 25cm
      // determine nearest parameterisation in discretisation
      int index = std::round((distance - fmin_d) / fstep_size);
      // check whether required parameterisation has been generated, generating if not
       if (VUV_timing[angle_bin][index].GetNdim() == 0) { generateParam(index, angle_bin); }
      // randomly sample parameterisation for each photon
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = VUV_timing[angle_bin][index].GetRandom(VUV_min[angle_bin][index], VUV_max[angle_bin][index]);
      }
    }
  }

  //......................................................................
  // VUV arrival times calculation function - pure geometric approximation for use in Xenon doped scenarios
  void 
  PDFastSimPAR::getVUVTimesGeo(std::vector<double>& arrivalTimes, const double distance) 
  {
    // times are fixed shift i.e. direct path only
    double t_prop_correction = distance / fvuv_vgroup_mean;
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = t_prop_correction;
    }
  }

  //......................................................................
  // VIS arrival times calculation functions
  void
  PDFastSimPAR::getVISTimes(std::vector<double>& arrivalTimes,
                            const TVector3 &ScintPoint,
                            const TVector3 &OpDetPoint)
  {
    // *************************************************************************************************
    //     Calculation of earliest arrival times and corresponding unsmeared
    //     distribution
    // *************************************************************************************************

    // set plane_depth for correct TPC:
    double plane_depth;
    if (ScintPoint[0] < 0) { plane_depth = -fplane_depth; }
    else {
      plane_depth = fplane_depth;
    }

    // calculate point of reflection for shortest path
    TVector3 bounce_point(plane_depth,ScintPoint[1],ScintPoint[2]);

    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point - ScintPoint).Mag();
    double Visdist = (OpDetPoint - bounce_point).Mag();

    // calculate times taken by VUV part of path
    int angle_bin_vuv = 0; // on-axis by definition
    getVUVTimes(arrivalTimes, VUVdist, angle_bin_vuv);

    // sum parts to get total transport times times
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] += Visdist / fvis_vmean;
    }

    // *************************************************************************************************
    //      Smearing of arrival time distribution
    // *************************************************************************************************
    // calculate fastest time possible
    // vis part
    double vis_time = Visdist / fvis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < fmin_d) {
      vuv_time = VUVdist / fvuv_vgroup_max;
    }
    else {
      // find index of required parameterisation
      const size_t index = std::round((VUVdist - fmin_d) / fstep_size);
      // find shortest time
      vuv_time = VUV_min[angle_bin_vuv][index];
    }
    // sum
    double fastest_time = vis_time + vuv_time;

    // calculate angle theta between bound_point and optical detector
    double cosine_theta = std::abs(OpDetPoint[0] - bounce_point[0]) / Visdist;
    double theta = fast_acos(cosine_theta) * 180. / CLHEP::pi;

    // determine smearing parameters using interpolation of generated points:
    // 1). tau = exponential smearing factor, varies with distance and angle
    // 2). cutoff = largest smeared time allowed, preventing excessively large
    //     times caused by exponential distance to cathode
    double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
    // angular bin
    size_t theta_bin = theta / fangle_bin_timing_vis;
    // radial distance from centre of TPC (y,z plane)
    double r = std::sqrt(std::pow(ScintPoint[1] - fcathode_centre[1], 2) + std::pow(ScintPoint[2] - fcathode_centre[2], 2));

    // cut-off and tau
    // cut-off
    // interpolate in d_c for each r bin
    std::vector<double> interp_vals(fcut_off_pars[theta_bin].size(), 0.0);
    for (size_t i = 0; i < fcut_off_pars[theta_bin].size(); i++){
        interp_vals[i] = interpolate(fdistances_refl, fcut_off_pars[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double cutoff = interpolate(fradial_distances_refl, interp_vals, r, true);

    // tau
    // interpolate in x for each r bin
    std::vector<double> interp_vals_tau(ftau_pars[theta_bin].size(), 0.0);
    for (size_t i = 0; i < ftau_pars[theta_bin].size(); i++){
        interp_vals_tau[i] = interpolate(fdistances_refl, ftau_pars[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double tau = interpolate(fradial_distances_refl, interp_vals_tau, r, true);

    // apply smearing:
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      double arrival_time_smeared;
      // if time is already greater than cutoff, do not apply smearing
      if (arrivalTimes[i] >= cutoff) { continue; }
      // otherwise smear
      else {
        unsigned int counter = 0;
        // loop until time generated is within cutoff limit
        // most are within single attempt, very few take more than two
        do {
          // don't attempt smearings too many times
          if (counter >= 10) {                      // TODO: unhardcode
            arrival_time_smeared = arrivalTimes[i]; // don't smear
            break;
          }
          else {
            // generate random number in appropriate range
            double x = CLHEP::RandFlat::shoot(&fScintTimeEngine, 0.5, 1.0);
            // apply the exponential smearing
            arrival_time_smeared =
              arrivalTimes[i] + (arrivalTimes[i] - fastest_time) * (std::pow(x, -tau) - 1);
          }
          counter++;
        } while (arrival_time_smeared > cutoff);
      }
      arrivalTimes[i] = arrival_time_smeared;
    }
  }

  //......................................................................
  void
  PDFastSimPAR::generateParam(const size_t index, const size_t angle_bin)
  {
    // get distance
    double distance_in_cm = (index * fstep_size) + fmin_d;

    // time range
    const double signal_t_range = 5000.; // TODO: unhardcode

    // parameterisation TF1
    TF1 fVUVTiming;

    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
    double t_direct_min = distance_in_cm / fvuv_vgroup_max;

    // Defining the model function(s) describing the photon transportation timing vs distance
    // Getting the landau parameters from the time parametrization
    std::array<double, 3> pars_landau;
    interpolate3(pars_landau,
                 fparameters[0][0],
                 fparameters[2][angle_bin],
                 fparameters[3][angle_bin],
                 fparameters[1][angle_bin],
                 distance_in_cm,
                 true);
    // Deciding which time model to use (depends on the distance)
    // defining useful times for the VUV arrival time shapes
    if (distance_in_cm >= finflexion_point_distance) {
      double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
      // Set model: Landau
      fVUVTiming = TF1("fVUVTiming", model_far, 0, signal_t_range, 4);
      fVUVTiming.SetParameters(pars_far);
    }
    else {
      // Set model: Landau + Exponential
      fVUVTiming = TF1("fVUVTiming", model_close, 0, signal_t_range, 7);
      // Exponential parameters
      double pars_expo[2];
      // Getting the exponential parameters from the time parametrization
      pars_expo[1] = interpolate(fparameters[4][0], fparameters[5][angle_bin], distance_in_cm, true);
      pars_expo[0] = interpolate(fparameters[4][0], fparameters[6][angle_bin], distance_in_cm, true);
      pars_expo[0] *= pars_landau[2];
      pars_expo[0] = std::log(pars_expo[0]);
      // this is to find the intersection point between the two functions:
      TF1 fint = TF1("fint", finter_d, pars_landau[0], 4 * t_direct_mean, 5);
      double parsInt[5] = {
        pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
      fint.SetParameters(parsInt);
      double t_int = fint.GetMinimumX();
      double minVal = fint.Eval(t_int);
      // the functions must intersect - output warning if they don't
      if (minVal > 0.015) {
        std::cout << "WARNING: Parametrization of VUV light discontinuous for distance = "
                  << distance_in_cm << std::endl;
        std::cout << "WARNING: This shouldn't be happening " << std::endl;
      }
      double parsfinal[7] = {t_int,
                             pars_landau[0],
                             pars_landau[1],
                             pars_landau[2],
                             pars_expo[0],
                             pars_expo[1],
                             t_direct_min};
      fVUVTiming.SetParameters(parsfinal);
    }

    // set the number of points used to sample parameterisation
    // for shorter distances, peak is sharper so more sensitive sampling required
    int fsampling; // TODO: unhardcode
    if (distance_in_cm < 50) { fsampling = 10000; }
    else if (distance_in_cm < 100) {
      fsampling = 5000;
    }
    else {
      fsampling = 1000;
    }
    fVUVTiming.SetNpx(fsampling);

    // calculate max and min distance relevant to sample parameterisation
    // max
    const size_t nq_max = 1;
    double xq_max[nq_max];
    double yq_max[nq_max];
    xq_max[0] = 0.995; // include 99.5%
    fVUVTiming.GetQuantiles(nq_max, yq_max, xq_max);
    double max = yq_max[0];
    // min
    double min = t_direct_min;

    // store TF1 and min/max, this allows identical TF1 to be used every time sampling
    // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
    // all subsequent calls check if it has been generated previously and are ~100+ times quicker
    VUV_timing[angle_bin][index] = fVUVTiming;
    VUV_max[angle_bin][index] = max;
    VUV_min[angle_bin][index] = min;
  }

  //======================================================================
  //   Returns interpolated value at x from parallel arrays ( xData, yData )
  //   Assumes that xData has at least two elements, is sorted and is strictly
  //   monotonic increasing boolean argument extrapolate determines behaviour
  //   beyond ends of array (if needed)
  double
  PDFastSimPAR::interpolate(const std::vector<double>& xData,
                            const std::vector<double>& yData,
                            double x,
                            bool extrapolate,
                            size_t i)
  {
    if (i == 0) {
      size_t size = xData.size();
      if (x >= xData[size - 2]) { // special case: beyond right end
        i = size - 2;
      }
      else {
        while (x > xData[i + 1])
          i++;
      }
    }
    double xL = xData[i];
    double xR = xData[i + 1];
    double yL = yData[i];
    double yR = yData[i + 1]; // points on either side (unless beyond ends)
    if (!extrapolate) {       // if beyond ends of array and not extrapolating
      if (x < xL) return yL;
      if (x > xR) return yL;
    }
    const double dydx = (yR - yL) / (xR - xL); // gradient
    return yL + dydx * (x - xL);               // linear interpolation
  }

  void
  PDFastSimPAR::interpolate3(std::array<double, 3>& inter,
                             const std::vector<double>& xData,
                             const std::vector<double>& yData1,
                             const std::vector<double>& yData2,
                             const std::vector<double>& yData3,
                             double x,
                             bool extrapolate)
  {
    size_t size = xData.size();
    size_t i = 0;               // find left end of interval for interpolation
    if (x >= xData[size - 2]) { // special case: beyond right end
      i = size - 2;
    }
    else {
      while (x > xData[i + 1])
        i++;
    }
    double xL = xData[i];
    double xR = xData[i + 1]; // points on either side (unless beyond ends)
    double yL1 = yData1[i];
    double yR1 = yData1[i + 1];
    double yL2 = yData2[i];
    double yR2 = yData2[i + 1];
    double yL3 = yData3[i];
    double yR3 = yData3[i + 1];

    if (!extrapolate) { // if beyond ends of array and not extrapolating
      if (x < xL) {
        inter[0] = yL1;
        inter[1] = yL2;
        inter[2] = yL3;
        return;
      }
      if (x > xR) {
        inter[0] = yL1;
        inter[1] = yL2;
        inter[2] = yL3;
        return;
      }
    }
    const double m = (x - xL) / (xR - xL);
    inter[0] = m * (yR1 - yL1) + yL1;
    inter[1] = m * (yR2 - yL2) + yL2;
    inter[2] = m * (yR3 - yL3) + yL3;
  }

  //......................................................................
  // solid angle of circular aperture
  double
  PDFastSimPAR::Disk_SolidAngle(const double d, const double h, const double b)
  {
    if (b <= 0. || d < 0. || h <= 0.) return 0.;
    const double leg2 = (b + d) * (b + d);
    const double aa = std::sqrt(h * h / (h * h + leg2));
    if (isApproximatelyZero(d)) { return 2. * CLHEP::pi * (1. - aa); }
    double bb = 2. * std::sqrt(b * d / (h * h + leg2));
    double cc = 4. * b * d / leg2;

    if (isDefinitelyGreaterThan(d, b)) {
      try {
        return 2. * aa *
               (std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()) -
                boost::math::ellint_1(bb, noLDoublePromote()));
      }
      catch (std::domain_error& e) {
        if (isApproximatelyEqual(d, b, 1e-9)) {
          mf::LogWarning("PDFastSimPAR")
            << "Elliptic Integral in Disk_SolidAngle() given parameters "
               "outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("PDFastSimPAR")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters "
               "outside domain.\n"
            << "\nbb: " << bb << "\ncc: " << cc << "Exception message: " << e.what();
          return 0.;
        }
      }
    }
    if (isDefinitelyLessThan(d, b)) {
      try {
        return 2. * CLHEP::pi -
               2. * aa *
                 (boost::math::ellint_1(bb, noLDoublePromote()) +
                  std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()));
      }
      catch (std::domain_error& e) {
        if (isApproximatelyEqual(d, b, 1e-9)) {
          mf::LogWarning("PDFastSimPAR")
            << "Elliptic Integral in Disk_SolidAngle() given parameters "
               "outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("PDFastSimPAR")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters "
               "outside domain.\n"
            << "\nbb: " << bb << "\ncc: " << cc << "Exception message: " << e.what();
          return 0.;
        }
      }
    }
    if (isApproximatelyEqual(d, b)) {
      return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
    }
    return 0.;
  }

  //......................................................................
  // solid angle of rectangular aperture
  double
  PDFastSimPAR::Rectangle_SolidAngle(const double a, const double b, const double d)
  {
    double aa = a / (2. * d);
    double bb = b / (2. * d);
    double aux = (1. + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
    // return 4 * std::acos(std::sqrt(aux));
    return 4. * fast_acos(std::sqrt(aux));
  }

  double
  PDFastSimPAR::Rectangle_SolidAngle(Dims const& o, geo::Vector_t const& v, double OpDetOrientation)
  {
    // v is the position of the track segment with respect to
    // the center position of the arapuca window

    // solid angle calculation depends on orientation of PD, set correct distances to use
    double d1;
    double d2; 
    if (OpDetOrientation == 1) {
      // lateral PD, arapuca plane fixed in y direction
      d1 = std::abs(v.X());
      d2 = std::abs(v.Y()); 
    }
    else {
      // anode/cathode PD, arapuca plane fixed in x direction [default]
      d1 = std::abs(v.Y());
      d2 = std::abs(v.X());
    }
    // arapuca plane fixed in x direction
    if (isApproximatelyZero(d1) && isApproximatelyZero(v.Z())) {
      return Rectangle_SolidAngle(o.h, o.w, d2);
    }
    if (isDefinitelyGreaterThan(d1, o.h * .5) && isDefinitelyGreaterThan(std::abs(v.Z()), o.w * .5)) {
      double A = d1 - o.h * .5;
      double B = std::abs(v.Z()) - o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (A + o.h), 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * A, 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * (A + o.h), 2. * B, d2) +
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    if ((d1 <= o.h * .5) && (std::abs(v.Z()) <= o.w * .5)) {
      double A = -d1 + o.h * .5;
      double B = -std::abs(v.Z()) + o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (o.h - A), 2. * (o.w - B), d2) +
                          Rectangle_SolidAngle(2. * A, 2. * (o.w - B), d2) +
                          Rectangle_SolidAngle(2. * (o.h - A), 2. * B, d2) +
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    if (isDefinitelyGreaterThan(d1, o.h * .5) && (std::abs(v.Z()) <= o.w * .5)) {
      double A = d1 - o.h * .5;
      double B = -std::abs(v.Z()) + o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (A + o.h), 2. * (o.w - B), d2) -
                          Rectangle_SolidAngle(2. * A, 2. * (o.w - B), d2) +
                          Rectangle_SolidAngle(2. * (A + o.h), 2. * B, d2) -
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    if ((d1 <= o.h * .5) && isDefinitelyGreaterThan(std::abs(v.Z()), o.w * .5)) {
      double A = -d1 + o.h * .5;
      double B = std::abs(v.Z()) - o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (o.h - A), 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * (o.h - A), 2. * B, d2) +
                          Rectangle_SolidAngle(2. * A, 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    // error message if none of these cases, i.e. something has gone wrong!
    // std::cout << "Warning: invalid solid angle call." << std::endl;
    return 0.;
  }

  //......................................................................
  // solid angle of dome aperture
  double
  PDFastSimPAR::Omega_Dome_Model(const double distance, const double theta) const {
    // this function calculates the solid angle of a semi-sphere of radius b,
    // as a correction to the analytic formula of the on-axix solid angle,
    // as we move off-axis an angle theta. We have used 9-angular bins
    // with delta_theta width.

    // par0 = Radius correction close
    // par1 = Radius correction far
    // par2 = breaking distance betwween "close" and "far"

    double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221};
    double par1[9] = {0, 0, 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
    const double delta_theta = 10.;
    int j = int(theta/delta_theta);
    // PMT radius
    const double b = fradius; // cm
    // distance form which the model parameters break (empirical value)
    const double d_break = 5*b; //par2

    if(distance >= d_break) {
      double R_apparent_far = b - par1[j];
      return  (2*CLHEP::pi * (1 - std::sqrt(1 - std::pow(R_apparent_far/distance,2))));

    }
    else {
      double R_apparent_close = b - par0[j];
      return (2*CLHEP::pi * (1 - std::sqrt(1 - std::pow(R_apparent_close/distance,2))));
    }
  }


  // ---------------------------------------------------------------------------

} // namespace phot

DEFINE_ART_MODULE(phot::PDFastSimPAR)
