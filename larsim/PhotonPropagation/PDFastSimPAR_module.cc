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

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
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
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityTypes.h" // phot::MappedT0s_t
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"
#include "larsim/Simulation/LArG4Parameters.h"

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

  template <typename TReal>
  inline constexpr double
  dist(const TReal* x, const TReal* y, const unsigned int dimension)
  {
    double d = 0.;
    for (unsigned int p = 0; p < dimension; ++p) {
      d += (*(x + p) - *(y + p)) * (*(x + p) - *(y + p));
    }
    return std::sqrt(d);
  }
  template <typename TVector3>
  inline constexpr double
  dist(const std::array<double, 3> x, const TVector3 y, const unsigned int dimension, const unsigned int start)
  {
    double d = 0.;
    for (unsigned int p = start; p < dimension; ++p) {
      d += (x[p] - y[p]) * (x[p] - y[p]);
    }
    return std::sqrt(d);
  }
} // namespace

namespace phot {
  class PDFastSimPAR : public art::EDProducer {
  public:
    explicit PDFastSimPAR(fhicl::ParameterSet const&);
    void produce(art::Event&) override;

    void Initialization();

    void getVUVTimes(std::vector<double>& arrivalTimes, double distance_in_cm);
    void getVISTimes(std::vector<double>& arrivalTimes, TVector3 ScintPoint, TVector3 OpDetPoint);

    void generateParam(const size_t index);

    void AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                     std::map<size_t, int>& ChannelMap,
                     sim::OpDetBacktrackerRecord btr);

  private:
    // structure definition for solid angle of rectangle function
    struct Dims {
      double h, w; // height, width
    };

    struct OpticalDetector {
      double h; // height
      double w; // width
      geo::Point_t OpDetPoint;
      int type;
    };

    int VUVHits(const double Nphotons_created,   geo::Point_t const& ScintPoint, OpticalDetector const& opDet);
    int VISHits(const double Nphotons_created,   geo::Point_t const& ScintPoint, OpticalDetector const& opDet);

    std::vector<double> PropTime(geo::Point_t ScintPoint, int OpChannel, int NPhotons, bool is_refl); // const;

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
    double Rectangle_SolidAngle_Cal(const double a, const double b, const double d);
    double Rectangle_SolidAngle(Dims const& o, geo::Vector_t const& v);
    // solid angle of circular aperture calculation functions
    double Disk_SolidAngle(const double d, const double h, const double b);
    double Omega_Dome_Model(const double distance, const double theta);

    bool fDoSlowComponent;
    art::InputTag simTag;
    std::unique_ptr<ScintTime> fScintTime; // Tool to retrive timinig of scintillation
    CLHEP::HepRandomEngine& fPhotonEngine;
    std::unique_ptr<CLHEP::RandPoissonQ> fRandPoissPhot;
    CLHEP::HepRandomEngine& fScintTimeEngine;

    std::map<size_t, int> PDChannelToSOCMapDirect; // Where each OpChan is.
    std::map<size_t, int> PDChannelToSOCMapReflect; // Where each OpChan is.

    // For new VUV time parametrization
    double fstep_size, fmax_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance;
    std::vector<double> fparameters[9];
    // vector containing generated VUV timing parameterisations
    std::vector<TF1> VUV_timing;
    // vector containing min and max range VUV timing parameterisations are
    // sampled to
    std::vector<double> VUV_max;
    std::vector<double> VUV_min;

    // For new VIS time parameterisation
    double fvis_vmean, fn_LAr_vis, fn_LAr_vuv;
    std::vector<double> fdistances_refl;
    std::vector<std::vector<double>> fcut_off_pars;
    std::vector<std::vector<double>> ftau_pars;

    // For VUV semi-analytic hits
    // array of correction for the VUV Nhits estimation
    std::vector<std::vector<double>> fGHvuvpars;
    // To account for the border effects
    std::vector<double> fborder_corr;
    double fYactive_corner, fZactive_corner, fReference_to_corner;

    std::vector<geo::BoxBoundedGeo> const fActiveVolumes;
    // array of corrections for VIS Nhits estimation
//    std::vector<std::vector<double>> fvispars;
    std::vector<double> fvis_border_distances_x;
    std::vector<double> fvis_border_distances_r;
    std::vector<std::vector<std::vector<double>>> fvis_border_correction;
    bool fApplyVisBorderCorrection;
    std::string fVisBorderCorrectionType;

    double fplane_depth;
    geo::Point_t fCathode_centre;

    // Optical detector properties for semi-analytic hits
    double fradius;
    Dims fcathode_plane;
    
    int fdelta_angle;
    int fL_abs_vuv;
    std::vector<geo::Point_t> fOpDetCenter;
    std::vector<int> fOpDetType;
    std::vector<double> fOpDetLength;
    std::vector<double> fOpDetHeight;


    // For VUV semi-analytic hits
    double fdelta_angulo_vuv;
    // flat PDs
    bool fIsFlatPDCorr;
    std::vector<std::vector<double>> fGHvuvpars_flat;
    std::vector<double> fborder_corr_angulo_flat;
    std::vector<std::vector<double>> fborder_corr_flat;
    // dome PDs
    bool fIsDomePDCorr;
    std::vector<std::vector<double>> fGHvuvpars_dome;
    std::vector<double> fborder_corr_angulo_dome;
    std::vector<std::vector<double>> fborder_corr_dome;

    // For VIS semi-analytic hits
    bool fStoreReflected;
    // correction parameters for VIS Nhits estimation
    double fdelta_angulo_vis;
    // flat PDs
    std::vector<double> fvis_distances_x_flat;
    std::vector<double> fvis_distances_r_flat;
    std::vector<std::vector<std::vector<double>>> fvispars_flat;
    // dome PDs
    std::vector<double> fvis_distances_x_dome;
    std::vector<double> fvis_distances_r_dome;
    std::vector<std::vector<std::vector<double>>> fvispars_dome;
    
    double fcathode_zdimension, fcathode_ydimension;

    // Photon visibility service instance.
    PhotonVisibilityService const* const fPVS;

    MappedFunctions_t ParPropTimeTF1;

    /// Whether photon propagation is performed only from active volumes
    bool const fOnlyActiveVolume = false;
    /// Allows running even if light on cryostats `C:1` and higher is not supported.
    /// Currently hard coded "no"
    bool const fOnlyOneCryostat = false;
    /// Whether the cathodes are fully opaque; currently hard coded "no".
    bool const fOpaqueCathode = false;

    bool isOpDetInSameTPC(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) const;
    bool isScintInActiveVolume(geo::Point_t const& ScintPoint);

    static std::vector<geo::BoxBoundedGeo> extractActiveVolumes(geo::GeometryCore const& geom);
  };

  //......................................................................
  PDFastSimPAR::PDFastSimPAR(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fDoSlowComponent{pset.get<bool>("DoSlowComponent")}
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fScintTime{art::make_tool<ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService> {}
                      ->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,
                                                                                 "HepJamesRandom",
                                                                                 "scinttime",
                                                                                 pset,
                                                                                 "SeedScintTime"))
    , fActiveVolumes{extractActiveVolumes(*(lar::providerFrom<geo::Geometry>()))}
    , fPVS(art::ServiceHandle<PhotonVisibilityService const>().get())
    {
        std::cout << "PDFastSimPAR Module Construct" << std::endl;
        
        Initialization();
        art::ServiceHandle<sim::LArG4Parameters const> lgp;        
        if (lgp->UseLitePhotons())
        {
            std::cout << "Use Lite Photon." << std::endl;
            produces< std::vector<sim::SimPhotonsLite> >();
            produces< std::vector<sim::OpDetBacktrackerRecord> >();
            
            if(fPVS->StoreReflected())
            {
                std::cout << "Store Reflected Photons" << std::endl;
                produces< std::vector<sim::SimPhotonsLite> >("Reflected");
                produces< std::vector<sim::OpDetBacktrackerRecord> >("Reflected");     
            }            
        }
        else
        {
            std::cout << "Use Sim Photon." << std::endl;
            produces< std::vector<sim::SimPhotons> >();
            if(fPVS->StoreReflected())
            {
                std::cout << "Store Reflected Photons" << std::endl;            
                produces< std::vector<sim::SimPhotons> >("Reflected");     
            }            
        }        
    }
  //......................................................................
  void
  PDFastSimPAR::produce(art::Event& event)
  {
      std::cout << "PDFastSimPAR Module Producer" << std::endl;
      
      art::ServiceHandle<sim::LArG4Parameters const> lgp;
      
      auto const nOpChannels = fPVS->NOpChannels();
      
      std::unique_ptr< std::vector< sim::SimPhotons > >             phot   (new std::vector<sim::SimPhotons>);
      std::unique_ptr< std::vector< sim::SimPhotonsLite > >         phlit  (new std::vector<sim::SimPhotonsLite>);
      std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr  (new std::vector<sim::OpDetBacktrackerRecord>);
      
      std::unique_ptr< std::vector< sim::SimPhotons > >             phot_ref   (new std::vector<sim::SimPhotons>);
      std::unique_ptr< std::vector< sim::SimPhotonsLite > >         phlit_ref  (new std::vector<sim::SimPhotonsLite>);
      std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr_ref  (new std::vector<sim::OpDetBacktrackerRecord>);
      
      auto& dir_photcol(*phot);
      auto& ref_photcol(*phot_ref);
      auto& dir_phlitcol(*phlit);
      auto& ref_phlitcol(*phlit_ref);
      dir_photcol.resize(nOpChannels);
      ref_photcol.resize(nOpChannels);
      dir_phlitcol.resize(nOpChannels);
      ref_phlitcol.resize(nOpChannels);
      for (unsigned int i = 0; i < nOpChannels; i ++)
      {
          dir_photcol[i].fOpChannel  = i;
          ref_photcol[i].fOpChannel  = i;
          dir_phlitcol[i].OpChannel  = i;
          ref_phlitcol[i].OpChannel  = i;
      }

      art::Handle<std::vector<sim::SimEnergyDeposit>> edepHandle;
      if (!event.getByLabel(simTag, edepHandle))
      {
          std::cout << "PDFastSimPAR Module Cannot getByLabel: " << simTag << std::endl;
          return;
      }

      auto const& edeps = edepHandle;
      
      int num_points    = 0;
      for (auto const& edepi : *edeps) 
      {
          num_points++;
    
          int trackID = edepi.TrackID();
          double nphot = edepi.NumPhotons();
          double edeposit = edepi.Energy() / nphot;
          double pos[3] = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
          geo::Point_t const ScintPoint = {pos[0], pos[1], pos[2]};
          if (fOnlyActiveVolume && !isScintInActiveVolume(ScintPoint))
          {
              continue;
          }
          
          double nphot_fast = edepi.NumFPhotons();
          double nphot_slow = edepi.NumSPhotons();
          
          for(unsigned int channel = 0; channel< nOpChannels; channel++)
          {
              if(fOpaqueCathode && !isOpDetInSameTPC(ScintPoint, fOpDetCenter.at(channel)))
              {
                  continue;
              }
              const PDFastSimPAR::OpticalDetector op{fOpDetHeight.at(channel),
                                                     fOpDetLength.at(channel),
                                                     fOpDetCenter.at(channel),
                                                     fOpDetType.at(channel)};                          
              if (lgp->UseLitePhotons())
              {
                  sim::OpDetBacktrackerRecord tmpbtr(channel);
                  if (nphot_fast > 0)
                  {
                      auto n = VUVHits(nphot_fast, ScintPoint, op);
                      std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, false);
                      for (long i = 0; i < n; ++i)
                      {
                          //calculates the time at which the photon was produced
                          fScintTime->GenScintTime(true, fScintTimeEngine);
                          auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                          ++ dir_phlitcol[channel].DetectedPhotons[time];
                          tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                      }
                  }
                  
                  if ((nphot_slow > 0) && fDoSlowComponent)
                  {
                      auto n = VUVHits(nphot_slow, ScintPoint, op);
                      std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, false);
                      for (long i = 0; i < n; ++i)
                      {
                          fScintTime->GenScintTime(false, fScintTimeEngine);
                          auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                          ++ dir_phlitcol[channel].DetectedPhotons[time];
                          tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                      }
                  }
                  
                  AddOpDetBTR(*opbtr, PDChannelToSOCMapDirect, tmpbtr);
                  
                  if(fPVS->StoreReflected())
                  {
                      sim::OpDetBacktrackerRecord tmpbtr_ref(channel); 
                      if (nphot_fast > 0)
                      {
                          auto n = VISHits(nphot_fast, ScintPoint, op);
                          std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, true);
                          for (long i = 0; i < n; ++i)
                          {
                              //calculates the time at which the photon was produced
                              fScintTime->GenScintTime(true, fScintTimeEngine);
                              auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                              ++ ref_phlitcol[channel].DetectedPhotons[time];
                              tmpbtr_ref.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                          }
                      }
                      
                      if ((nphot_slow > 0) && fDoSlowComponent)
                      {
                          auto n = VISHits(nphot_slow, ScintPoint, op);
                          std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, true);
                          for (long i = 0; i < n; ++i)
                          {
                              fScintTime->GenScintTime(false, fScintTimeEngine);
                              auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                              ++ ref_phlitcol[channel].DetectedPhotons[time];
                              tmpbtr_ref.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                          }
                      }
                      
                      AddOpDetBTR(*opbtr_ref, PDChannelToSOCMapReflect, tmpbtr_ref);
                  }
              }
              else
              {
                  sim::OnePhoton photon;
                  photon.SetInSD         = false;
                  photon.InitialPosition = edepi.End();
                  photon.Energy          = 9.7e-6;
                  
                  if (nphot_fast > 0)
                  {
                      auto n = VUVHits(nphot_fast, ScintPoint, op);
                      std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, false);
                      for (long i = 0; i < n; ++i)
                      {
                          fScintTime->GenScintTime(true, fScintTimeEngine);
                          auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                          photon.Time = time;
                          dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
                      }
                  }
                  
                  if ((nphot_slow > 0) && fDoSlowComponent)
                  {
                      auto n = VUVHits(nphot_slow, ScintPoint, op);
                      std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, false);
                      for (long i = 0; i < n; ++i)
                      {
                          fScintTime->GenScintTime(false, fScintTimeEngine);
                          auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                          photon.Time = time;
                          dir_photcol[channel].insert(dir_photcol[channel].end(), 1, photon);
                      }
                  }
                  
                  if(fPVS->StoreReflected())
                  {
                      sim::OnePhoton photon;
                      photon.SetInSD         = false;
                      photon.InitialPosition = edepi.End();
                      photon.Energy          = 2.7e-6;
                      
                      if (nphot_fast > 0)
                      {
                          auto n = VISHits(nphot_fast, ScintPoint, op);
                          std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, true);
                          for (long i = 0; i < n; ++i)
                          {
                              fScintTime->GenScintTime(true, fScintTimeEngine);
                              auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                              photon.Time = time;
                              ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                          }
                      }
                      
                      if ((nphot_slow > 0) && fDoSlowComponent)
                      {
                          auto n = VISHits(nphot_slow, ScintPoint, op);
                          std::vector<double> arrivaltime = PropTime(ScintPoint, channel, n, true);
                          for (long i = 0; i < n; ++i)
                          {
                              fScintTime->GenScintTime(false, fScintTimeEngine);
                              auto time   = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime() + arrivaltime[i]);
                              photon.Time = time;
                              ref_photcol[channel].insert(ref_photcol[channel].end(), 1, photon);
                          }
                      }
                  }
              }
          }
      }    

      PDChannelToSOCMapDirect.clear();
      PDChannelToSOCMapReflect.clear();
      
      if (lgp->UseLitePhotons())
      {
          event.put(move(phlit));
          event.put(move(opbtr));
          if (fPVS->StoreReflected())
          {
              event.put(move(phlit_ref), "Reflected");
              event.put(move(opbtr_ref), "Reflected");
          }
      }
      else
      {
          event.put(move(phot));
          if (fPVS->StoreReflected())
          {
              event.put(move(phot_ref), "Reflected");
          }
      }
      
      return;
  }

  //......................................................................
  void
  PDFastSimPAR::AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                            std::map<size_t, int>& ChannelMap,
                            sim::OpDetBacktrackerRecord btr)
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

    assert(fPVS);
    fRandPoissPhot = std::make_unique<CLHEP::RandPoissonQ>(fPhotonEngine);
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

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

    // Load Gaisser-Hillas corrections for VUV semi-analytic hits
    std::cout << "Loading the GH corrections" << std::endl;
    fPVS->LoadVUVSemiAnalyticProperties(fIsFlatPDCorr,
                                        fIsDomePDCorr,
                                        fdelta_angulo_vuv,
                                        fradius
                                      );
    if (!fIsFlatPDCorr && !fIsDomePDCorr) {
      throw cet::exception("OpFastScintillation")
          << "Both isFlatPDCorr and isDomePDCorr parameters are false, at least one type of parameterisation is required for the semi-analytic light simulation." << "\n";
    }
    if (fIsFlatPDCorr) {
      fPVS->LoadGHFlat(fGHvuvpars_flat,
                       fborder_corr_angulo_flat,
                       fborder_corr_flat
                      );
    }
    if (fIsDomePDCorr) {
      fPVS->LoadGHDome(fGHvuvpars_dome,
                       fborder_corr_angulo_dome,
                       fborder_corr_dome
                      );
    }


    fCathode_centre = {geom.TPC(0, 0).GetCathodeCenter().X(),
                       fActiveVolumes[0].CenterY(),
                       fActiveVolumes[0].CenterZ()};
    mf::LogTrace("PDFastSimPAR") << "fCathode_centre: " << fCathode_centre << " cm";

    for (size_t const i : util::counter(fPVS->NOpChannels())) {
      geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpDet(i);
      fOpDetCenter.push_back(opDet.GetCenter());
      if (opDet.isBar()) {
        fOpDetType.push_back(0); // Arapucas
        fOpDetLength.push_back(opDet.Length());
        fOpDetHeight.push_back(opDet.Height());
      }
      else {
        fOpDetType.push_back(1); // PMTs
        //    std::cout<<"Radio: "<<geom.OpDetGeoFromOpDet(i).RMax()<<std::endl;
        fOpDetLength.push_back(-1);
        fOpDetHeight.push_back(-1);
      }
      // std::cout << "OpChannel: " << i << "  Optical_Detector_Type: " << fOpDetType.at(i)
      //           << "  APERTURE_height: " << opDet.Height()
      //           << "  APERTURE_width: "  << opDet.Length()
      //           << std::endl;
    }

    if (fPVS->IncludePropTime()) {
      std::cout << "Using parameterisation of timings." << std::endl;
      // VUV time parapetrization
      fPVS->LoadTimingsForVUVPar(fparameters,
                                 fstep_size,
                                 fmax_d,
                                 fvuv_vgroup_mean,
                                 fvuv_vgroup_max,
                                 finflexion_point_distance);

      // create vector of empty TF1s that will be replaces with the
      // parameterisations that are generated as they are required
      // default TF1() constructor gives function with 0 dimensions,
      // can then check numDim to qucikly see if a parameterisation
      // has been generated
      // TODO: warning: conversion from double to int
      int num_params = (fmax_d - 25) / fstep_size; // for d < 25cm, no parameterisaton, a delta
                                                   // function is used instead // TODO: unhardcode
      std::vector<TF1> VUV_timing_temp(num_params, TF1());
      VUV_timing = VUV_timing_temp;

      // initialise vectors to contain range parameterisations sampled
      // to in each case when using TF1->GetRandom(xmin,xmax), must be
      // in same range otherwise sampling is regenerated, this is the
      // slow part!
      std::vector<double> VUV_empty(num_params, 0);
      VUV_max = VUV_empty;
      VUV_min = VUV_empty;

      // VIS time parameterisation
      if (fPVS->StoreReflected())
      {
          fPVS->LoadTimingsForVISPar(fdistances_refl, fcut_off_pars, ftau_pars,fvis_vmean, fn_LAr_vis, fn_LAr_vuv);
          
          std::cout << "Loading visible light corrections" << std::endl;
          fPVS->LoadVisSemiAnalyticProperties(fdelta_angulo_vis,fradius);
          if (fIsFlatPDCorr)
          {
              fPVS->LoadVisParsFlat(fvis_distances_x_flat, fvis_distances_r_flat, fvispars_flat);
          }
          if (fIsDomePDCorr)
          {
              fPVS->LoadVisParsDome(fvis_distances_x_dome,fvis_distances_r_dome,fvispars_dome);
          }
          
          // cathode dimensions
          fcathode_ydimension = fActiveVolumes[0].SizeY();
          fcathode_zdimension = fActiveVolumes[0].SizeZ();
          
          // set cathode plane struct for solid angle function
          fcathode_plane.h = fcathode_ydimension;
          fcathode_plane.w = fcathode_zdimension;
          fplane_depth = std::abs(fCathode_centre.X());
          
      }
    }

    // LAr absorption length in cm
    std::map<double, double> abs_length_spectrum =
      lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
    std::vector<double> x_v, y_v;
    for (auto elem : abs_length_spectrum) {
      x_v.push_back(elem.first);
      y_v.push_back(elem.second);
    }
    fL_abs_vuv =
      std::round(interpolate(x_v, y_v, 9.7, false)); // 9.7 eV: peak of VUV emission spectrum

    // Load Gaisser-Hillas corrections for VUV semi-analytic hits
    std::cout << "Loading the GH corrections" << std::endl;
    fPVS->LoadGHForVUVCorrection(fGHvuvpars, fborder_corr, fradius);
    assert(!fGHvuvpars.empty());
    assert(!fborder_corr.empty());

    fdelta_angle = 10.; // angle bin size

    // Needed for Nhits-model border corrections (in cm)
    fYactive_corner = fActiveVolumes[0].HalfSizeY();
    fZactive_corner = fActiveVolumes[0].HalfSizeZ();

    fReference_to_corner =
      std::sqrt(fYactive_corner * fYactive_corner + fZactive_corner * fZactive_corner);

    std::cout << "For border corrections: " << fborder_corr[0] << "  " << fborder_corr[1]
              << std::endl;
    std::cout << "Photocathode-plane centre (z,y) = (" << fCathode_centre.Z() << ", "
              << fCathode_centre.Y() << ") "
              << "and corner (z, y) = (" << fZactive_corner << ", " << fYactive_corner << ")"
              << std::endl;
    std::cout << "Reference_to_corner: " << fReference_to_corner << std::endl;
  }

  //......................................................................
  // VUV semi-analytic hits calculation
  int
  PDFastSimPAR::VUVHits(const double Nphotons_created,
                        geo::Point_t const& ScintPoint,
                        OpticalDetector const& opDet)
  {
    // distance and angle between ScintPoint and OpDetPoint
    geo::Vector_t const relative = ScintPoint - opDet.OpDetPoint;
    const double distance = relative.R();
    const double cosine = std::abs(relative.X()) / distance;
    // const double theta = std::acos(cosine) * 180. / CLHEP::pi;
    const double theta = fast_acos(cosine) * 180. / CLHEP::pi;

    double solid_angle = 0.;
    // ARAPUCAS
    if (opDet.type == 0) {
      // get scintillation point coordinates relative to arapuca window centre
      geo::Vector_t const abs_relative{
        std::abs(relative.X()), std::abs(relative.Y()), std::abs(relative.Z())};
      solid_angle = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_relative);
    }
    // PMTs
    else if (opDet.type == 1) {
      const double zy_offset = std::sqrt(relative.Y() * relative.Y() + relative.Z() * relative.Z());
      const double x_distance = std::abs(relative.X());
      solid_angle = Disk_SolidAngle(zy_offset, x_distance, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" << std::endl;
    }

    // calculate number of photons hits by geometric acceptance: accounting for
    // solid angle and LAr absorbtion length
    double hits_geo =
      std::exp(-1. * distance / fL_abs_vuv) * (solid_angle / (4 * CLHEP::pi)) * Nphotons_created;

    // apply Gaisser-Hillas correction for Rayleigh scattering distance
    // and angular dependence offset angle bin
    const size_t j = (theta / fdelta_angle); // TODO: warning: conversion from double to int

    // Accounting for border effects
    double z_to_corner = std::abs(ScintPoint.Z() - fZactive_corner) - fZactive_corner;
    double y_to_corner = std::abs(ScintPoint.Y()) - fYactive_corner;
    double distance_to_corner =
      std::sqrt(y_to_corner * y_to_corner + z_to_corner * z_to_corner); // in the ph-cathode plane
    double pars_ini_[4] = {
      fGHvuvpars[0][j] + fborder_corr[0] * (distance_to_corner - fReference_to_corner),
      fGHvuvpars[1][j] + fborder_corr[1] * (distance_to_corner - fReference_to_corner),
      fGHvuvpars[2][j],
      fGHvuvpars[3][j]};
    double GH_correction = Gaisser_Hillas(distance, pars_ini_);
    int hits_vuv = fRandPoissPhot->fire(GH_correction * hits_geo / cosine);

    return hits_vuv;
  }

  //......................................................................
  // VIS hits semi-analytic model calculation
  int PDFastSimPAR::VISHits(const double Nphotons_created,   geo::Point_t const& ScintPoint, OpticalDetector const& opDet)
  {
      double plane_depth;
      if (ScintPoint.X() < 0.)
      {
          plane_depth = -fplane_depth;
      }
      else
      {
          plane_depth = fplane_depth;
      }
            
      // get scintpoint coords relative to centre of cathode plane
      geo::Vector_t const ScintPoint_relative = {std::abs(ScintPoint.X() - plane_depth),
                                                 std::abs(ScintPoint.Y() - fCathode_centre.Y()),
                                                 std::abs(ScintPoint.Z() - fCathode_centre.Z())};
      
      // calculate solid angle of cathode from the scintillation point
      double solid_angle_cathode = Rectangle_SolidAngle(fcathode_plane, ScintPoint_relative);
      
      // calculate distance and angle between ScintPoint and hotspot
      // vast majority of hits in hotspot region directly infront of scintpoint,
      // therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
      double distance_cathode = std::abs(plane_depth - ScintPoint.X());
      // calculate hits on cathode plane via geometric acceptance
      double cathode_hits_geo = std::exp(-1. * distance_cathode / fL_abs_vuv) *(solid_angle_cathode / (4. * CLHEP::pi)) * Nphotons_created;
      // determine Gaisser-Hillas correction including border effects
      // use flat correction
      double r = std::sqrt(std::pow(ScintPoint.Y() - fCathode_centre.Y(), 2) + std::pow(ScintPoint.Z() - fCathode_centre.Z(), 2));
      double pars_ini[4] = {0, 0, 0, 0};
      double s1 = 0; double s2 = 0; double s3 = 0;
      if(fIsFlatPDCorr) 
      {
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
      const double cathode_hits_rec = GH_correction * cathode_hits_geo;
      
      // detemine hits on each PD
      geo::Point_t const hotspot_v = {plane_depth, ScintPoint.Y(), ScintPoint.Z()};
      
      std::array<double, 3> hotspot;
      geo::vect::fillCoords(hotspot, hotspot_v);

    // 1). calculate total number of hits of VUV photons on reflective
    // foils via solid angle + Gaisser-Hillas corrections.
    // Done outside as it doesn't depend on OpDetPoint

    // 2). calculate number of these hits which reach the optical
    // detector from the hotspot via solid angle

//    geo::Vector_t const emission_relative = hotspot - opDet.OpDetPoint;

//    // calculate solid angle of optical channel
//    double solid_angle_detector = 0.;
//    // ARAPUCAS
//    if (opDet.type == 0) {
//      // get hotspot coordinates relative to opDet
//      geo::Vector_t const abs_emission_relative{std::abs(emission_relative.X()),
//                                                std::abs(emission_relative.Y()),
//                                                std::abs(emission_relative.Z())};
//      solid_angle_detector = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_emission_relative);
//    }
//    // PMTs
//    else if (opDet.type == 1) {
//      const double zy_offset = std::sqrt(emission_relative.Y() * emission_relative.Y() +
//                                         emission_relative.Z() * emission_relative.Z());
//      const double x_distance = std::abs(emission_relative.X());
//      solid_angle_detector = Disk_SolidAngle(zy_offset, x_distance, fradius);
//    }
//    else {
//      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" << std::endl;
//    }

//    // calculate number of hits via geometeric acceptance
//    double hits_geo = (solid_angle_detector / (2. * CLHEP::pi)) * cathode_hits_rec; // 2*pi due to presence of reflective foils

//    // calculate distances and angles for application of corrections
//    // distance to hotspot
//    const double distance_vuv = (ScintPoint - hotspot).R();
//    // distance from hotspot to optical detector
//    const double distance_vis = emission_relative.R();
//    //  angle between hotspot and optical detector
//    const double cosine_vis = std::abs(emission_relative.X()) / distance_vis;
//    // const double theta_vis = std::acos(cosine_vis) * 180. / CLHEP::pi;
//    const double theta_vis = fast_acos(cosine_vis) * 180. / CLHEP::pi;
//    
//    const size_t k = (theta_vis / fdelta_angle); // TODO: warning: conversion from double to int
//    // apply geometric correction
//    double pars_ini_vis[6] = {fvispars[0][k],
//                              fvispars[1][k],
//                              fvispars[2][k],
//                              fvispars[3][k],
//                              fvispars[4][k],
//                              fvispars[5][k]};
//    double geo_correction = Pol_5(distance_vuv, pars_ini_vis);
//    int hits_rec = fRandPoissPhot->fire(geo_correction * hits_geo / cosine_vis);
//    
//    // apply border correction
//    int hits_vis = 0;
//    if (fApplyVisBorderCorrection) {
//      // calculate distance for interpolation depending on model
//      double r = 0.;
//      if (fVisBorderCorrectionType == "Radial") {
//        r = std::sqrt(std::pow(ScintPoint.Y() - fCathode_centre.Y(), 2) +
//                      std::pow(ScintPoint.Z() - fCathode_centre.Z(), 2));
//      }
//      else if (fVisBorderCorrectionType == "Vertical") {
//        r = std::abs(ScintPoint.Y());
//      }
//      else {
//        std::cout << "Invalid border correction type - defaulting to using "
//                     "central value"
//                  << std::endl;
//      }
//      // interpolate in x for each r bin
//      const size_t nbins_r = fvis_border_correction[k].size();
//      std::vector<double> interp_vals(nbins_r, 0.0);
//      {
//        size_t idx = 0;
//        size_t size = fvis_border_distances_x.size();
//        if (std::abs(ScintPoint.X()) >= fvis_border_distances_x[size - 2])
//          idx = size - 2;
//        else {
//          while (std::abs(ScintPoint.X()) > fvis_border_distances_x[idx + 1])
//            idx++;
//        }
//        for (size_t i = 0; i < nbins_r; ++i) {
//          interp_vals[i] = interpolate(fvis_border_distances_x,
//                                       fvis_border_correction[k][i],
//                                       std::abs(ScintPoint.X()),
//                                       false,
//                                       idx);
//        }
//      }
//      std::cout << "funct 7"  << std::endl;
//      // interpolate in r
//      double border_correction = interpolate(fvis_border_distances_r, interp_vals, r, false);
//      // apply border correction
//      double hits_rec_borders = border_correction * hits_rec / cosine_vis;

//      // round final result
//      hits_vis = std::round(hits_rec_borders);
//    }
//    else {
//      hits_vis = hits_rec;
//    }

      // the interface has been converted into geo::Point_t, the implementation not yet
      std::array<double, 3U> SPoint;
      std::array<double, 3U> DPoint;
      geo::vect::fillCoords(SPoint, ScintPoint);
      geo::vect::fillCoords(DPoint, opDet.OpDetPoint);
      
      // calculate distances and angles for application of corrections
      // distance from hotspot to optical detector
      double distance_vis = dist(&hotspot[0], &DPoint[0], 3);
      //  angle between hotspot and optical detector
      double cosine_vis = std::abs(hotspot[0] - DPoint[0]) / distance_vis;
      // double theta_vis = std::acos(cosine_vis) * 180. / CLHEP::pi;
      double theta_vis = fast_acos(cosine_vis) * 180. / CLHEP::pi;
      
      // calculate solid angle of optical detector
      double solid_angle_detector = 0;
      // rectangular aperture
      if (opDet.type == 0) 
      {
          // get hotspot coordinates relative to opDet
          std::array<double, 3> emission_relative = {std::abs(hotspot[0] - DPoint[0]),
                                                     std::abs(hotspot[1] - DPoint[1]),
                                                     std::abs(hotspot[2] - DPoint[2])};

           geo::Vector_t const abs_emission_relative{std::abs(emission_relative[0]),
                                                     std::abs(emission_relative[1]),
                                                     std::abs(emission_relative[2])};
          // calculate solid angle
          solid_angle_detector = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_emission_relative);
      }
      // dome aperture
      else if (opDet.type == 1) 
      {
          solid_angle_detector = Omega_Dome_Model(distance_vis, theta_vis);
      }
      // disk aperture
      else if (opDet.type == 2) 
      {
           // offset in z-y plane
           double d = dist(&hotspot[1], &DPoint[1], 2);
           // drift distance (in x)
           double h = std::abs(hotspot[0] - DPoint[0]);
           // calculate solid angle
           solid_angle_detector = Disk_SolidAngle(d, h, fradius);
      }
      else 
      {
          std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk" << std::endl;
      }
      
      // calculate number of hits via geometeric acceptance
      double hits_geo = (solid_angle_detector / (2. * CLHEP::pi)) * cathode_hits_rec; // 2*pi due to presence of reflective foils
      
      // determine correction factor, depending on PD type
      const size_t k = (theta_vis / fdelta_angulo_vis);         // off-set angle bin
      TVector3 fcathode_centre;
      fcathode_centre[0] = fCathode_centre.X();
      fcathode_centre[1] = fCathode_centre.Y();
      fcathode_centre[2] = fCathode_centre.Z();
      double radial = dist(SPoint, fcathode_centre, 2, 1);      // radial distance from centre of detector (Y-Z)
      double d_c = std::abs(SPoint[0] - plane_depth);      // distance to cathode
      double border_correction = 0;
      // flat PDs
      if ((opDet.type == 0 || opDet.type == 2) && fIsFlatPDCorr)
      {
          // interpolate in d_c for each r bin
          const size_t nbins_r = fvispars_flat[k].size();
          std::vector<double> interp_vals(nbins_r, 0.0);

          size_t idx = 0;
          size_t size = fvis_distances_x_flat.size();
          if (d_c >= fvis_distances_x_flat[size - 2])
          {
              idx = size - 2;
          }
          
          else 
          {
              while (d_c > fvis_distances_x_flat[idx + 1])
              idx++;
          }
          for (size_t i = 0; i < nbins_r; ++i) 
          {
              interp_vals[i] = interpolate(fvis_distances_x_flat,
              fvispars_flat[k][i],
              d_c,
              false,
              idx);
          }
          // interpolate in r
          border_correction = interpolate(fvis_distances_r_flat, interp_vals, radial, false);
      }
      // dome PDs
      else if (opDet.type == 1 && fIsDomePDCorr) 
      {
          // interpolate in d_c for each r bin
          const size_t nbins_r = fvispars_dome[k].size();
          std::vector<double> interp_vals(nbins_r, 0.0);
          size_t idx = 0;
          size_t size = fvis_distances_x_dome.size();
          if (d_c >= fvis_distances_x_dome[size - 2])
          {
              idx = size - 2;
          }
          else
          {
              while (d_c > fvis_distances_x_dome[idx + 1])
              idx++;
          }
          for (size_t i = 0; i < nbins_r; ++i) 
          {
              interp_vals[i] = interpolate(fvis_distances_x_dome,
              fvispars_dome[k][i],
              d_c,
              false,
              idx);
          }
    
          // interpolate in r
          border_correction = interpolate(fvis_distances_r_dome, interp_vals, radial, false);
      }
      else 
      {
          std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk. Or corrections for chosen optical detector type missing." << std::endl;
      }
      
      // apply correction factor
      double hits_rec = border_correction * hits_geo / cosine_vis;
      double hits_vis = std::round(G4Poisson(hits_rec));
      
      return hits_vis;
  }

  bool
  PDFastSimPAR::isOpDetInSameTPC(geo::Point_t const& ScintPoint,
                                 geo::Point_t const& OpDetPoint) const
  {
    // check optical channel is in same TPC as scintillation light, if not return 0 hits
    // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero, otherwise return 0 hits
    if (((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) &&
        std::abs(OpDetPoint.X()) > 10.) { // TODO: unhardcode
      return false;
    }
    return true;
  }

  bool
  PDFastSimPAR::isScintInActiveVolume(geo::Point_t const& ScintPoint)
  {
    //semi-analytic approach only works in the active volume
    return fActiveVolumes[0].ContainsPosition(ScintPoint);
  }

  //......................................................................
  std::vector<double> PDFastSimPAR::PropTime(geo::Point_t ScintPoint, int OpChannel, int NPhotons, bool is_refl)
  {
      std::vector<double>  arrival_time_dist;
      arrival_time_dist.resize(NPhotons);
      
    if (fPVS->IncludeParPropTime() && fPVS->IncludePropTime()) {
      throw cet::exception("PDFastSimPAR")
        << "Cannot have both propagation time models simultaneously.";
    }
    else if (fPVS->IncludeParPropTime() &&
             !(ParPropTimeTF1 && (ParPropTimeTF1[OpChannel].GetNdim() == 1))) {
      // Warning: TF1::GetNdim()==1 will tell us if the TF1 is really defined or
      // it is the default one. This will fix a segfault when using timing and
      // interpolation.
      G4cout << "WARNING: Requested parameterized timing, but no function found. "
                "Not applying propagation time."
             << G4endl;
    }
    else if (fPVS->IncludeParPropTime()) {
      if (fPVS->StoreReflected())
        throw cet::exception("PDFastSimPAR")
          << "No parameterized propagation time for reflected light";
      for (size_t i = 0; i < arrival_time_dist.size(); ++i) {
        arrival_time_dist[i] = ParPropTimeTF1[OpChannel].GetRandom();
      }
    }
    else if (fPVS->IncludePropTime()) {
      // Get VUV photons arrival time distribution from the parametrization
      geo::Point_t const& opDetCenter = fOpDetCenter.at(OpChannel);
      if (!(fPVS->StoreReflected())) {
        const G4ThreeVector OpDetPoint(
          opDetCenter.X() * CLHEP::cm, opDetCenter.Y() * CLHEP::cm, opDetCenter.Z() * CLHEP::cm);

        const G4ThreeVector ScintPointG4TV (ScintPoint.X() * CLHEP::cm, ScintPoint.Y() * CLHEP::cm, ScintPoint.Z() * CLHEP::cm);
        double distance_in_cm = (ScintPointG4TV - OpDetPoint).mag() / CLHEP::cm; // this must be in CENTIMETERS!
        getVUVTimes(arrival_time_dist, distance_in_cm);              // in ns
      }
      else 
      {
            getVISTimes(arrival_time_dist, geo::vect::toTVector3(ScintPoint), geo::vect::toTVector3(opDetCenter)); // in ns
      }
    }

    return arrival_time_dist;
  }

  //......................................................................
  // VUV arrival times calculation function
  void
  PDFastSimPAR::getVUVTimes(std::vector<double>& arrivalTimes, double distance)
  {
    if (distance < 25) { // TODO: unhardcode
      // times are fixed shift i.e. direct path only
      double t_prop_correction = distance / fvuv_vgroup_mean;
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = t_prop_correction;
      }
    }
    else { // distance >= 25cm
      // determine nearest parameterisation in discretisation
      int index = std::round((distance - 25) / fstep_size); // TODO: unhardcode
      // check whether required parameterisation has been generated, generating if
      // not
      if (VUV_timing[index].GetNdim() == 0) { generateParam(index); }
      // randomly sample parameterisation for each photon
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = VUV_timing[index].GetRandom(VUV_min[index], VUV_max[index]);
      }
    }
  }

  //......................................................................
  // VIS arrival times calculation functions
  void
  PDFastSimPAR::getVISTimes(std::vector<double>& arrivalTimes,
                            TVector3 ScintPoint,
                            TVector3 OpDetPoint)
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

    // calculate point of reflection for shortest path accounting for difference
    // in refractive indices vectors for storing results
    TVector3 image(0, 0, 0);
    TVector3 bounce_point(0, 0, 0);

    // distance to wall
    TVector3 v_to_wall(plane_depth - ScintPoint[0], 0, 0);

    // hotspot is point on wall where TPB is activated most intensely by the
    // scintillation
    TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);

    // define "image" by reflecting over plane
    image = hotspot + v_to_wall * (fn_LAr_vis / fn_LAr_vuv);

    // find point of intersection with plane j of ray from the PMT to the image
    TVector3 tempvec = (OpDetPoint - image).Unit();
    double tempnorm = ((image - hotspot).Mag()) / std::abs(tempvec[0]);
    bounce_point = image + tempvec * tempnorm;

    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point - ScintPoint).Mag();
    double Visdist = (OpDetPoint - bounce_point).Mag();

    // calculate times taken by each part
    getVUVTimes(arrivalTimes, VUVdist);

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
    if (VUVdist < 25) { // TODO: unhardcode
      vuv_time = VUVdist / fvuv_vgroup_mean;
    }
    else {
      // find index of required parameterisation
      const size_t index = std::round((VUVdist - 25) / fstep_size); // TODO: unhardcode
      // find shortest time
      vuv_time = VUV_min[index];
    }
    // sum
    double fastest_time = vis_time + vuv_time;

    // calculate angle alpha between scintillation point and reflection point
    double cosine_alpha = std::abs(ScintPoint[0] - bounce_point[0]) / VUVdist;
    // double alpha = std::acos(cosine_alpha) * 180. / CLHEP::pi;
    double alpha = fast_acos(cosine_alpha) * 180. / CLHEP::pi;

    // determine smearing parameters using interpolation of generated points:
    // 1). tau = exponential smearing factor, varies with distance and angle
    // 2). cutoff = largest smeared time allowed, preventing excessively large
    //     times caused by exponential distance to cathode
    double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
    // angular bin
    size_t alpha_bin =
      alpha / 10; // TODO: unhardcode // TODO: warning: conversion from double to int
    if (alpha_bin >= ftau_pars.size()) {
      alpha_bin = ftau_pars.size() - 1; // default to the largest available bin if alpha
                                        // larger than parameterised region; i.e. last bin
                                        // effectively [last bin start value, 90] deg bin
    }
    // cut-off and tau
    double cutoff =
      interpolate(fdistances_refl, fcut_off_pars[alpha_bin], distance_cathode_plane, true);
    double tau = interpolate(fdistances_refl, ftau_pars[alpha_bin], distance_cathode_plane, true);
    if (tau < 0) { tau = 0; } // if tau extrapolate goes wrong

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
  PDFastSimPAR::generateParam(const size_t index)
  {
    // get distance
    double distance_in_cm = (index * fstep_size) + 25; // TODO: unhardcode

    // time range
    const double signal_t_range = 5000.; // TODO: unhardcode

    // parameterisation TF1
    TF1 fVUVTiming;

    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
    double t_direct_min = distance_in_cm / fvuv_vgroup_max;

    // Defining the model function(s) describing the photon transportation timing
    // vs distance Getting the landau parameters from the time parametrization
    std::array<double, 3> pars_landau;
    interpolate3(pars_landau,
                 fparameters[0],
                 fparameters[2],
                 fparameters[3],
                 fparameters[1],
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
      pars_expo[1] = interpolate(fparameters[4], fparameters[5], distance_in_cm, true);
      // For simplicity, not considering the small dependency with the offset
      // angle in pars_expo[0] Using the value for the [30,60deg] range.
      // fparameters[6] and fparameters[8] are the values for [0,30deg] range and
      // [60,90deg] range respectively
      pars_expo[0] = fparameters[7].at(0) + fparameters[7].at(1) * distance_in_cm;
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
        std::cout << "WARNING: Parametrization of VUV light discontinuous for "
                     "distance = "
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
    xq_max[0] = 0.99; // include 99%
    fVUVTiming.GetQuantiles(nq_max, yq_max, xq_max);
    double max = yq_max[0];
    // min
    double min = t_direct_min;

    // generate the sampling
    // the first call of GetRandom generates the timing sampling and stores it in
    // the TF1 object, this is the slow part all subsequent calls check if it has
    // been generated previously and are ~100+ times quicker add timing to the
    // vector of timings and range to vectors of ranges
    VUV_timing[index] = fVUVTiming;
    VUV_max[index] = max;
    VUV_min[index] = min;
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
  PDFastSimPAR::Rectangle_SolidAngle_Cal(const double a, const double b, const double d)
  {
    double aa = a / (2. * d);
    double bb = b / (2. * d);
    double aux = (1. + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
    // return 4 * std::acos(std::sqrt(aux));
    return 4. * fast_acos(std::sqrt(aux));
  }

  double
  PDFastSimPAR::Rectangle_SolidAngle(Dims const& o, geo::Vector_t const& v)
  {
    // v is the position of the track segment with respect to
    // the center position of the arapuca window

    // arapuca plane fixed in x direction
    if (isApproximatelyZero(v.Y()) && isApproximatelyZero(v.Z())) {
      return Rectangle_SolidAngle_Cal(o.h, o.w, v.X());
    }
    if (isDefinitelyGreaterThan(v.Y(), o.h * .5) && isDefinitelyGreaterThan(v.Z(), o.w * .5)) {
      double A = v.Y() - o.h * .5;
      double B = v.Z() - o.w * .5;
      double to_return = (Rectangle_SolidAngle_Cal(2. * (A + o.h), 2. * (B + o.w), v.X()) -
                          Rectangle_SolidAngle_Cal(2. * A, 2. * (B + o.w), v.X()) -
                          Rectangle_SolidAngle_Cal(2. * (A + o.h), 2. * B, v.X()) +
                          Rectangle_SolidAngle_Cal(2. * A, 2. * B, v.X())) *
                         .25;
      return to_return;
    }
    if ((v.Y() <= o.h * .5) && (v.Z() <= o.w * .5)) {
      double A = -v.Y() + o.h * .5;
      double B = -v.Z() + o.w * .5;
      double to_return = (Rectangle_SolidAngle_Cal(2. * (o.h - A), 2. * (o.w - B), v.X()) +
                          Rectangle_SolidAngle_Cal(2. * A, 2. * (o.w - B), v.X()) +
                          Rectangle_SolidAngle_Cal(2. * (o.h - A), 2. * B, v.X()) +
                          Rectangle_SolidAngle_Cal(2. * A, 2. * B, v.X())) *
                         .25;
      return to_return;
    }
    if (isDefinitelyGreaterThan(v.Y(), o.h * .5) && (v.Z() <= o.w * .5)) {
      double A = v.Y() - o.h * .5;
      double B = -v.Z() + o.w * .5;
      double to_return = (Rectangle_SolidAngle_Cal(2. * (A + o.h), 2. * (o.w - B), v.X()) -
                          Rectangle_SolidAngle_Cal(2. * A, 2. * (o.w - B), v.X()) +
                          Rectangle_SolidAngle_Cal(2. * (A + o.h), 2. * B, v.X()) -
                          Rectangle_SolidAngle_Cal(2. * A, 2. * B, v.X())) *
                         .25;
      return to_return;
    }
    if ((v.Y() <= o.h * .5) && isDefinitelyGreaterThan(v.Z(), o.w * .5)) {
      double A = -v.Y() + o.h * .5;
      double B = v.Z() - o.w * .5;
      double to_return = (Rectangle_SolidAngle_Cal(2. * (o.h - A), 2. * (B + o.w), v.X()) -
                          Rectangle_SolidAngle_Cal(2. * (o.h - A), 2. * B, v.X()) +
                          Rectangle_SolidAngle_Cal(2. * A, 2. * (B + o.w), v.X()) -
                          Rectangle_SolidAngle_Cal(2. * A, 2. * B, v.X())) *
                         .25;
      return to_return;
    }
    // error message if none of these cases, i.e. something has gone wrong!
    // std::cout << "Warning: invalid solid angle call." << std::endl;
    return 0.;
  }

  // ---------------------------------------------------------------------------
  std::vector<geo::BoxBoundedGeo>
  PDFastSimPAR::extractActiveVolumes(geo::GeometryCore const& geom)
  {
    std::vector<geo::BoxBoundedGeo> activeVolumes;
    activeVolumes.reserve(geom.Ncryostats());

    for (geo::CryostatGeo const& cryo : geom.IterateCryostats()) {

      // can't use it default-constructed since it would always include origin
      geo::BoxBoundedGeo box{cryo.TPC(0).ActiveBoundingBox()};

      for (geo::TPCGeo const& TPC : cryo.IterateTPCs())
        box.ExtendToInclude(TPC.ActiveBoundingBox());

      activeVolumes.push_back(std::move(box));

    } // for cryostats

    return activeVolumes;
  } // PDFastSimPAR::extractActiveVolumes()

  // ---------------------------------------------------------------------------
    double PDFastSimPAR::Omega_Dome_Model(const double distance, const double theta)
    {
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
        
        if(distance >= d_break)
        {
            double R_apparent_far = b - par1[j];
            return  (2*CLHEP::pi * (1 - std::sqrt(1 - std::pow(R_apparent_far/distance,2))));
        }
        else
        {
            double R_apparent_close = b - par0[j];
            return (2*CLHEP::pi * (1 - std::sqrt(1 - std::pow(R_apparent_close/distance,2))));
        }
    }

} // namespace phot

DEFINE_ART_MODULE(phot::PDFastSimPAR)
