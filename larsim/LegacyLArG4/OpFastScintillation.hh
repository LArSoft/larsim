// Class adapted for LArSoft by Ben Jones, MIT 10/10/12
//
// This class is a physics process based on the standard Geant4
// scintillation process.
//
// It has been stripped down and adapted to form the backbone of
// the LArG4 fast optical simulation.  Photons, instead of being
// produced and added to the geant4 particle stack, are logged
// and used to predict the visibility of this step to each PMT in
// the detector.
//
// The photonvisibilityservice looks up the visibility of the relevant
// xyz point, and if a photon is detected at a given PMT, one OnePhoton
// object is logged in the OpDetPhotonTable
//
// At the end of the event, the OpDetPhotonTable is read out
// by LArG4, and detected photons are stored in the event.
//
// This process can be used alongside the standard Cerenkov process,
// which does step geant4 opticalphotons.  Both the fast scintillation
// table and the geant4 sensitive detectors are read out by LArG4 to
// produce a combined SimPhoton collection.
//
//
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        OpFastScintillation.hh
// Description:	Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//                         of energy deposited by particle type
//                         Thanks to Zach Hartwig (Department of Nuclear
//                         Science and Engineeering - MIT)
//              2005-07-28 add G4ProcessType to constructor
//              2002-11-21 change to user G4Poisson for small MeanNumPotons
//              2002-11-07 allow for fast and slow scintillation
//              2002-11-05 make use of constant material properties
//              2002-05-16 changed to inherit from VRestDiscreteProcess
//              2002-05-09 changed IsApplicable method
//              1999-10-29 add method and class descriptors
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef OpFastScintillation_h
#define OpFastScintillation_h 1

/////////////
// Includes
/////////////

#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t
#include "larsim/PhotonPropagation/PhotonVisibilityTypes.h" // phot::MappedT0s_t

#include "Geant4/G4ForceCondition.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4PhysicsOrderedFreeVector.hh"
#include "Geant4/G4PhysicsTable.hh"
#include "Geant4/G4ProcessType.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Types.hh"
#include "Geant4/G4VRestDiscreteProcess.hh"

#include "TF1.h"
#include "TVector3.h"

#include <memory> // std::unique_ptr

class G4EmSaturation;
class G4Step;
class G4Track;
class G4VParticleChange;
namespace CLHEP {
  class RandGeneral;
}
namespace geo {
  class GeometryCore;
}
namespace phot {
  class PhotonVisibilityService;
}

// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VRestDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

namespace larg4 {

  class OpFastScintillation : public G4VRestDiscreteProcess {

  private:
    //////////////
    // Operators
    //////////////

    // OpFastScintillation& operator=(const OpFastScintillation &right);

  public: // Without description
    ////////////////////////////////
    // Constructors and Destructor
    ////////////////////////////////

    OpFastScintillation(const G4String& processName = "Scintillation",
                        G4ProcessType type = fElectromagnetic);

    ~OpFastScintillation();

    ////////////
    // Methods
    ////////////

  public: // With description
    // OpFastScintillation Process has both PostStepDoIt (for energy
    // deposition of particles in flight) and AtRestDoIt (for energy
    // given to the medium by particles at rest)

    virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
    // Returns true -> 'is applicable', for any particle type except
    // for an 'opticalphoton' and for short-lived particles

    G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);
    // Returns infinity; i. e. the process does not limit the step,
    // but sets the 'StronglyForced' condition for the DoIt to be
    // invoked at every step.

    G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*);
    // Returns infinity; i. e. the process does not limit the time,
    // but sets the 'StronglyForced' condition for the DoIt to be
    // invoked at every step.

    virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
    virtual G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep);

    // These are the methods implementing the scintillation process.

    void SetTrackSecondariesFirst(const G4bool state);
    // If set, the primary particle tracking is interrupted and any
    // produced scintillation photons are tracked next. When all
    // have been tracked, the tracking of the primary resumes.

    void SetFiniteRiseTime(const G4bool state);
    // If set, the OpFastScintillation process expects the user to have
    // set the constant material property FAST/SLOWSCINTILLATIONRISETIME.

    G4bool GetTrackSecondariesFirst() const;
    // Returns the boolean flag for tracking secondaries first.

    G4bool GetFiniteRiseTime() const;
    // Returns the boolean flag for a finite scintillation rise time.

    void SetScintillationYieldFactor(const G4double yieldfactor);
    // Called to set the scintillation photon yield factor, needed when
    // the yield is different for different types of particles. This
    // scales the yield obtained from the G4MaterialPropertiesTable.

    G4double GetScintillationYieldFactor() const;
    // Returns the photon yield factor.

    void SetScintillationExcitationRatio(const G4double excitationratio);
    // Called to set the scintillation exciation ratio, needed when
    // the scintillation level excitation is different for different
    // types of particles. This overwrites the YieldRatio obtained
    // from the G4MaterialPropertiesTable.

    G4double GetScintillationExcitationRatio() const;
    // Returns the scintillation level excitation ratio.

    G4PhysicsTable* GetFastIntegralTable() const;
    // Returns the address of the fast scintillation integral table.

    G4PhysicsTable* GetSlowIntegralTable() const;
    // Returns the address of the slow scintillation integral table.

    void AddSaturation(G4EmSaturation* sat) { emSaturation = sat; }
    // Adds Birks Saturation to the process.

    void RemoveSaturation() { emSaturation = NULL; }
    // Removes the Birks Saturation from the process.

    G4EmSaturation* GetSaturation() const { return emSaturation; }
    // Returns the Birks Saturation.

    void SetScintillationByParticleType(const G4bool);
    // Called by the user to set the scintillation yield as a function
    // of energy deposited by particle type

    G4bool GetScintillationByParticleType() const { return scintillationByParticleType; }
    // Return the boolean that determines the method of scintillation
    // production

    void DumpPhysicsTable() const;
    // Prints the fast and slow scintillation integral tables.

    /*std::vector<double> GetVUVTime(double, int);
      std::vector<double> GetVisibleTimeOnlyCathode(double, int);*/
    // old timings -- to be deleted

    void getVUVTimes(std::vector<double>& arrivalTimes,
                     const double distance_in_cm,
                     const size_t angle_bin);
    void generateParam(const size_t index, const size_t angle_bin);
    // Functions for vuv component Landau + Exponential timing parameterisation, updated method

    void getVISTimes(std::vector<double>& arrivalTimes,
                     const TVector3& ScintPoint,
                     const TVector3& OpDetPoint);
    // Visible component timing parameterisation

    void detectedDirectHits(std::map<size_t, int>& DetectedNum,
                            const double Num,
                            geo::Point_t const& ScintPoint) const;
    void detectedReflecHits(std::map<size_t, int>& ReflDetectedNum,
                            const double Num,
                            geo::Point_t const& ScintPoint) const;

  protected:
    void BuildThePhysicsTable();
    // It builds either the fast or slow scintillation integral table;
    // or both.

    bool RecordPhotonsProduced(const G4Step& aStep, double N);
    // Note the production of N photons in at point xyz.
    //  pass on to generate detector response, etc.

    ///////////////////////
    // Class Data Members
    ///////////////////////

    std::unique_ptr<G4PhysicsTable> theSlowIntegralTable;
    std::unique_ptr<G4PhysicsTable> theFastIntegralTable;

    G4bool fTrackSecondariesFirst;
    G4bool fFiniteRiseTime;

    G4double YieldFactor;

    G4double ExcitationRatio;

    G4bool scintillationByParticleType;

  private:
    struct OpticalDetector {
      double h; // height
      double w; // width
      geo::Point_t OpDetPoint;
      int type;
    };

    /// Returns whether the semi-analytic visibility parametrization is being used.
    bool usesSemiAnalyticModel() const;

    int VUVHits(const double Nphotons_created,
                geo::Point_t const& ScintPoint,
                OpticalDetector const& opDet) const;
    // Calculates semi-analytic model number of hits for vuv component

    int VISHits(geo::Point_t const& ScintPoint,
                OpticalDetector const& opDet,
                const double cathode_hits_rec,
                const std::array<double, 3> hotspot) const;
    // Calculates semi-analytic model number of hits for visible component

    G4double single_exp(const G4double t, const G4double tau2) const;
    G4double bi_exp(const G4double t, const G4double tau1, const G4double tau2) const;

    G4double scint_time(const G4Step& aStep,
                        G4double ScintillationTime,
                        G4double ScintillationRiseTime) const;
    void propagationTime(std::vector<double>& arrival_time_dist,
                         G4ThreeVector x0,
                         const size_t OpChannel,
                         bool Reflected = false); //const;

    // emission time distribution when there is a finite rise time
    G4double sample_time(const G4double tau1, const G4double tau2) const;

    // Facility for TPB emission energies
    double reemission_energy() const;
    std::map<double, double> tpbemission;
    std::unique_ptr<CLHEP::RandGeneral> fTPBEm;

    void average_position(G4Step const& aStep, double* xzyPos) const;

    G4EmSaturation* emSaturation;
    // functions and parameters for the propagation time parametrization
    phot::MappedFunctions_t ParPropTimeTF1;
    phot::MappedT0s_t ReflT0s;

    /*TF1 const* functions_vuv[8];
     TF1 const* functions_vis[5];
     double fd_break;
     double fd_max;
     double ftf1_sampling_factor;
     double ft0_max, ft0_break_point;*/

    //For new VUV time parametrization
    double fstep_size, fmin_d, fmax_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance,
      fangle_bin_timing_vuv;
    std::vector<std::vector<double>> fparameters[7];
    // vector containing generated VUV timing parameterisations
    std::vector<std::vector<TF1>> VUV_timing;
    // vector containing min and max range VUV timing parameterisations are sampled to
    std::vector<std::vector<double>> VUV_max;
    std::vector<std::vector<double>> VUV_min;

    // For new VIS time parameterisation
    double fvis_vmean, fangle_bin_timing_vis;
    std::vector<double> fdistances_refl;
    std::vector<double> fradial_distances_refl;
    std::vector<std::vector<std::vector<double>>> fcut_off_pars;
    std::vector<std::vector<std::vector<double>>> ftau_pars;

    struct Dims {
      double h, w; // height, width
    };

    // solid angle of rectangular aperture calculation functions
    double Rectangle_SolidAngle(const double a, const double b, const double d) const;
    double Rectangle_SolidAngle(Dims const& o, const std::array<double, 3> v) const;
    // solid angle of circular aperture calculation functions
    double Disk_SolidAngle(const double d, const double h, const double b) const;
    // solid angle of a dome aperture calculation functions
    double Omega_Dome_Model(const double distance, const double theta) const;

    // For VUV semi-analytic hits
    // Gaisser-Hillas correction parameters for VUV Nhits estimation
    G4double Gaisser_Hillas(const double x, const double* par) const;
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

    // geometry properties
    double fplane_depth, fcathode_zdimension, fcathode_ydimension;
    TVector3 fcathode_centre;
    std::vector<geo::BoxBoundedGeo> const fActiveVolumes;

    // Optical detector properties for semi-analytic hits
    // int foptical_detector_type;  // unused
    double fradius;
    Dims fcathode_plane;
    int fL_abs_vuv;
    std::vector<geo::Point_t> fOpDetCenter;
    std::vector<int> fOpDetType;
    std::vector<double> fOpDetLength;
    std::vector<double> fOpDetHeight;
    //double fGlobalTimeOffset;

    void ProcessStep(const G4Step& step);

    bool const bPropagate; ///< Whether propagation of photons is enabled.

    /// Photon visibility service instance.
    phot::PhotonVisibilityService const* const fPVS;

    /// Whether the semi-analytic model is being used for photon visibility.
    bool const fUseNhitsModel = false;
    /// Whether photon propagation is performed only from active volumes
    bool const fOnlyActiveVolume = false;
    /// Allows running even if light on cryostats `C:1` and higher is not supported.
    /// Currently hard coded "no"
    bool const fOnlyOneCryostat = false;
    /// Whether the cathodes are fully opaque; currently hard coded "no".
    bool const fOpaqueCathode = false;

    bool isOpDetInSameTPC(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) const;
    bool isScintInActiveVolume(geo::Point_t const& ScintPoint);
    double interpolate(const std::vector<double>& xData,
                       const std::vector<double>& yData,
                       double x,
                       bool extrapolate,
                       size_t i = 0) const;
    void interpolate3(std::array<double, 3>& inter,
                      const std::vector<double>& xData,
                      const std::vector<double>& yData1,
                      const std::vector<double>& yData2,
                      const std::vector<double>& yData3,
                      double x,
                      bool extrapolate) const;

    static std::vector<geo::BoxBoundedGeo> extractActiveVolumes(geo::GeometryCore const& geom);

  }; // class OpFastScintillation

  double finter_d(double*, double*);
  double LandauPlusExpoFinal(double*, double*);
  double model_close(double*, double*);
  double model_far(double*, double*);

  static const size_t acos_bins = 2000000;
  constexpr double acos_table(const double x);
  double fast_acos(const double x);

  ////////////////////
  // Inline methods
  ////////////////////
  inline G4bool OpFastScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
  {
    if (aParticleType.GetParticleName() == "opticalphoton") return false;
    if (aParticleType.IsShortLived()) return false;

    return true;
  }

  inline void OpFastScintillation::SetTrackSecondariesFirst(const G4bool state)
  {
    fTrackSecondariesFirst = state;
  }

  inline void OpFastScintillation::SetFiniteRiseTime(const G4bool state)
  {
    fFiniteRiseTime = state;
  }

  inline G4bool OpFastScintillation::GetTrackSecondariesFirst() const
  {
    return fTrackSecondariesFirst;
  }

  inline G4bool OpFastScintillation::GetFiniteRiseTime() const { return fFiniteRiseTime; }

  inline void OpFastScintillation::SetScintillationYieldFactor(const G4double yieldfactor)
  {
    YieldFactor = yieldfactor;
  }

  inline G4double OpFastScintillation::GetScintillationYieldFactor() const { return YieldFactor; }

  inline void OpFastScintillation::SetScintillationExcitationRatio(const G4double excitationratio)
  {
    ExcitationRatio = excitationratio;
  }

  inline G4double OpFastScintillation::GetScintillationExcitationRatio() const
  {
    return ExcitationRatio;
  }

  inline G4PhysicsTable* OpFastScintillation::GetSlowIntegralTable() const
  {
    return theSlowIntegralTable.get();
  }

  inline G4PhysicsTable* OpFastScintillation::GetFastIntegralTable() const
  {
    return theFastIntegralTable.get();
  }

  inline void OpFastScintillation::DumpPhysicsTable() const
  {
    if (theFastIntegralTable) {
      G4int PhysicsTableSize = theFastIntegralTable->entries();
      G4PhysicsOrderedFreeVector* v;
      for (G4int i = 0; i < PhysicsTableSize; i++) {
        v = (G4PhysicsOrderedFreeVector*)(*theFastIntegralTable)[i];
        v->DumpValues();
      }
    }
    if (theSlowIntegralTable) {
      G4int PhysicsTableSize = theSlowIntegralTable->entries();
      G4PhysicsOrderedFreeVector* v;
      for (G4int i = 0; i < PhysicsTableSize; i++) {
        v = (G4PhysicsOrderedFreeVector*)(*theSlowIntegralTable)[i];
        v->DumpValues();
      }
    }
  }

  template <typename TReal>
  inline constexpr double dist(const TReal* x, const TReal* y, const unsigned int dimension)
  {
    double d = 0.;
    for (unsigned int p = 0; p < dimension; ++p) {
      d += (*(x + p) - *(y + p)) * (*(x + p) - *(y + p));
    }
    return std::sqrt(d);
  }

  template <typename TVector3>
  inline constexpr double dist(const std::array<double, 3> x,
                               const TVector3 y,
                               const unsigned int dimension,
                               const unsigned int start)
  {
    double d = 0.;
    for (unsigned int p = start; p < dimension; ++p) {
      d += (x[p] - y[p]) * (x[p] - y[p]);
    }
    return std::sqrt(d);
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
  // for example, default tolerance may not work if you are comparing double with float
  template <typename TReal>
  inline constexpr static bool isApproximatelyZero(
    TReal a,
    TReal tolerance = std::numeric_limits<TReal>::epsilon())
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

  // template<typename Function, typename... Args>
  // auto OpFastScintillation::invoke_memoized(Function function, Args... args)
  // {
  //   using key_type   = std::tuple<Args...>;
  //   using value_type = std::invoke_result_t<Function, Args...>;
  //   static_assert(! std::is_same_v<Function, std::function<value_type(Args...)>>,
  //                 "cannot memoize on std::function (use a lambda instead)");
  //   static_assert(! std::is_same_v<Function, value_type(*)(Args...)>,
  //                 "cannot memoize on function pointer (use a lambda instead)");
  //   static std::mutex mutex;
  //   static std::map<key_type, value_type> cache;
  //   auto key  = std::tuple(args...);
  //   auto lock = std::lock_guard<std::mutex>(mutex);
  //   if (cache.count(key)){
  //     return cache[key];
  //   }
  //   return cache[key] = std::apply(function, key);
  // }

} // namespace larg4

#endif /* OpFastScintillation_h */
