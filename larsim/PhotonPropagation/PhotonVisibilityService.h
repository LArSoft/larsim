////////////////////////////////////////////////////////////////////////
// \file PhotonVisibilityService.h
//
// \brief Service to report opdet visibility to different points in
//         the system
//
// \author bjpjones@mit.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef PHOTONVISIBILITYSERVICE_H
#define PHOTONVISIBILITYSERVICE_H


#include "larcorealg/Geometry/geo_vectors_utils.h"          // geo::vect namespace
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t
#include "larsim/PhotonPropagation/IPhotonLibrary.h"
//#include "larsim/PhotonPropagation/LibraryMappingTools/IPhotonMappingTransformations.h"
//#include "larsim/PhotonPropagation/PhotonVisibilityTypes.h"
#include "larsim/Simulation/PhotonVoxels.h"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace fhicl {
  class ParameterSet;
}

class TF1;

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <string>
#include <vector>

///General LArSoft Utilities
namespace phot{
  
  class PhotonVisibilityService {
  public:
    
    ~PhotonVisibilityService();
    PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    
    void reconfigure(fhicl::ParameterSet const& p);
    
    double GetQuenchingFactor(double dQdx) const;
    
    static double DistanceToOpDet(          double const* xyz, unsigned int OpDet );
    static double SolidAngleFactor(         double const* xyz, unsigned int OpDet );
    float GetVisibility(                    double const* xyz, unsigned int OpChannel, bool wantReflected=false ) const;

    float const* GetAllVisibilities( double const* xyz, bool wantReflected=false ) const;
    
    void LoadLibrary() const;
    void StoreLibrary();
    
    
    void StoreLightProd(    int  VoxID,  double  N );
    void RetrieveLightProd( int& VoxID,  double& N ) const;
    
    void SetLibraryEntry(  int VoxID, int OpChannel, float N, bool wantReflected=false );
    float GetLibraryEntry( int VoxID, int OpChannel, bool wantReflected=false ) const;
    float const* GetLibraryEntries( int VoxID, bool wantReflected=false ) const;

    float const* GetReflT0s( double const* xyz ) const;
    void SetLibraryReflT0Entry( int VoxID, int OpChannel, float value );
    float const* GetLibraryReflT0Entries( int VoxID ) const;
    float GetLibraryReflT0Entry( int VoxID, int Channel ) const;
 
    const std::vector<float>* GetTimingPar( double const* xyz ) const;
    void SetLibraryTimingParEntry( int VoxID, int OpChannel, float value, size_t parnum );
    const std::vector<float>* GetLibraryTimingParEntries( int VoxID ) const;
    float GetLibraryTimingParEntry( int VoxID, int Channel, size_t npar ) const;

    TF1* GetTimingTF1( double const* xyz ) const;
    void SetLibraryTimingTF1Entry( int VoxID, int OpChannel, TF1 func );
    TF1* GetLibraryTimingTF1Entries( int VoxID ) const;
 
    void SetDirectLightPropFunctions(TF1 const* functions[8], double& d_break, double& d_max, double& tf1_sampling_factor) const;
    void SetReflectedCOLightPropFunctions(TF1 const* functions[5], double& t0_max, double& t0_break_point) const;
    
    void LoadVUVSemiAnalyticProperties(bool& isFlatPDCorr,
                                       double& delta_angulo_vuv,
                                       double& radius) const;
    void LoadGHFlat(std::vector<std::vector<double>>& GHvuvpars_flat,
                    std::vector<double>& border_corr_angulo_flat,
                    std::vector<std::vector<double>>& border_corr_flat) const;
    void LoadVisSemiAnalyticProperties(double& delta_angulo_vis, double& radius) const;

    bool IsBuildJob() const { return fLibraryBuildJob; }
    bool UseParameterization() const {return fParameterization;}
    bool StoreReflected() const { return fStoreReflected; }
    bool StoreReflT0() const { return fStoreReflT0; }
    bool IncludeParPropTime() const { return fParPropTime; }
    size_t ParPropTimeNpar() const { return fParPropTime_npar; }
    std::string ParPropTimeFormula() const { return fParPropTime_formula; }

    bool IncludePropTime() const { return fIncludePropTime; }
    bool UseNhitsModel() const { return fUseNhitsModel; }

    const sim::PhotonVoxelDef& GetVoxelDef() const {return fVoxelDef; }
    size_t NOpChannels() const;
    
  private:

    const TVector3 LibLocation(const double * xyz) const;

    int    fCurrentVoxel;
    double fCurrentValue;
    // for c2: fCurrentReflValue is unused
    //double fCurrentReflValue;

    float  fXmin, fXmax;
    float  fYmin, fYmax;
    float  fZmin, fZmax;
    int    fNx, fNy, fNz;

    bool fUseCryoBoundary;
    
    bool                 fLibraryBuildJob;
    bool                 fDoNotLoadLibrary;
    bool                 fParameterization;
    bool                 fHybrid;
    bool                 fStoreReflected;
    bool                 fStoreReflT0;
    bool                 fIncludePropTime;
    bool                 fUseNhitsModel;

    bool                 fParPropTime;
    size_t               fParPropTime_npar;
    std::string		 fParPropTime_formula;
    int                  fParPropTime_MaxRange;
    bool                 fInterpolate;
    bool                 fReflectOverZeroX;

    TF1 *fparslogNorm = nullptr;
    TF1 *fparslogNorm_far = nullptr;
    TF1 *fparsMPV = nullptr;
    TF1 *fparsMPV_far = nullptr;
    TF1 *fparsWidth = nullptr;
    TF1 *fparsCte = nullptr;
    TF1 *fparsCte_far = nullptr;
    TF1 *fparsSlope = nullptr;
    double fD_break, fD_max, fTF1_sampling_factor;
    TF1 *fparslogNorm_refl = nullptr;
    TF1 *fparsMPV_refl = nullptr;
    TF1 *fparsWidth_refl = nullptr;
    TF1 *fparsCte_refl = nullptr;
    TF1 *fparsSlope_refl = nullptr;
    double fT0_max, fT0_break_point;
   
    std::string          fLibraryFile;      
    mutable IPhotonLibrary* fTheLibrary;
    sim::PhotonVoxelDef  fVoxelDef;

    //for the semi-analytic vuv/direct light signal (number of hits) correction
    bool fIsFlatPDCorr;

    double fdelta_angulo_vuv;
    // flat PDs
    std::vector<std::vector<double>> fGHvuvpars_flat;
    std::vector<double> fborder_corr_angulo_flat;
    std::vector<std::vector<double>> fborder_corr_flat;

    // optical detector information, rest using geometry service
    double fradius;

    // --- BEGIN Implementation functions --------------------------------------
    /// @name Implementation functions
    /// @{

    static double DistanceToOpDetImpl(geo::Point_t const& p, unsigned int OpDet);

    static double SolidAngleFactorImpl(geo::Point_t const& p, unsigned int OpDet);
    
  }; // class PhotonVisibilityService
} //namespace phot
DECLARE_ART_SERVICE(phot::PhotonVisibilityService, LEGACY)
#endif // UTIL_DETECTOR_PROPERTIES_H
