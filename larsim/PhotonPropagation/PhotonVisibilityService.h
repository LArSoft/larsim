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


#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "larsim/PhotonPropagation/IPhotonLibrary.h"
#include "larsim/Simulation/PhotonVoxels.h"

class TF1;

///General LArSoft Utilities
namespace phot{
  
  class PhotonVisibilityService {
  public:
    
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
    
    bool IsBuildJob() const { return fLibraryBuildJob; }
    bool UseParameterization() const {return fParameterization;}
    bool StoreReflected() const { return fStoreReflected; }
    bool StoreReflT0() const { return fStoreReflT0; }
    bool IncludeParPropTime() const { return fParPropTime; }
    size_t ParPropTimeNpar() const { return fParPropTime_npar; }
    std::string ParPropTimeFormula() const { return fParPropTime_formula; }

    bool IncludePropTime() const { return fIncludePropTime; }

    const sim::PhotonVoxelDef& GetVoxelDef() const {return fVoxelDef; }
    size_t NOpChannels() const;
    
  private:
    
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

    bool                 fParPropTime;
    size_t               fParPropTime_npar;
    std::string		 fParPropTime_formula;

    bool                 fInterpolate;

    TF1 *fparslogNorm;
    TF1 *fparslogNorm_far;
    TF1 *fparsMPV;
    TF1 *fparsMPV_far;
    TF1 *fparsWidth;
    TF1 *fparsCte;
    TF1 *fparsCte_far;
    TF1 *fparsSlope;
    double fD_break, fD_max, fTF1_sampling_factor;
    TF1 *fparslogNorm_refl;
    TF1 *fparsMPV_refl;
    TF1 *fparsWidth_refl;
    TF1 *fparsCte_refl;
    TF1 *fparsSlope_refl;
    double fT0_max, fT0_break_point;
   
    std::string          fLibraryFile;      
    mutable IPhotonLibrary* fTheLibrary;
    sim::PhotonVoxelDef  fVoxelDef;
    
    
  }; // class PhotonVisibilityService
} //namespace phot
DECLARE_ART_SERVICE(phot::PhotonVisibilityService, LEGACY)
#endif // UTIL_DETECTOR_PROPERTIES_H
