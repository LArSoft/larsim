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
#include "Simulation/PhotonVoxels.h"

///General LArSoft Utilities
namespace phot{
  class PhotonLibrary;
  
  class PhotonVisibilityService {
  public:
    
    PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~PhotonVisibilityService();
    
    void reconfigure(fhicl::ParameterSet const& p);
    
    double GetQuenchingFactor(double dQdx);
    
    double DistanceToOpDet(                 double* xyz, unsigned int OpChannel );
    double SolidAngleFactor(                double* xyz, unsigned int OpChannel );
    float GetVisibility(                    double* xyz, unsigned int OpChannel );         

    std::vector<float>* GetAllVisibilities( double* xyz );
    
    void StoreLibrary();
    
    
    void StoreLightProd(    int  VoxID,  double  N );
    void RetrieveLightProd( int& VoxID,  double& N );
    
    void SetLibraryEntry(   int VoxID, int OpChannel, float N);
    float GetLibraryEntry( int VoxID, int OpChannel);
    std::vector<float>* GetLibraryEntries( int VoxID );

    
    bool IsBuildJob() { return fLibraryBuildJob; }
    bool UseParameterization() {return fParameterization;}

    sim::PhotonVoxelDef GetVoxelDef() {return fVoxelDef; }

  private:
    
    int    fCurrentVoxel;
    double fCurrentValue;
    
    float  fXmin, fXmax;
    float  fYmin, fYmax;
    float  fZmin, fZmax;
    int    fNx, fNy, fNz;

    bool fUseCryoBoundary;
    
    bool                 fLibraryBuildJob;
    bool                 fDoNotLoadLibrary;
    bool                 fParameterization;
    std::string          fLibraryFile;      
    PhotonLibrary *      fTheLibrary;
    sim::PhotonVoxelDef  fVoxelDef;
    
    
  }; // class PhotonVisibilityService
} //namespace utils
DECLARE_ART_SERVICE(phot::PhotonVisibilityService, LEGACY)
#endif // UTIL_DETECTOR_PROPERTIES_H
