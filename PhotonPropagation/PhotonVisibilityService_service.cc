////////////////////////////////////////////////////////////////////////
//
//  \file PhotonVisibilityService_service.cc
//
////////////////////////////////////////////////////////////////////////
//
//  Ben Jones, MIT 2012
//
//  This service reports the visibility of a particular point in
//  the detector to each OpDet.  This is used by the fast
//  optical simulation and by track-light association algorithms.
//
//  Visibility is defined as the fraction of isotropically produced
//  photons from a detector voxel which are expected to reach the 
//  OpDet in question.
//
//  This information is lookup up from a previousely generated
//  optical library file, whose path is specified to this service.
//
//  Note that it is important that the voxelization schemes match
//  between the library and the service instance for sensible results.
// 
//
// Framework includes

// LArSoft includes
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "PhotonPropagation/PhotonLibrary.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/OpDetGeo.h"

namespace phot{

  //--------------------------------------------------------------------
  PhotonVisibilityService::PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg) 
  {
    this->reconfigure(pset);
    mf::LogInfo("PhotonVisibilityService")<<"PhotonVisbilityService initializing"<<std::endl;
    
    fTheLibrary = new PhotonLibrary();
    
    art::ServiceHandle<geo::Geometry> geom;

    size_t NVoxels = GetVoxelDef().GetNVoxels();
    size_t NOpChannels = geom->NOpChannels();
    

    
    if((!fLibraryBuildJob)&&(!fDoNotLoadLibrary))
      {
	std::string LibraryFileWithPath;
	cet::search_path sp("FW_SEARCH_PATH");

	if( !sp.find_file(fLibraryFile, LibraryFileWithPath) )
	  throw cet::exception("PhotonVisibilityService") << "Unable to find photon library in "  << sp.to_string();

    if(!fParameterization) fTheLibrary->LoadLibraryFromFile(LibraryFileWithPath, NVoxels, NOpChannels);
      }
    else
      {
	mf::LogInfo("PhotonVisibilityService") << " Vis service running library build job.  Please ensure " 
					       << " job contains LightSource, LArG4, SimPhotonCounter"<<std::endl;
	fTheLibrary->CreateEmptyLibrary(NVoxels, NOpChannels);
      }
  }

  //--------------------------------------------------------------------
  PhotonVisibilityService::~PhotonVisibilityService() 
  {
  }
  
  void PhotonVisibilityService::StoreLibrary()
  {
    if(fLibraryBuildJob )
      {
	mf::LogInfo("PhotonVisibilityService") << " Vis service "
					       << " Storing Library entries to file..." <<std::endl;
	fTheLibrary->StoreLibraryToFile(fLibraryFile);
      }
  }
  

  //--------------------------------------------------------------------
  void PhotonVisibilityService::reconfigure(fhicl::ParameterSet const& p)
  {

    art::ServiceHandle<geo::Geometry> geom;
    
    // Library details
    fLibraryBuildJob      = p.get< bool        >("LibraryBuildJob"     );
    if(geom->DetId() == 3)
    {fParameterization     = p.get< bool        >("LBNE10ktParameterization"    );}
    else
    {fParameterization = false;}
    fLibraryFile          = p.get< std::string >("LibraryFile"         );
    fDoNotLoadLibrary     = p.get< bool        >("DoNotLoadLibrary"    );

    // Voxel parameters
    fUseCryoBoundary      = p.get< bool        >("UseCryoBoundary"     );
  	
    
    if(fUseCryoBoundary)
      {
	double CryoBounds[6];
	geom->CryostatBoundaries(CryoBounds);
	fXmin = CryoBounds[0];
	fXmax = CryoBounds[1];
	fYmin = CryoBounds[2];
	fYmax = CryoBounds[3];
	fZmin = CryoBounds[4];
	fZmax = CryoBounds[5];
      }
    else
      {
	fXmin      = p.get< double       >("XMin"     );
	fXmax      = p.get< double       >("XMax"     );
	fYmin      = p.get< double       >("YMin"     );
	fYmax      = p.get< double       >("YMax"     );
	fZmin      = p.get< double       >("ZMin"     );
	fZmax      = p.get< double       >("ZMax"     );
      }

    fNx        = p.get< int          >("NX"       );
    fNy        = p.get< int          >("NY"       );
    fNz        = p.get< int          >("NZ"       );
    
    fVoxelDef = sim::PhotonVoxelDef(fXmin, fXmax, fNx, fYmin, fYmax, fNy, fZmin, fZmax, fNz);

    return;
	
  }



  //------------------------------------------------------

  // Eventually we will calculate the light quenching factor here
  double PhotonVisibilityService::GetQuenchingFactor(double dQdx)
  {
    // for now, no quenching
    return 1.0;

  }


  //------------------------------------------------------

  // Get a vector of the relative visibilities of each OpDet
  //  in the event to a point xyz

  std::vector<float>* PhotonVisibilityService::GetAllVisibilities(double * xyz)
  {
    int VoxID = fVoxelDef.GetVoxelID(xyz);
    return GetLibraryEntries(VoxID);
  }


  //------------------------------------------------------

  // Get distance to optical detector OpDet
  double PhotonVisibilityService::DistanceToOpDet( double* xyz, unsigned int OpChannel )
  {
    art::ServiceHandle<geo::Geometry> geom;
 
    // Find the right OpDet
    unsigned int c=0, o=0;
    geom->OpChannelToCryoOpDet(OpChannel, o, c);

    // Get its coordinates
    return geom->Cryostat(c).OpDet(o).DistanceToPoint(xyz);
      
  }


  //------------------------------------------------------


  // Get the solid angle reduction factor for planar optical detector OpDet
  double PhotonVisibilityService::SolidAngleFactor( double* xyz, unsigned int OpChannel )
  {
    art::ServiceHandle<geo::Geometry> geom;
   
    // Find the right OpDet
    unsigned int c=0, o=0;
    geom->OpChannelToCryoOpDet(OpChannel, o, c);
    
    double CosTheta = geom->Cryostat(c).OpDet(o).CosThetaFromNormal(xyz);

    return CosTheta;
  }

  //------------------------------------------------------

  float PhotonVisibilityService::GetVisibility(double * xyz, unsigned int OpChannel)
  {
    int VoxID = fVoxelDef.GetVoxelID(xyz);  
    return GetLibraryEntry(VoxID, OpChannel);
  }


  //------------------------------------------------------

  void PhotonVisibilityService::StoreLightProd(int VoxID, double N)
  {
    fCurrentVoxel = VoxID;
    fCurrentValue = N;
    mf::LogInfo("PhotonVisibilityService") << " PVS notes production of " << N << " photons at Vox " << VoxID<<std::endl; 
  }


  //------------------------------------------------------

  
  void PhotonVisibilityService::RetrieveLightProd(int& VoxID, double& N)
  {
    N     = fCurrentValue;
    VoxID = fCurrentVoxel;
  }
  
  //------------------------------------------------------

  void PhotonVisibilityService::SetLibraryEntry(int VoxID, int OpChannel, float N)
  {
    fTheLibrary->SetCount(VoxID,OpChannel, N);
    mf::LogInfo("PhotonVisibilityService") << " PVS logging " << VoxID << " " << OpChannel<<std::endl;
  }

  //------------------------------------------------------

  

  std::vector<float>* PhotonVisibilityService::GetLibraryEntries(int VoxID)
  {
    return fTheLibrary->GetCounts(VoxID);
  }

  //------------------------------------------------------

  float PhotonVisibilityService::GetLibraryEntry(int VoxID, int Channel)
  {
    return fTheLibrary->GetCount(VoxID, Channel);
  }


} // namespace

namespace phot{
 
  DEFINE_ART_SERVICE(PhotonVisibilityService)

} // namespace phot
