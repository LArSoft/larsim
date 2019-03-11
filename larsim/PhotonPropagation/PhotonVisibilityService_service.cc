// ////////////////////////////////////////////////////////////////////////
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
//  This information is lookup up from a previously generated
//  optical library file, whose path is specified to this service.
//
//  Note that it is important that the voxelization schemes match
//  between the library and the service instance for sensible results.
// 
//
// Framework includes

// LArSoft includes
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PhotonLibrary.h"
#include "larsim/Simulation/PhotonVoxels.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "larsim/PhotonPropagation/PhotonLibrary.h"
#include "larsim/PhotonPropagation/PhotonLibraryHybrid.h"

#include "TF1.h"

#include <math.h>

namespace phot{

  PhotonVisibilityService::~PhotonVisibilityService()
  {
    delete fparslogNorm;
    delete fparslogNorm_far;
    delete fparsMPV;
    delete fparsMPV_far;
    delete fparsWidth;
    delete fparsCte;
    delete fparsCte_far;
    delete fparsSlope;
    delete fparslogNorm_refl;
    delete fparsMPV_refl;
    delete fparsWidth_refl;
    delete fparsCte_refl;
    delete fparsSlope_refl;
    delete fTheLibrary;
  }

  //--------------------------------------------------------------------
  PhotonVisibilityService::PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry &/*reg*/) :
    
    fCurrentVoxel(0),
    fCurrentValue(0.),
    fXmin(0.),
    fXmax(0.),
    fYmin(0.),
    fYmax(0.),
    fZmin(0.),
    fZmax(0.),
    fNx(0),
    fNy(0),
    fNz(0),
    fUseCryoBoundary(false),
    fLibraryBuildJob(false),
    fDoNotLoadLibrary(false),
    fParameterization(false),
    fHybrid(false),
    fStoreReflected(false),
    fStoreReflT0(false),
    fIncludePropTime(false),
    fUseNhitsModel(false),
    fParPropTime(false),
    fParPropTime_npar(0),
    fParPropTime_formula(),
    fParPropTime_MaxRange(),
    fInterpolate(false),
    fReflectOverZeroX(false),
    fparslogNorm(nullptr),
    fparslogNorm_far(nullptr),
    fparsMPV(nullptr),
    fparsMPV_far(nullptr),
    fparsWidth(nullptr),
    fparsCte(nullptr),
    fparsCte_far(nullptr),
    fparsSlope(nullptr),
    fD_break(0.0),
    fD_max(0.0),
    fTF1_sampling_factor(0.0),
    fparslogNorm_refl(nullptr),
    fparsMPV_refl(nullptr),
    fparsWidth_refl(nullptr),
    fparsCte_refl(nullptr),
    fparsSlope_refl(nullptr),
    fT0_max(0.0),
    fT0_break_point(0.0),
    fLibraryFile(),
    fTheLibrary(nullptr),
    fVoxelDef()
  {
    this->reconfigure(pset);
    mf::LogInfo("PhotonVisibilityService")<<"PhotonVisbilityService initializing"<<std::endl;
  }

  //--------------------------------------------------------------------
  void PhotonVisibilityService::LoadLibrary() const
  {
    // Don't do anything if the library has already been loaded.

    if(fTheLibrary == 0) {

      if((!fLibraryBuildJob)&&(!fDoNotLoadLibrary)) {
	std::string LibraryFileWithPath;
	cet::search_path sp("FW_SEARCH_PATH");

	if( !sp.find_file(fLibraryFile, LibraryFileWithPath) )
	  throw cet::exception("PhotonVisibilityService") << "Unable to find photon library in "  << sp.to_string() << "\n";

	if(!fParameterization) {
          art::ServiceHandle<geo::Geometry> geom;

	  mf::LogInfo("PhotonVisibilityService") << "PhotonVisibilityService Loading photon library from file "
						 << LibraryFileWithPath
                                                 << " for "
                                                 << GetVoxelDef().GetNVoxels()
                                                 << " voxels and "
                                                 << geom->NOpDets()
                                                 << " optical detectors."
						 << std::endl;

          if(fHybrid){
            fTheLibrary = new PhotonLibraryHybrid(LibraryFileWithPath,
                                                  GetVoxelDef());
          }
          else{
            PhotonLibrary* lib = new PhotonLibrary;
            fTheLibrary = lib;

            size_t NVoxels = GetVoxelDef().GetNVoxels();
            lib->LoadLibraryFromFile(LibraryFileWithPath, NVoxels, fStoreReflected, fStoreReflT0, fParPropTime_npar, fParPropTime_MaxRange);
          }
	}
      }
      else {
        art::ServiceHandle<geo::Geometry> geom;

        size_t NOpDets = geom->NOpDets();
        size_t NVoxels = GetVoxelDef().GetNVoxels();
	mf::LogInfo("PhotonVisibilityService") << " Vis service running library build job.  Please ensure " 
					       << " job contains LightSource, LArG4, SimPhotonCounter"<<std::endl;
        PhotonLibrary* lib = new PhotonLibrary;
        fTheLibrary = lib;

        lib->CreateEmptyLibrary(NVoxels, NOpDets, fStoreReflected, fStoreReflT0, fParPropTime_npar);
      }

    }
  }

  //--------------------------------------------------------------------
  void PhotonVisibilityService::StoreLibrary()
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    if(fLibraryBuildJob )
      {
          
          if(fHybrid){
              std::cout<< "This is would be building a Hybrid Library. Not defined. "<<std::endl;
          }
          mf::LogInfo("PhotonVisibilityService") << " Vis service "
					       << " Storing Library entries to file..." <<std::endl;
          PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
          lib->StoreLibraryToFile(fLibraryFile, fStoreReflected, fStoreReflT0, fParPropTime_npar);
      }
  }
  

  //--------------------------------------------------------------------
  void PhotonVisibilityService::reconfigure(fhicl::ParameterSet const& p)
  {

    art::ServiceHandle<geo::Geometry> geom;
    
    // Library details
    fLibraryBuildJob      = p.get< bool        >("LibraryBuildJob", false);
    fParameterization     = p.get< bool        >("DUNE10ktParameterization", false);
    fHybrid               = p.get< bool        >("HybridLibrary", false);
    fLibraryFile          = p.get< std::string >("LibraryFile", "");
    fDoNotLoadLibrary     = p.get< bool        >("DoNotLoadLibrary");
    fStoreReflected       = p.get< bool        >("StoreReflected", false);
    fStoreReflT0          = p.get< bool        >("StoreReflT0",  false);
    // Parametrizations (time and Nhits)
    fIncludePropTime      = p.get< bool        >("IncludePropTime", false);
    fUseNhitsModel        = p.get< bool        >("UseNhitsModel", false);
    // Voxel parameters
    fUseCryoBoundary      = p.get< bool        >("UseCryoBoundary", false);
    fInterpolate          = p.get< bool        >("Interpolate", false);
    fReflectOverZeroX     = p.get< bool        >("ReflectOverZeroX", false);

    fParPropTime          = p.get< bool        >("ParametrisedTimePropagation", false);
    fParPropTime_npar     = p.get< size_t      >("ParametrisedTimePropagationNParameters", 0);
    fParPropTime_formula  = p.get< std::string >("ParametrisedTimePropagationFittedFormula","");
    fParPropTime_MaxRange = p.get< int         >("ParametrisedTimePropagationMaxRange", 200);
    
    if (!fParPropTime)
    {
      fParPropTime_npar=0;
    }

    if(!fUseNhitsModel) {

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

    }

    if(fIncludePropTime)
      {

	// load VUV arrival time distribution parametrization (no detector dependent at first order)
	std::cout<<"Loading the VUV time parametrization"<<std::endl;
	fDistances_all     = p.get<std::vector<double> >("Distances_landau");
	fNorm_over_entries = p.get<std::vector<double> >("Norm_over_entries");
	fMpv = p.get<std::vector<double> >("Mpv");
	fWidth = p.get<std::vector<double> >("Width");
	fDistances = p.get<std::vector<double> >("Distances_exp");
	fSlope = p.get<std::vector<double> >("Slope");
	fExpo_over_Landau_norm[0] = p.get<std::vector<double> >("Expo_over_Landau_norm_0");
	fExpo_over_Landau_norm[1] = p.get<std::vector<double> >("Expo_over_Landau_norm_30");
	fExpo_over_Landau_norm[2] = p.get<std::vector<double> >("Expo_over_Landau_norm_60");
	fstep_size = p.get<double>("step_size");
	fmax_d= p.get<double>("max_d");
	fvuv_vgroup_mean= p.get<double>("vuv_vgroup_mean");
	fvuv_vgroup_max= p.get<double>("vuv_vgroup_max");
	finflexion_point_distance= p.get<double>("inflexion_point_distance");

	if (fStoreReflected) 
	{

	  // load VIS arrival time distribution paramterisation
	  std::cout << "Loading the VIS time paramterisation" << std::endl;
	  fDistances_refl = p.get<std::vector<double>> ("Distances_refl");
	  fCut_off = p.get<std::vector<std::vector<double>>> ("Cut_off");
	  fTau = p.get<std::vector<std::vector<double>>> ("Tau");
	  fvis_vmean = p.get<double> ("vis_vmean");
	  fn_LAr_VUV = p.get<double> ("n_LAr_VUV");
	  fn_LAr_vis = p.get<double> ("n_LAr_vis");
	  fPlane_Depth = p.get<double>("Plane_Depth");

	}
	//ALL BELOW IS OLD PARAMETRIZATION. TO BE REMOVED SOON (reflected component time needs to be updated too)
	// Construct parameterized model parameter functions.         
	/*std::cout<< "Getting direct light parameters from .fcl file"<<std::endl;
	std::vector<std::string> direct_functions = p.get<std::vector<std::string> >("Direct_functions");
	//range of distances where the parametrization is valid                                            
	fD_break = p.get<double>("D_break");
	fD_max = p.get<double>("D_max");

	fTF1_sampling_factor = p.get<double>("TF1_sampling_factor");	

	std::vector<double> direct_landauNormpars = p.get<std::vector<double> >("Direct_landauNormpars");
	fparslogNorm = new TF1("fparslogNorm", direct_functions[0].c_str(), 0., fD_break);
	for(unsigned int i=0; i<direct_landauNormpars.size(); ++i)
	  fparslogNorm->SetParameter(i, direct_landauNormpars[i]);
	
	std::vector<double> direct_landauMPVpars = p.get<std::vector<double> >("Direct_landauMPVpars");
	fparsMPV = new TF1("fparsMPV", direct_functions[1].c_str(), 0., fD_break);
	for(unsigned int i=0; i<direct_landauMPVpars.size(); ++i)       
	  fparsMPV->SetParameter(i, direct_landauMPVpars[i]);

	std::vector<double> direct_landauWidthpars = p.get<std::vector<double> >("Direct_landauWidthpars");
	fparsWidth = new TF1("fparsWidth", direct_functions[2].c_str(), 0., fD_break);
	for(unsigned int i=0; i<direct_landauWidthpars.size(); ++i)                
	  fparsWidth->SetParameter(i, direct_landauWidthpars[i]);

	std::vector<double> direct_expoCtepars = p.get<std::vector<double> >("Direct_expoCtepars");
	fparsCte = new TF1("fparsCte", direct_functions[3].c_str(), 0., fD_break);
	for(unsigned int i=0; i<direct_expoCtepars.size(); ++i)  
	  fparsCte->SetParameter(i, direct_expoCtepars[i]);

	std::vector<double> direct_expoSlopepars = p.get<std::vector<double> >("Direct_expoSlopepars");
	fparsSlope = new TF1("fparsSlope", direct_functions[4].c_str(), 0., fD_break);
	for(unsigned int i=0; i<direct_expoSlopepars.size(); ++i)
	  fparsSlope->SetParameter(i, direct_expoSlopepars[i]);

	std::vector<double> direct_landauNormpars_far = p.get<std::vector<double> >("Direct_landauNormpars_far");
        fparslogNorm_far = new TF1("fparslogNorm_far", direct_functions[5].c_str(), fD_break, fD_max);
        for(unsigned int i=0; i<direct_landauNormpars_far.size(); ++i)
          fparslogNorm_far->SetParameter(i, direct_landauNormpars_far[i]);

	std::vector<double> direct_landauMPVpars_far = p.get<std::vector<double> >("Direct_landauMPVpars_far");
        fparsMPV_far = new TF1("fparsMPV_far", direct_functions[6].c_str(), fD_break, fD_max);
        for(unsigned int i=0; i<direct_landauMPVpars_far.size(); ++i)
          fparsMPV_far->SetParameter(i, direct_landauMPVpars_far[i]);

	std::vector<double> direct_expoCtepars_far = p.get<std::vector<double> >("Direct_expoCtepars_far");
	fparsCte_far = new TF1("fparsCte_far", direct_functions[7].c_str(), fD_break - 50., fD_max);
	for(unsigned int i=0; i<direct_expoCtepars_far.size(); ++i)
          fparsCte_far->SetParameter(i, direct_expoCtepars_far[i]);

	std::vector<std::string> reflected_functions = p.get<std::vector<std::string> >("Reflected_functions");
        //times where the parametrizations are valid or change
	fT0_max = p.get<double>("T0_max");
	fT0_break_point = p.get<double>("T0_break_point");

	std::vector<double> reflected_landauNormpars = p.get<std::vector<double> >("Reflected_landauNormpars");
	fparslogNorm_refl = new TF1("fparslogNorm_refl", reflected_functions[0].c_str(), 0., fT0_max);
	for(unsigned int i=0; i<reflected_landauNormpars.size(); ++i)
	  fparslogNorm_refl->SetParameter(i, reflected_landauNormpars[i]);

	std::vector<double> reflected_landauMPVpars = p.get<std::vector<double> >("Reflected_landauMPVpars");
	fparsMPV_refl = new TF1("fparsMPV_refl", reflected_functions[1].c_str(), 0., fT0_max);
	for(unsigned int i=0; i<reflected_landauMPVpars.size(); ++i)
	  fparsMPV_refl->SetParameter(i, reflected_landauMPVpars[i]);

	std::vector<double> reflected_landauWidthpars = p.get<std::vector<double> >("Reflected_landauWidthpars");
	fparsWidth_refl = new TF1("fparsWidth_refl", reflected_functions[2].c_str(), 0., fT0_max);
	for(unsigned int i=0; i<reflected_landauWidthpars.size(); ++i)
	  fparsWidth_refl->SetParameter(i, reflected_landauWidthpars[i]);

	std::vector<double> reflected_expoCtepars = p.get<std::vector<double> >("Reflected_expoCtepars");
	fparsCte_refl = new TF1("fparsCte_refl", reflected_functions[3].c_str(), 0., fT0_max);
	for(unsigned int i=0; i<reflected_expoCtepars.size(); ++i)
	  fparsCte_refl->SetParameter(i, reflected_expoCtepars[i]);

	std::vector<double> reflected_expoSlopepars = p.get<std::vector<double> >("Reflected_expoSlopepars");
	fparsSlope_refl = new TF1("fparsSlope_refl", reflected_functions[4].c_str(), 0., fT0_max);
	for(unsigned int i=0; i<reflected_expoSlopepars.size(); ++i)
	  fparsSlope_refl->SetParameter(i, reflected_expoSlopepars[i]);
	*/

      }

    if(fUseNhitsModel) {
      	// VUV
      	fGH_PARS = p.get<std::vector<std::vector<double> > >("GH_PARS");
      	
        if (fStoreReflected) 
	{
        // VIS
      	fVIS_PARS = p.get<std::vector<std::vector<double>>>("VIS_PARS");
	fCATHODE_height = p.get<double>("CATHODE_height");
	fCATHODE_width = p.get<double>("CATHODE_width");
	fCATHODE_centre = p.get<std::vector<double>>("CATHODE_centre");
	fPlane_Depth = fCATHODE_centre[0];
	}
      	
	// Optical channel dimensions
	fOptical_Detector_Type = p.get<double>("Optical_Detector_Type");
      	fAPERTURE_height = p.get<double>("APERTURE_height");
      	fAPERTURE_width = p.get<double>("APERTURE_width");
      	fPMT_radius = p.get<double>("PMT_radius");
    }


    return;
	
  }



  //------------------------------------------------------

  // Eventually we will calculate the light quenching factor here
  double PhotonVisibilityService::GetQuenchingFactor(double /* dQdx */) const
  {
    // for now, no quenching
    return 1.0;

  }


  //------------------------------------------------------

  // Get a vector of the relative visibilities of each OpDet
  //  in the event to a point xyz

  float const* PhotonVisibilityService::GetAllVisibilities(double const* xyz, bool wantReflected) const
  {
    if(fInterpolate){
      static std::vector<float> ret;
      if(ret.size() != NOpChannels()) ret.resize(NOpChannels());
      for(unsigned int i = 0; i < NOpChannels(); ++i) ret[i] = GetVisibility(xyz, i, wantReflected);
      return &ret.front();
    }
    else{
      size_t VoxID = fVoxelDef.GetVoxelID(LibLocation(xyz));
      return GetLibraryEntries(VoxID, wantReflected);
    }
  }


  //------------------------------------------------------

  // Get distance to optical detector OpDet
  double PhotonVisibilityService::DistanceToOpDet( double const* xyz, unsigned int OpDet )
  {
    art::ServiceHandle<geo::Geometry> geom;
    return geom->OpDetGeoFromOpDet(OpDet).DistanceToPoint(xyz);
      
  }


  //------------------------------------------------------


  // Get the solid angle reduction factor for planar optical detector OpDet
  double PhotonVisibilityService::SolidAngleFactor( double const* xyz, unsigned int OpDet )
  {
    art::ServiceHandle<geo::Geometry> geom;
    return geom->OpDetGeoFromOpDet(OpDet).CosThetaFromNormal(xyz);
  }

  //------------------------------------------------------

  float PhotonVisibilityService::GetVisibility(double const* xyz, unsigned int OpChannel, bool wantReflected) const
  {
    // Static to avoid reallocating this buffer between calls
    // (GetNeighbouringVoxelIDs makes sure to clear it).
    static std::vector<sim::PhotonVoxelDef::NeiInfo> neis;

    if(fInterpolate){
      // In case we're outside the bounding box we'll get an empty vector here
      // and return visibility 0, which seems OK.
      fVoxelDef.GetNeighboringVoxelIDs(LibLocation(xyz), neis);
    }
    else{
      // For no interpolation, use a single entry with weight 1
      neis.clear();
      neis.emplace_back(fVoxelDef.GetVoxelID(LibLocation(xyz)), 1);
    }

    // Sum up all the weighted neighbours to get interpolation behaviour
    float vis = 0;
    for(const sim::PhotonVoxelDef::NeiInfo& n: neis){
      vis += n.weight * GetLibraryEntry(n.id, OpChannel, wantReflected);
    }
    return vis;
  }


  //------------------------------------------------------

  void PhotonVisibilityService::StoreLightProd(int VoxID, double N)
  {
    fCurrentVoxel = VoxID;
    fCurrentValue = N;
    mf::LogInfo("PhotonVisibilityService") << " PVS notes production of " << N << " photons at Vox " << VoxID<<std::endl; 
  }


  //------------------------------------------------------

  
  void PhotonVisibilityService::RetrieveLightProd(int& VoxID, double& N) const
  {
    N     = fCurrentValue;
    VoxID = fCurrentVoxel;
  }
  
  //------------------------------------------------------

  void PhotonVisibilityService::SetLibraryEntry(int VoxID, int OpChannel, float N, bool wantReflected)
  {   
    if(fTheLibrary == 0)
      LoadLibrary();
    
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);

    if(!wantReflected)
    lib->SetCount(VoxID,OpChannel, N);
    
    else
    lib->SetReflCount(VoxID,OpChannel, N);
    
    //std::cout<< " PVS logging " << VoxID << " " << OpChannel<<std::endl;
    MF_LOG_DEBUG("PhotonVisibilityService") << " PVS logging " << VoxID << " " << OpChannel<<std::endl;
  }

  //------------------------------------------------------

  float const* PhotonVisibilityService::GetLibraryEntries(int VoxID, bool wantReflected) const
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    if(!wantReflected)
      return fTheLibrary->GetCounts(VoxID);
    else
      return fTheLibrary->GetReflCounts(VoxID);
  }

  //------------------------------------------------------

  float PhotonVisibilityService::GetLibraryEntry(int VoxID, int Channel, bool wantReflected) const
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    if(!wantReflected)
      return fTheLibrary->GetCount(VoxID, Channel);
    else
      return fTheLibrary->GetReflCount(VoxID, Channel);
  }
  // Methodes to handle the extra library parameter refl T0
  //------------------------------------------------------

  // Get a vector of the refl <tfirst> of each OpDet
  //  in the event to a point xyz

  float const* PhotonVisibilityService::GetReflT0s(double const* xyz) const
  {
    int VoxID = fVoxelDef.GetVoxelID(LibLocation(xyz));
    return GetLibraryReflT0Entries(VoxID);
  }

  //------------------------------------------------------

  float const* PhotonVisibilityService::GetLibraryReflT0Entries(int VoxID) const
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    return fTheLibrary->GetReflT0s(VoxID);
  }

  //------------------------------------------------------     

  void PhotonVisibilityService::SetLibraryReflT0Entry(int VoxID, int OpChannel, float T0)
  {
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
    if(fTheLibrary == 0)
      LoadLibrary();

    lib->SetReflT0(VoxID,OpChannel,T0);

    MF_LOG_DEBUG("PhotonVisibilityService") << " PVS logging " << VoxID << " " << OpChannel<<std::endl;
  }

  //------------------------------------------------------      

  float PhotonVisibilityService::GetLibraryReflT0Entry(int VoxID, int Channel) const
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    return fTheLibrary->GetReflT0(VoxID, Channel);
  }

  //------------------------------------------------------


/////////////****////////////

  const std::vector<float>* PhotonVisibilityService::GetTimingPar(double const* xyz) const
  {
    int VoxID = fVoxelDef.GetVoxelID(LibLocation(xyz));
    return GetLibraryTimingParEntries(VoxID);
  }

  TF1* PhotonVisibilityService::GetTimingTF1(double const* xyz) const
  {
    int VoxID = fVoxelDef.GetVoxelID(LibLocation(xyz));
    return GetLibraryTimingTF1Entries(VoxID);
  }


  //------------------------------------------------------

  const std::vector<float>* PhotonVisibilityService::GetLibraryTimingParEntries(int VoxID) const
  {
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
    if(fTheLibrary == 0)
      LoadLibrary();

    return lib->GetTimingPars(VoxID);
  }

  //------------------------------------------------------

  TF1* PhotonVisibilityService::GetLibraryTimingTF1Entries(int VoxID) const
  {
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
    if(fTheLibrary == 0)
      LoadLibrary();

    return lib->GetTimingTF1s(VoxID);
  }

  //------------------------------------------------------     

  void PhotonVisibilityService::SetLibraryTimingParEntry(int VoxID, int OpChannel, float par, size_t parnum)
  {
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
    if(fTheLibrary == 0)
      LoadLibrary();

    lib->SetTimingPar(VoxID,OpChannel,par, parnum);

    MF_LOG_DEBUG("PhotonVisibilityService") << " PVS logging " << VoxID << " " << OpChannel<<std::endl;
  }

  //------------------------------------------------------

  void PhotonVisibilityService::SetLibraryTimingTF1Entry(int VoxID, int OpChannel, TF1 func)
  {
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
    if(fTheLibrary == 0)
      LoadLibrary();

    lib->SetTimingTF1(VoxID,OpChannel,func);

    MF_LOG_DEBUG("PhotonVisibilityService") << " PVS logging " << VoxID << " " << OpChannel<<std::endl;
  }


  //------------------------------------------------------

  float PhotonVisibilityService::GetLibraryTimingParEntry(int VoxID, int Channel, size_t npar) const
  {
    PhotonLibrary* lib = dynamic_cast<PhotonLibrary*>(fTheLibrary);
    if(fTheLibrary == 0)
      LoadLibrary();

    return lib->GetTimingPar(VoxID, Channel,npar);
  }

  //------------------------------------------------------

  size_t PhotonVisibilityService::NOpChannels() const
  {
    if(fTheLibrary == 0)
      LoadLibrary();
    
    return fTheLibrary->NOpChannels();
  }

  //------------------------------------------------------
  void PhotonVisibilityService::SetDirectLightPropFunctions(TF1 const* functions[8], double& d_break, double& d_max, double& tf1_sampling_factor) const
  {
    functions[0] = fparslogNorm;
    functions[1] = fparsMPV;
    functions[2] = fparsWidth;
    functions[3] = fparsCte;
    functions[4] = fparsSlope; 
    functions[5] = fparslogNorm_far;
    functions[6] = fparsMPV_far;
    functions[7] = fparsCte_far;

    d_break = fD_break;
    d_max = fD_max;
    tf1_sampling_factor = fTF1_sampling_factor;
  }

  //------------------------------------------------------
  void PhotonVisibilityService::SetReflectedCOLightPropFunctions(TF1 const* functions[5], double& t0_max, double& t0_break_point) const
  {
    functions[0] = fparslogNorm_refl;
    functions[1] = fparsMPV_refl;
    functions[2] = fparsWidth_refl;
    functions[3] = fparsCte_refl;
    functions[4] = fparsSlope_refl;

    t0_max = fT0_max;
    t0_break_point = fT0_break_point;
  }



  //------------------------------------------------------                                                                             
  void PhotonVisibilityService::LoadTimingsForVUVPar(std::vector<double> v[9], double& step_size, double& max_d, 
						     double& vuv_vgroup_mean, double& vuv_vgroup_max, double& inflexion_point_distance) const
  {
    v[0] = fDistances_all;
    v[1] = fNorm_over_entries;
    v[2] = fMpv;
    v[3] = fWidth;
    v[4] = fDistances;
    v[5] = fSlope;
    v[6] = fExpo_over_Landau_norm[0];
    v[7] = fExpo_over_Landau_norm[1];
    v[8] = fExpo_over_Landau_norm[2];
    
    step_size = fstep_size;
    max_d = fmax_d;
    vuv_vgroup_mean = fvuv_vgroup_mean;
    vuv_vgroup_max = fvuv_vgroup_max;
    inflexion_point_distance = finflexion_point_distance;

  }

  void PhotonVisibilityService::LoadTimingsForVISPar(std::vector<double>& distances, std::vector<std::vector<double>>& cut_off, std::vector<std::vector<double>>& tau,
										 double& vis_vmean, double& n_vis, double& n_vuv, double& plane_depth) const
  {
  distances = fDistances_refl;
  cut_off = fCut_off;
  tau = fTau;

  vis_vmean = fvis_vmean;
  n_vis = fn_LAr_vis;
  n_vuv = fn_LAr_VUV;
  plane_depth = fPlane_Depth;

  }


  void PhotonVisibilityService::LoadGHForVUVCorrection(std::vector<std::vector<double>>& v, double& w, double& h, double& r, int& op_det_type) const
  {
    v = fGH_PARS;    
 
    op_det_type = fOptical_Detector_Type;
    h = fAPERTURE_height;
    w = fAPERTURE_width;
    r = fPMT_radius;

  }

  void PhotonVisibilityService::LoadParsForVISCorrection(std::vector<std::vector<double>>& v, double& plane_depth, double& w_cathode, double& h_cathode, std::vector<double>& cntr_cathode, double& w, double& h, double& r, int& op_det_type) const
  {
    v = fVIS_PARS;
    plane_depth = fPlane_Depth;
    w_cathode = fCATHODE_width;
    h_cathode = fCATHODE_height;
    cntr_cathode = fCATHODE_centre;

    op_det_type = fOptical_Detector_Type;
    h = fAPERTURE_height;
    w = fAPERTURE_width;
    r = fPMT_radius;

  }


  //------------------------------------------------------
  /***
   * Preform any necessary transformations on the coordinates before trying to access
   * a voxel ID.
   **/
  const TVector3 PhotonVisibilityService::LibLocation(const double * xyz) const
  {
    TVector3 location(xyz);
    
    // Always use postive X coordinate if set
    if (fReflectOverZeroX && location.x() < 0) {
      location.SetX( fabs(location.x() ) );
    }
    return location;
  }

} // namespace

namespace phot{
 
  DEFINE_ART_SERVICE(PhotonVisibilityService)

} // namespace phot
