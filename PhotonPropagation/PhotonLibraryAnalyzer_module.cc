#ifndef PHOTLIBANALYZER_H
#define PHOTLIBANALYZER_H

//
// Name: PhotonLibraryAnalyzer.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//
#include <iostream>
#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 

#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "Simulation/PhotonVoxels.h"

#include "TH1D.h"
#include "TH2D.h"

namespace phot {

  class PhotonLibraryAnalyzer : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit PhotonLibraryAnalyzer(fhicl::ParameterSet const& pset);
    virtual ~PhotonLibraryAnalyzer();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void analyze(const art::Event& evt);
    
    void endJob();
  
    

  private:


  };
  
}

#endif

namespace phot {

  //----------------------------------------------------------------------------
  PhotonLibraryAnalyzer::PhotonLibraryAnalyzer(const fhicl::ParameterSet& pset)
    : EDAnalyzer(pset)
  {
    reconfigure(pset);
  }

  //----------------------------------------------------------------------------
  PhotonLibraryAnalyzer::~PhotonLibraryAnalyzer()
  {
  }

  //----------------------------------------------------------------------------
  void PhotonLibraryAnalyzer::reconfigure(fhicl::ParameterSet const& pset)
  {
    
  }

  //----------------------------------------------------------------------------
  void PhotonLibraryAnalyzer::beginJob()
  
  {
    mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - begin"<<
      std::endl;
    
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<PhotonVisibilityService> pvs;
    art::ServiceHandle<geo::Geometry> geom;

    int NOpDet = geom->NOpDet();

    sim::PhotonVoxelDef TheVoxelDef = pvs->GetVoxelDef();
    TVector3 Steps = TheVoxelDef.GetSteps();
    TVector3 UpperCorner = TheVoxelDef.GetRegionUpperCorner();
    TVector3 LowerCorner = TheVoxelDef.GetRegionLowerCorner();
 
    int XSteps = int(Steps[0]);
    int YSteps = int(Steps[1]);
    int ZSteps = int(Steps[2]);
 
    
    mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - making historams"<<
      std::endl;
    
    TH2D* XProjection = tfs->make<TH2D>("XProjection","XProjection",YSteps,0,YSteps,ZSteps,0,ZSteps);
    TH2D* YProjection = tfs->make<TH2D>("YProjection","YProjection",XSteps,0,XSteps,ZSteps,0,ZSteps);
    TH2D* ZProjection = tfs->make<TH2D>("ZProjection","ZProjection",XSteps,0,XSteps,YSteps,0,YSteps);
   
    //    TH1D * PMTsNoVisibility = tfs->make<TH1D>("PMTsNoVisibility","PMTsNoVisibility", NOpDet,0,NOpDet); 

    TH1D* VisByN = tfs->make<TH1D>("VisByN","VisByN", NOpDet, 0, NOpDet);

    TH2D* XInvisibles = tfs->make<TH2D>("XInvisibles","XInvisibles",YSteps,0,YSteps,ZSteps,0,ZSteps);
    TH2D* YInvisibles = tfs->make<TH2D>("YInvisibles","YInvisibles",XSteps,0,XSteps,ZSteps,0,ZSteps);
    TH2D* ZInvisibles = tfs->make<TH2D>("ZInvisibles","ZInvisibles",XSteps,0,XSteps,YSteps,0,YSteps);


  
    std::vector<TH2D*> TheXCrossSections;
    std::vector<TH2D*> TheYCrossSections;
    std::vector<TH2D*> TheZCrossSections;

    
    for(int i=0; i!=XSteps; ++i)
      {
	std::stringstream ss("");
	ss.flush();
	ss<<"projX"<<i;
	TheXCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), YSteps, 0,YSteps, ZSteps, 0,ZSteps));

      }

  for(int i=0; i!=YSteps; ++i)
      {
	std::stringstream ss("");
	ss.flush();
	ss<<"projY"<<i;
	TheYCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), XSteps, 0,XSteps, ZSteps, 0, ZSteps));

      }

  for(int i=0; i!=ZSteps; ++i)
      {
	std::stringstream ss("");
	ss.flush();
	ss<<"projZ"<<i;
	TheZCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), XSteps, 0,XSteps, YSteps, 0,YSteps));

      }
    
  mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - running through voxels "<<
    std::endl;

  int reportnum=10000;

  for(int i=0; i!=TheVoxelDef.GetNVoxels(); ++i)
    {
      if(i%reportnum==0) std::cout<<"Photon library analyzer at voxel " << i<<std::endl;
    
      std::vector<int> Coords = TheVoxelDef.GetVoxelCoords(i);
                
      std::vector<float>* Visibilities = pvs->GetLibraryEntries(i);
      
      float TotalVis=0;
      for(size_t ichan=0; ichan!=Visibilities->size(); ++ichan)
      {
	TotalVis+=Visibilities->at(ichan);	
      }
      
      VisByN->Fill(Visibilities->size());
      
      if(TotalVis==0)
	{
	  XInvisibles->Fill(Coords[1],Coords[2]);
	  YInvisibles->Fill(Coords[0],Coords[2]);
	  ZInvisibles->Fill(Coords[0],Coords[1]);
	}

      TheXCrossSections.at(Coords.at(0))->Fill(Coords[1],Coords[2],TotalVis);
      TheYCrossSections.at(Coords.at(1))->Fill(Coords[0],Coords[2],TotalVis);
      TheZCrossSections.at(Coords.at(2))->Fill(Coords[0],Coords[1],TotalVis);
      XProjection->Fill(Coords[1], Coords[2], TotalVis);
      YProjection->Fill(Coords[0], Coords[2], TotalVis);
      ZProjection->Fill(Coords[0], Coords[1], TotalVis);
     
    }
  
  mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - end"<<
    std::endl;
  }

  //----------------------------------------------------------------------------
  void PhotonLibraryAnalyzer::analyze(const art::Event& evt)
  {
  
  }



  
  //----------------------------------------------------------------------------
  void PhotonLibraryAnalyzer::endJob()
  {

  }
}


namespace phot {
  DEFINE_ART_MODULE(PhotonLibraryAnalyzer)
}
