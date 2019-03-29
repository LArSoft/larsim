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

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/PhotonVoxels.h"
#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump::vector3D()

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

namespace phot {

  class PhotonLibraryAnalyzer : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit PhotonLibraryAnalyzer(fhicl::ParameterSet const& pset);

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);
    
  private:
    std::string fAltXAxis;
    int         fOpDet;
    bool        fEachSlice;
    bool        fEachDetector;
  };
  
}

namespace phot {

  //----------------------------------------------------------------------------
  PhotonLibraryAnalyzer::PhotonLibraryAnalyzer(const fhicl::ParameterSet& pset)
    : EDAnalyzer(pset)
    , fAltXAxis{pset.get<std::string>("alt_x_axis")}
    , fOpDet{pset.get<int>("opdet")}
    , fEachSlice{pset.get<bool>("each_slice")}
    , fEachDetector{pset.get<bool>("each_detector")}
  {
    std::cout<<"Photon library analyzer constructor "<<std::endl;
  }

  //----------------------------------------------------------------------------
  void PhotonLibraryAnalyzer::beginJob()
  
  {
    mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - begin"<< std::endl;

    
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<PhotonVisibilityService> pvs;
    art::ServiceHandle<geo::Geometry> geom;

    int NOpDet = pvs->NOpChannels();

    sim::PhotonVoxelDef TheVoxelDef = pvs->GetVoxelDef();
    auto const& UpperCorner = TheVoxelDef.GetRegionUpperCorner();
    auto const& LowerCorner = TheVoxelDef.GetRegionLowerCorner();

    mf::LogInfo("PhotonLibraryAnalyzer") << "UpperCorner: " << lar::dump::vector3D(UpperCorner) << "\n"
                                         << "LowerCorner: " << lar::dump::vector3D(LowerCorner);

    auto const [ XSteps, YSteps, ZSteps ] = TheVoxelDef.GetSteps(); // unsigned int

    // for c2: FullVolume is unused, just call tfs->make
    // TH3D *FullVolume = tfs->make<TH3D>("FullVolume","FullVolume", 
    tfs->make<TH3D>("FullVolume","FullVolume", 
                    XSteps,LowerCorner.X(),UpperCorner.X(),
                    YSteps,LowerCorner.Y(),UpperCorner.Y(),
                    ZSteps,LowerCorner.Z(),UpperCorner.Z());

    
    int reportnum=10000;

    int newX, newY;
    if (fAltXAxis == "Z") { 
      newX = 2; // Z
      newY = 1; // Y
    }
    else {
      newX = 1; // Y
      newY = 2; // Z
    }

 
    
    mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - making historams"<< std::endl;

    TH2D* XProjection;
    if (fAltXAxis == "Z") XProjection = tfs->make<TH2D>("XProjection","XProjection",ZSteps,0,ZSteps,YSteps,0,YSteps);
    else                  XProjection = tfs->make<TH2D>("XProjection","XProjection",YSteps,0,YSteps,ZSteps,0,ZSteps);
    TH2D* YProjection = tfs->make<TH2D>("YProjection","YProjection",XSteps,0,XSteps,ZSteps,0,ZSteps);
    TH2D* ZProjection = tfs->make<TH2D>("ZProjection","ZProjection",XSteps,0,XSteps,YSteps,0,YSteps);
   
    //    TH1D * PMTsNoVisibility = tfs->make<TH1D>("PMTsNoVisibility","PMTsNoVisibility", NOpDet,0,NOpDet); 

    TH1D* VisByN = tfs->make<TH1D>("VisByN","VisByN", NOpDet, 0, NOpDet);

    TH2D* XInvisibles;
    if (fAltXAxis == "Z") XInvisibles = tfs->make<TH2D>("XInvisibles","XInvisibles",ZSteps,0,ZSteps,YSteps,0,YSteps);
    else                  XInvisibles = tfs->make<TH2D>("XInvisibles","XInvisibles",YSteps,0,YSteps,ZSteps,0,ZSteps);
    TH2D* YInvisibles = tfs->make<TH2D>("YInvisibles","YInvisibles",XSteps,0,XSteps,ZSteps,0,ZSteps);
    TH2D* ZInvisibles = tfs->make<TH2D>("ZInvisibles","ZInvisibles",XSteps,0,XSteps,YSteps,0,YSteps);

    

    

    std::vector<TH2D*> TheXCrossSections;
    std::vector<TH2D*> TheYCrossSections;
    std::vector<TH2D*> TheZCrossSections;
    if (fEachSlice) {

    
      for(unsigned int i=0; i!=XSteps; ++i)
      {
        std::stringstream ss("");
        ss.flush();
        ss<<"projX"<<i;
        if (fAltXAxis == "Z") 	
          TheXCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), ZSteps, 0,ZSteps, YSteps, 0,YSteps));
        else	
          TheXCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), YSteps, 0,YSteps, ZSteps, 0,ZSteps));                  


      }

      for(unsigned int i=0; i!=YSteps; ++i)
      {
        std::stringstream ss("");
        ss.flush();
        ss<<"projY"<<i;
        TheYCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), XSteps, 0,XSteps, ZSteps, 0, ZSteps));

      }

      for(unsigned int i=0; i!=ZSteps; ++i)
      {
        std::stringstream ss("");
        ss.flush();
        ss<<"projZ"<<i;
        TheZCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), XSteps, 0,XSteps, YSteps, 0,YSteps));

      }
    }


    std::vector<TH2D*> TheXProjections;
    std::vector<TH2D*> TheYProjections;
    std::vector<TH2D*> TheZProjections;

    if (fEachDetector) {

      mf::LogInfo("PhotonLibraryAnalyzer")<<"Making projections for each of " << NOpDet << " photon detectors" << std::endl;
    
      for(int i=0; i<NOpDet; ++i)
      {
        char ss[99];

        sprintf(ss, "ProjXOpDet%d", i);
        if (fAltXAxis == "Z")
          TheXProjections.push_back(tfs->make<TH2D>(ss, ss, ZSteps, 0,ZSteps, YSteps, 0,YSteps));
        else
          TheXProjections.push_back(tfs->make<TH2D>(ss, ss, YSteps, 0,YSteps, ZSteps, 0,ZSteps));

        sprintf(ss, "ProjYOpDet%d", i);
        TheYProjections.push_back(tfs->make<TH2D>(ss, ss, XSteps, 0,XSteps, ZSteps, 0, ZSteps));

        sprintf(ss, "ProjZOpDet%d", i);
        TheZProjections.push_back(tfs->make<TH2D>(ss, ss, XSteps, 0,XSteps, YSteps, 0,YSteps));
      }
    }

    
    mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - running through voxels "<< std::endl;


    for(unsigned int i=0; i!=TheVoxelDef.GetNVoxels(); ++i)
    {
      if(i%reportnum==0) std::cout<<"Photon library analyzer at voxel " << i<<std::endl;
    
      auto const Coords = TheVoxelDef.GetVoxelCoords(i);
                
      const float* Visibilities = pvs->GetLibraryEntries(i);
      size_t NOpChannels = pvs->NOpChannels();
      
      float TotalVis=0;
      if (fOpDet < 0) {
        for(size_t ichan=0; ichan!=NOpChannels; ++ichan)
        {
          TotalVis+=Visibilities[ichan];	
        }
      }
      else {
        TotalVis = Visibilities[fOpDet];
      }
      
      VisByN->Fill(NOpChannels);
      
      if(TotalVis==0)
      {
        XInvisibles->Fill(Coords[newX],Coords[newY]);
        YInvisibles->Fill(Coords[0],Coords[2]);
        ZInvisibles->Fill(Coords[0],Coords[1]);
      }

      if (fEachSlice) {
        TheXCrossSections.at(Coords.at(0))->Fill(Coords[newX],Coords[newY],TotalVis);
        TheYCrossSections.at(Coords.at(1))->Fill(Coords[0],Coords[2],TotalVis);
        TheZCrossSections.at(Coords.at(2))->Fill(Coords[0],Coords[1],TotalVis);
      }

      if (fEachDetector) {
        for(size_t ichan=0; ichan!=NOpChannels; ++ichan) {
          TheXProjections.at(ichan)->Fill(Coords[newX],Coords[newY],Visibilities[ichan]);
          TheYProjections.at(ichan)->Fill(Coords[0],Coords[2],Visibilities[ichan]);
          TheZProjections.at(ichan)->Fill(Coords[0],Coords[1],Visibilities[ichan]);
        }
      }

      // Always make the summed projections
      XProjection->Fill(Coords[newX], Coords[newY], TotalVis);
      YProjection->Fill(Coords[0], Coords[2], TotalVis);
      ZProjection->Fill(Coords[0], Coords[1], TotalVis);
     
    }
  
    mf::LogInfo("PhotonLibraryAnalyzer")<<"Analyzing photon library - end"<< std::endl;
  }

  //----------------------------------------------------------------------------
  void PhotonLibraryAnalyzer::analyze(const art::Event& /*evt*/)
  {
  
  }

}


namespace phot {
  DEFINE_ART_MODULE(PhotonLibraryAnalyzer)
}
