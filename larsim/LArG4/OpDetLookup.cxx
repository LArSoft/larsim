////////////////////////////////////////////////////////////////////////
/// \file OpDetLookup.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpDetLookup class.
//
// See comments in the OpDetLookup.h file.
//
// Ben Jones, MIT, 06/04/2010
//


#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Geant4/G4VPhysicalVolume.hh"

#include "larsim/LArG4/OpDetLookup.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/OpDetGeo.h"


namespace larg4 {
  OpDetLookup * TheOpDetLookup;

  //--------------------------------------------------
  OpDetLookup::OpDetLookup()
  {
    fTheTopOpDet=0;
  }

  //--------------------------------------------------
  OpDetLookup * OpDetLookup::Instance()
  {
    if(!TheOpDetLookup){
      TheOpDetLookup = new OpDetLookup;
    }
    return TheOpDetLookup;
  }

  //--------------------------------------------------
  int OpDetLookup::GetOpDet(std::string TheName)
  {
    return fTheOpDetMap[TheName];
  }

  //--------------------------------------------------
  int OpDetLookup::GetOpDet(G4VPhysicalVolume* TheVolume)
  {
    std::string TheName = TheVolume->GetName();
    return GetOpDet(TheName);
  }


  //--------------------------------------------------

  int OpDetLookup::FindClosestOpDet(G4VPhysicalVolume* vol, double& distance)
  {
    art::ServiceHandle<geo::Geometry const> geom;
    int    OpDetCount = 0;

    double MinDistance = UINT_MAX;
    int    ClosestOpDet   = -1;

    for(size_t o=0; o!=geom->NOpDets(); o++) {
      double xyz[3];
      geom->OpDetGeoFromOpDet(o).GetCenter(xyz);

      CLHEP::Hep3Vector DetPos(xyz[0],xyz[1],xyz[2]);
      CLHEP::Hep3Vector ThisVolPos = vol->GetTranslation();

      ThisVolPos/=CLHEP::cm;

      //	    std::cout<<"Det: " << xyz[0]<< " " <<xyz[1]<< " " << xyz[2]<<std::endl;
      //    std::cout<<"Vol: " << ThisVolPos.x()<< " " <<ThisVolPos.y() << " " <<ThisVolPos.z()<<std::endl;

      double Distance = (DetPos-ThisVolPos).mag();
      if(Distance < MinDistance)
      {
        MinDistance = Distance;
        ClosestOpDet  =  o;
      }
      OpDetCount++;
    }
    if(ClosestOpDet<0)
      {
	throw cet::exception("OpDetLookup Error") << "No nearby OpDet found!\n";
      }

    distance = MinDistance;
    return ClosestOpDet;
  }


  //--------------------------------------------------
  void OpDetLookup::AddPhysicalVolume(G4VPhysicalVolume * volume)
  {

    // mf::LogInfo("Optical") <<"G4 placing sensitive opdet"<<std::endl;

    std::stringstream VolName("");
    double Distance     = 0;

    int NearestOpDet = FindClosestOpDet(volume, Distance);

    VolName.flush();
    VolName << volume->GetName() << "_" << NearestOpDet;
    volume->SetName(VolName.str().c_str());

    fTheOpDetMap[VolName.str()] = NearestOpDet;

    // mf::LogInfo("Optical") << "Found closest volume: " << VolName.str().c_str() << " OpDet : " << fTheOpDetMap[VolName.str()]<<"  distance : " <<Distance<<std::endl;

  }


  //--------------------------------------------------
  int OpDetLookup::GetN()
  {
    return fTheTopOpDet;
  }

}
