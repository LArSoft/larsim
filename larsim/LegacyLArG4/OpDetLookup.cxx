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

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larsim/LegacyLArG4/OpDetLookup.h"

namespace larg4 {
  OpDetLookup* TheOpDetLookup;

  //--------------------------------------------------
  OpDetLookup::OpDetLookup() { fTheTopOpDet = 0; }

  //--------------------------------------------------
  OpDetLookup* OpDetLookup::Instance()
  {
    if (!TheOpDetLookup) { TheOpDetLookup = new OpDetLookup; }
    return TheOpDetLookup;
  }

  //--------------------------------------------------
  int OpDetLookup::GetOpDet(std::string const& TheName) { return fTheOpDetMap[TheName]; }

  //--------------------------------------------------
  int OpDetLookup::GetOpDet(G4VPhysicalVolume const* TheVolume)
  {
    return GetOpDet(TheVolume->GetName());
  }

  //--------------------------------------------------

  int OpDetLookup::FindClosestOpDet(G4VPhysicalVolume* vol, double& distance)
  {
    art::ServiceHandle<geo::Geometry const> geom;
    int OpDetCount = 0;

    double MinDistance = UINT_MAX;
    int ClosestOpDet = -1;

    for (size_t o = 0; o != geom->NOpDets(); o++) {
      auto const xyz = geom->OpDetGeoFromOpDet(o).GetCenter();

      CLHEP::Hep3Vector DetPos(xyz.X(), xyz.Y(), xyz.Z());
      CLHEP::Hep3Vector ThisVolPos = vol->GetTranslation();

      ThisVolPos /= CLHEP::cm;

      double Distance = (DetPos - ThisVolPos).mag();
      if (Distance < MinDistance) {
        MinDistance = Distance;
        ClosestOpDet = o;
      }
      OpDetCount++;
    }
    if (ClosestOpDet < 0) {
      throw cet::exception("OpDetLookup Error") << "No nearby OpDet found!\n";
    }

    distance = MinDistance;
    return ClosestOpDet;
  }

  //--------------------------------------------------
  void OpDetLookup::AddPhysicalVolume(G4VPhysicalVolume* volume)
  {
    std::stringstream VolName("");
    double Distance = 0;

    int NearestOpDet = FindClosestOpDet(volume, Distance);

    VolName.flush();
    VolName << volume->GetName() << "_" << NearestOpDet;
    volume->SetName(VolName.str().c_str());

    fTheOpDetMap[VolName.str()] = NearestOpDet;
  }

  //--------------------------------------------------
  int OpDetLookup::GetN() const { return fTheTopOpDet; }

}
