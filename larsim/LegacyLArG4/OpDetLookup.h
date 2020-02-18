////////////////////////////////////////////////////////////////////////
/// \file OpDetLookup.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Provide a map between G4VPhysicalVolumes of OpDets and OpDet Numbers's.
//
// There are two places where optical detectors must be known about
// in larsoft: In the Geant4 volume store, where they are associated
// to G4SensitiveDetectors, and in the Geometry service, where their
// positions and cryostat associations are known for reconstruction.
//
// These geometries are built independently, and the OpDet's in each
// volume is not provided by the geometry specification.
//
// The main function of this class is to provide a link between the
// unlabelled G4PhysicalVolumes, and the OpDet objects organized into
// vectors in the Geomtry/CryostatGeo objects, accessible through the
// geo::Geometry service.
//
// Any physical volume in the gdml which has the specified opdet
// name (specified in the geometry service) is given a sensitive detector.
// This sensitive detector passes the volume to this service, which
// determines based on its position which element in the OpDetGeo collection
// it must correspond to.
//
// It is then renamed accordingly and the link between the two objects
// is stored in a map<string, int> which relates the new G4 name
// to a detector number in the geometry.
//
//
// Ben Jones, MIT, 06/04/2010
//

#ifndef OpDetLOOKUP_h
#define OpDetLOOKUP_h 1

#include <map>
#include <string>

class G4VPhysicalVolume;

namespace larg4 {
  class OpDetLookup
    {
    public:
      ~OpDetLookup(){}
      static OpDetLookup * Instance();

      void AddPhysicalVolume(G4VPhysicalVolume *);
      int GetOpDet(G4VPhysicalVolume *);
      int GetOpDet(std::string);
      int GetN();
      int FindClosestOpDet(G4VPhysicalVolume* vol,  double& Distance);

    protected:
      OpDetLookup();

    private:
      std::map<std::string, int> fTheOpDetMap;
      int fTheTopOpDet;

    };

}


#endif
