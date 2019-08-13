////////////////////////////////////////////////////////////////////////
/// \file OpDetReadoutGeometry.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// OpDetReadoutGeometry
//
// Ben Jones, MIT, 07/16/10
//
// OpDet's are defined with a particular geometry, and only some part of the OpDet
// is sensitive.  When several are placed, only the mother volume, vol_OpDet is
// replicated.  Each constituent volume, including the sensitive volume, only
// exists once in the volume store, as a daughter of the mother volume.
//
// Hence to know which OpDet a photon steps into when it is inside a sensitive
// volume, we must define a readout geometry to identify which volumes are
// contained within each OpDet.
//
// This class is heavily based on LArG4ReadoutGeometry by Bill Seligman,
// which is very well commented.  See that file for further reference.
//

#ifndef OpDetReadoutGeometry_h
#define OpDetReadoutGeometry_h

#include "Geant4/G4String.hh"
#include "Geant4/G4Transform3D.hh"
#include "Geant4/G4VUserParallelWorld.hh"

#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;

namespace larg4
{

  class OpDetReadoutGeometry : public G4VUserParallelWorld
    {
    public:
      OpDetReadoutGeometry(G4String OpDetSensitiveName, const G4String name = "OpDetReadoutGeometry");
      virtual ~OpDetReadoutGeometry();

      virtual void Construct();
    private:
      void                            FindVolumes(G4VPhysicalVolume *, G4String, std::vector<G4Transform3D>, std::vector<G4LogicalVolume*>&, std::vector<G4Transform3D>&);
      std::vector<G4LogicalVolume*>   fOpDetVolumes;
      std::vector<G4Transform3D>      fOpDetTransformations;
      G4String fOpDetSensitiveName;

    };

}

#endif
