////////////////////////////////////////////////////////////////////////
/// \file OpDetSensitiveDetector.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This is the sensitive detector class for the OpDet detectors.  It is
// associated with the relevant detector volumes in DetectorConstruction.cxx
// and is called via the ProcessHits method every time a particle steps
// within the volume.
//
// The detector owns a hit collection which is passed back to LArG4 at
// the end of the event.  One OpDetSensitiveDetector corresponds to a set
// of OpDet's, which are looked up by their G4PhysicalVolume in the OpDetLookup class.
//
// Photons stepping into the volume are stopped and killed and their trackID,
// 4 position and 4 momentum are stored in the relevant SimPhotons.
//
// Ben Jones, MIT, 06/04/2010
//

#ifndef OpDetSensitiveDetector_h
#define OpDetSensitiveDetector_h 1

#include "Geant4/G4String.hh"
#include "Geant4/G4Types.hh"
#include "Geant4/G4VSensitiveDetector.hh"
class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

namespace larg4 {

  class OpDetLookup;
  class OpDetPhotonTable;

  class OpDetSensitiveDetector : public G4VSensitiveDetector
  {


  public:
    OpDetSensitiveDetector(G4String name);
    virtual ~OpDetSensitiveDetector(){}


    // Beginning and end of event
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*){}

    // Tidy up event in abort
    virtual void clear(){}

    // Run per step in sensitive volume
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory*);


    // Required but empty
    virtual void DrawAll(){}
    virtual void PrintAll(){}

  private:
    OpDetLookup              * fTheOpDetLookup;
    OpDetPhotonTable         * fThePhotonTable;

    //double                     fGlobalTimeOffset;
  };
}

#endif
