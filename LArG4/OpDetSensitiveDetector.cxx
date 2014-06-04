////////////////////////////////////////////////////////////////////////
/// \file OpDetSensitiveDetector.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpDetSensitiveDetector
//
// See comments in OpDetSensitiveDetector.h
//
// Ben Jones, MIT, 06/04/2010
//



#include "LArG4/OpDetSensitiveDetector.h"
#include "LArG4/OpDetPhotonTable.h"
#include "LArG4/OpDetLookup.h"
#include "Simulation/SimPhotons.h"
#include "Geant4/G4SDManager.hh"
#include "Utilities/DetectorProperties.h"

namespace larg4{


  OpDetSensitiveDetector::OpDetSensitiveDetector(G4String DetectorUniqueName)
    : G4VSensitiveDetector(DetectorUniqueName)
  {
    // Register self with sensitive detector manager
    G4SDManager::GetSDMpointer()->AddNewDetector(this);

    // Get instances of singleton classes
    fTheOpDetLookup        = OpDetLookup::Instance();
    fThePhotonTable        = OpDetPhotonTable::Instance();
    
    art::ServiceHandle<util::DetectorProperties> det;
    //fGlobalTimeOffset = det->ConvertTicksToTDC(0) * det->SamplingRate();
  }

  
  //--------------------------------------------------------

  G4bool OpDetSensitiveDetector::ProcessHits(G4Step * aStep, G4TouchableHistory *)
  {
    sim::OnePhoton ThePhoton;
    
    
    // Get photon data to store in the hit

    ThePhoton.SetInSD      = true;

    ThePhoton.InitialPosition     = TVector3(
					      aStep->GetTrack()->GetVertexPosition().x(),
					      aStep->GetTrack()->GetVertexPosition().y(),
					      aStep->GetTrack()->GetVertexPosition().z()	
					      );
    
    //ThePhoton.Time                = aStep->GetTrack()->GetGlobalTime() - fGlobalTimeOffset;
    ThePhoton.Time                = aStep->GetTrack()->GetGlobalTime();
    

    ThePhoton.Energy              =  aStep->GetTrack()->GetVertexKineticEnergy();
    
    // Lookup which OpDet we are in

    int OpDetChannel = fTheOpDetLookup->GetChannel(aStep->GetPreStepPoint()->GetPhysicalVolume());

    // Add this photon to the detected photons table
    fThePhotonTable->AddPhoton(OpDetChannel, std::move(ThePhoton));
    
    // Kill this photon track
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    return true;
    

  }

  //--------------------------------------------------------
  
  void OpDetSensitiveDetector::Initialize(G4HCofThisEvent *)
  {

  }

}
