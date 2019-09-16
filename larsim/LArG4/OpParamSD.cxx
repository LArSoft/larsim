////////////////////////////////////////////////////////////////////////
/// \file OpParamSD.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpParamSD
//
// See comments in OpParamSD.h
//
// Ben Jones, MIT, 2013
//

#include "cetlib_except/exception.h"
#include "larsim/LArG4/OpParamAction.h"
#include "larsim/LArG4/OpParamSD.h"

#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4TrackStatus.hh"

namespace larg4{


  OpParamSD::OpParamSD(G4String DetectorUniqueName, std::string ModelName, int Orientation, std::vector<std::vector<double> > ModelParameters)
    : G4VSensitiveDetector(DetectorUniqueName)
  {
    // Register self with sensitive detector manager
    G4SDManager::GetSDMpointer()->AddNewDetector(this);

    if(ModelName == "OverlaidWireplanes")
      fOpa = new OverlaidWireplanesAction(ModelParameters, Orientation);

    else if(ModelName == "TransparentPlaneAction")
      fOpa = new TransparentPlaneAction();

//     else if(ModelName == "SimpleWireplane")
//       fOpa = new SimpleWireplaneAction(ModelParameters, Orientation);


    //   else if( your model here )

    else
      {
        throw cet::exception("OpParamSD")<<"Error: Optical parameterization model " << ModelName <<" not found.\n";
      }

  }


  //--------------------------------------------------------

  G4bool OpParamSD::ProcessHits(G4Step * aStep, G4TouchableHistory *)
  {

    const G4Track*           aTrack    = aStep->GetTrack();
    const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();


    G4ThreeVector mom = aParticle->GetMomentumDirection();
    G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();
    if(!fPhotonAlreadyCrossed[aTrack->GetTrackID()])
      {
	if(G4BooleanRand(fOpa->GetAttenuationFraction(mom,pos)))
	  {
	    // photon survives - let it carry on
	    fPhotonAlreadyCrossed[aTrack->GetTrackID()]=true;
	  }
	else
	  {
	    // photon is absorbed
	    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	  }
      }
    return true;
  }

  //--------------------------------------------------------

  void OpParamSD::Initialize(G4HCofThisEvent *)
  {
    fPhotonAlreadyCrossed.clear();
  }

}
