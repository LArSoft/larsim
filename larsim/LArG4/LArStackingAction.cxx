//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "larsim/LArG4/LArStackingAction.h"
#include "larcore/Geometry/Geometry.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Geant4/G4MuonMinus.hh"
#include "Geant4/G4MuonPlus.hh"
#include "Geant4/G4StackManager.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4TrackStatus.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4VProcess.hh"

// ROOT includes
#include "RtypesCore.h"
#include "TString.h"
#include "TVector3.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"

LArStackingAction::LArStackingAction(G4int dum)
 : fstage(0)
 , freqMuon(2)
 , freqIsoMuon(0)
 , freqIso(10)
 , fangRoI(30.*CLHEP::deg)
{
  //theMessenger = new LArStackingActionMessenger(this);
  fStack = dum;
    // Positive values effect action in this routine. Negative values
    // effect action in G4BadIdeaAction.

}

LArStackingAction::~LArStackingAction()
{ //delete theMessenger;
}

G4ClassificationOfNewTrack
LArStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fWaiting;
  art::ServiceHandle<geo::Geometry const> geom;
  TString volName(InsideTPC(aTrack));
  Double_t buffer = 500; // Keep muNucl neutrals within 5m (for now) of larTPC.

  // These 3 for now for investigation in gdb.
  //int pdg = aTrack->GetDefinition()->GetPDGEncoding();
  int ppdg = aTrack->GetParentID();
  TString process("NA");
  if (ppdg) process = (TString)aTrack->GetCreatorProcess()->GetProcessName();


  switch(fstage){
  case 0: // Fstage 0 : Primary muons only
    if (aTrack->GetParentID()==0)
    {
      G4ParticleDefinition *particleType = aTrack->GetDefinition();
      if( ((particleType==G4MuonPlus::MuonPlusDefinition())
	  || (particleType==G4MuonMinus::MuonMinusDefinition())
	   )
	  && !volName.Contains("unknown")
	  ){
	classification = fUrgent;
      }
    }
    if (volName.Contains("unknown"))   classification = fKill;
    break;

  case 1: // Stage 1 : K0,Lambda,n's made urgent here.
          //           Suspended tracks will be sent to the waiting stack
    if(aTrack->GetTrackStatus()==fSuspend) { break; }

    if ((aTrack->GetDefinition()->GetPDGEncoding()==2112 || aTrack->GetDefinition()->GetPDGEncoding()==130 || aTrack->GetDefinition()->GetPDGEncoding()==310 || aTrack->GetDefinition()->GetPDGEncoding()==311 || aTrack->GetDefinition()->GetPDGEncoding()==3122 ) && (aTrack->GetParentID()==1) && !volName.Contains("unknown"))
      {

	const G4ThreeVector tr4Pos = aTrack->GetPosition();
	// G4 returns positions in mm, have to convert to cm for LArSoft coordinate systems
	const TVector3 trPos(tr4Pos.x()/CLHEP::cm,tr4Pos.y()/CLHEP::cm,tr4Pos.z()/CLHEP::cm);
	//double locNeut = trPos.Mag();
	classification = fUrgent;
	// std::cout << "LArStackingAction: DetHalfWidth, Height, FullLength: " << geom->DetHalfWidth() << ", " << geom->DetHalfHeight() << ", " << geom->DetLength() << std::endl;

	if (
	    trPos.X() < (geom->DetHalfWidth()*2.0 + buffer) && trPos.X() > (-buffer) &&
	    trPos.Y() < (geom->DetHalfHeight()*2.0 + buffer) && trPos.Y() > (-geom->DetHalfHeight()*2.0 - buffer) &&
	    trPos.Z() < (geom->DetLength() + buffer) && trPos.Z() > (-buffer)
	    )

	  { classification = fUrgent; break; }
	// These tracks need to be "scored" cuz every now and then they
	// might get to the LAr.
	else
	  { classification = fKill; break; }


      }

	//    if(aTrack->GetDefinition()->GetPDGCharge()==0.) { break; }
    if (volName.Contains("unknown"))   classification = fKill;
    break;

  default:
    // Track all other Primaries. Accept all secondaries in TPC.
    // Kill muon ionization electrons outside TPC
    // ignore primaries since they have no creator process

    if(aTrack->GetParentID() == 0 && !volName.Contains("unknown")){
      classification = fUrgent;
      break;
    }

    if(volName.Contains(geom->GetLArTPCVolumeName()) && aTrack->GetParentID()!=0)
      {
	classification = fUrgent;
	if (fStack & 0x4 &&
	    aTrack->GetCreatorProcess()->GetProcessName().contains("muIoni")
	    )
	  {
	    classification = fKill;
	  }
	break;
      }
    else if (volName.Contains("unknown") ){
      classification = fKill;
      break;
    }
    // Leave this here, even though I claim we've Killed these in stage 2.
    if(aTrack->GetDefinition()->GetPDGEncoding()==11
       && aTrack->GetCreatorProcess()->GetProcessName().contains("muIoni") )
      {
	classification = fKill;
	break;
      }
    // For now, kill every other thing, no matter where it is.
    classification = fKill;

  } // end switch

  return classification;
}


std::string LArStackingAction::InsideTPC(const G4Track * aTrack)
{

  art::ServiceHandle<geo::Geometry const> geom;
  const G4ThreeVector tr4Pos = aTrack->GetPosition();

  // G4 returns positions in mm, have to convert to cm for LArSoft coordinate systems
  const TVector3 trPos(tr4Pos.x()/CLHEP::cm,tr4Pos.y()/CLHEP::cm,tr4Pos.z()/CLHEP::cm);

  const std::string volName(geom->VolumeName(trPos));

  return volName;
}


void LArStackingAction::NewStage()
{

  // Here when Urgent stack is empty. Waiting stack about to be made Urgent,
  // upon saying ReClassify().
  fstage++;

    // I yanked the ExN04's use here of stackManager->clear(), which clears stack
    // and prepares to end the event. Think I may wanna do something like this if
    // muon has been tracked, its doca is large, there are no hit voxels in the TPC,
    // and I'm in further need of optimization.

  if(fstage==1){
  // Stage 0->1 : check if at least "reqMuon" hits on muon chamber
  //              otherwise abort current event

    stackManager->ReClassify();
    return;
  }

  else if(fstage==2){
  // Stage 1->2 : check the isolation of muon tracks
  //              at least "reqIsoMuon" isolated muons
  //              otherwise abort current event.
  //              Isolation requires "reqIso" or less hits
  //              (including own hits) in the RoI region
  //              in the tracker layers.
    stackManager->ReClassify();
    return;
  }

  else{
  // Other stage change : just re-classify
    stackManager->ReClassify();
  }
}

void LArStackingAction::PrepareNewEvent()
{
  fstage = 0;
  //trkHits = 0;
  //muonHits = 0;
}
