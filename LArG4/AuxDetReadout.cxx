////////////////////////////////////////////////////////////////////////
/// \file   AuxDetReadout.cxx
/// \brief  A Geant4 sensitive detector that accumulates information.
/// \author miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "LArG4/AuxDetReadout.h"
#include "LArG4/ParticleListAction.h"
#include "Utilities/DetectorProperties.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4HCofThisEvent.hh"
#include "Geant4/G4TouchableHistory.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/Randomize.hh"

#include <iostream>
#include <ctime>

namespace larg4 {
  
  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by ?.
  AuxDetReadout::AuxDetReadout(std::string const& name)
  : G4VSensitiveDetector(name)
  , fAuxDet(atoi( name.substr(15,name.size()).c_str() ))
  , fAuxDetSimChannel(sim::AuxDetSimChannel(fAuxDet))
  {
  }
  
  //---------------------------------------------------------------------------------------
  AuxDetReadout::~AuxDetReadout() {}
  
  //---------------------------------------------------------------------------------------
  // Called at the start of each event.
  void AuxDetReadout::Initialize(G4HCofThisEvent*)
  {
    //fAuxDetSimChannel.clear();
  }
  //---------------------------------------------------------------------------------------
  // Called at the end of each event.
  void AuxDetReadout::EndOfEvent(G4HCofThisEvent*)
  {
  }
  //---------------------------------------------------------------------------------------
  void AuxDetReadout::clear()
  {
  }

  //---------------------------------------------------------------------------------------
  // Called for each step. Create a vector of AuxDetSimTrack objects. One for each new TrackID.
  // update the exit position points if this TrackID is already on the vector.
  G4bool AuxDetReadout::ProcessHits( G4Step* step, G4TouchableHistory* )
  {
		// collect the info for this step

    const int trackID = ParticleListAction::GetCurrentTrackID();

    G4double energyDeposited =   step->GetTotalEnergyDeposit()/GeV;

		G4ThreeVector startG4(step->GetPreStepPoint()->GetPosition() );
		double startWorld[3]={startG4.getX()/cm,startG4.getY()/cm,startG4.getZ()/cm};
		
		double startTime = step->GetPreStepPoint()->GetGlobalTime()/ns;

    G4ThreeVector stopG4( step->GetPostStepPoint()->GetPosition());
		double stopWorld[3]={stopG4.getX()/cm,stopG4.getY()/cm,stopG4.getZ()/cm};
    
    G4ThreeVector stopG4Momentum( step->GetPostStepPoint()->GetMomentum());
    double stopWorldMomVector[3]={stopG4Momentum.getX()/cm,stopG4Momentum.getY()/cm,stopG4Momentum.getZ()/cm};

		double stopTime = step->GetPostStepPoint()->GetGlobalTime()/ns;
		
    fAuxDetSimChannel.AddParticleStep( trackID,
                                       energyDeposited,
                                       startWorld[0],
                                       startWorld[1],
                                       startWorld[2],
                                       startTime,
                                       stopWorld[0],
                                       stopWorld[1],
                                       stopWorld[2],
                                       stopTime,
                                       stopWorldMomVector[0],
                                       stopWorldMomVector[1],
                                       stopWorldMomVector[2]
                                     );
    
    return true;
  }
  
  //---------------------------------------------------------------------------------------
  // Never used but still have to be defined for G4
  void AuxDetReadout::DrawAll()  {}
  void AuxDetReadout::PrintAll() {}
  
} // namespace larg4
