////////////////////////////////////////////////////////////////////////
/// \file   AuxDetReadout.cxx
/// \brief  A Geant4 sensitive detector that accumulates information.
/// \author miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larsim/LArG4/AuxDetReadout.h"
#include "larsim/LArG4/ParticleListAction.h"

#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"

#include <utility> // std::move()
#include <algorithm> // std::find()

namespace larg4 {

  AuxDetReadout::AuxDetReadout(std::string const& name,
			       unsigned int       adNum,
			       unsigned int       svNum)
  : G4VSensitiveDetector(name)
  , fAuxDet(adNum)
  , fAuxDetSensitive(svNum)
  {
  }

  //----------------------------------------------------------------------
  AuxDetReadout::~AuxDetReadout() {}

  //---------------------------------------------------------------------------------------
  // Called at the start of each event.
  void AuxDetReadout::Initialize(G4HCofThisEvent*)
  {}
  //---------------------------------------------------------------------------------------
  // Called at the end of each event.
  void AuxDetReadout::EndOfEvent(G4HCofThisEvent*)
  {
    fAuxDetSimChannel = sim::AuxDetSimChannel(fAuxDet, std::move(fAuxDetIDEs), fAuxDetSensitive);
  }
  //---------------------------------------------------------------------------------------
  void AuxDetReadout::clear()
  {
    fAuxDetIDEs.clear();
  }

  //---------------------------------------------------------------------------------------
  // Called for each step. Create a vector of AuxDetSimTrack objects. One for each new TrackID.
  // update the exit position points if this TrackID is already on the vector.
  G4bool AuxDetReadout::ProcessHits( G4Step* step, G4TouchableHistory* )
  {
		// collect the info for this step

    const int trackID = ParticleListAction::GetCurrentTrackID();

    G4double energyDeposited = step->GetTotalEnergyDeposit()/CLHEP::GeV;

    G4ThreeVector startG4(step->GetPreStepPoint()->GetPosition() );
    double startWorld[3] = {startG4.getX()/CLHEP::cm,
                            startG4.getY()/CLHEP::cm,
                            startG4.getZ()/CLHEP::cm};

    double startTime = step->GetPreStepPoint()->GetGlobalTime()/CLHEP::ns;

    G4ThreeVector stopG4( step->GetPostStepPoint()->GetPosition());
    double stopWorld[3] = {stopG4.getX()/CLHEP::cm,
                           stopG4.getY()/CLHEP::cm,
                           stopG4.getZ()/CLHEP::cm};

    G4ThreeVector stopG4Momentum( step->GetPostStepPoint()->GetMomentum());
    double stopWorldMomVector[3] = {stopG4Momentum.getX()/CLHEP::GeV,
                                    stopG4Momentum.getY()/CLHEP::GeV,
                                    stopG4Momentum.getZ()/CLHEP::GeV};

    double stopTime = step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;

    this->AddParticleStep( trackID,
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

  // Moved here from AuxDetSimChannel.cxx
  void AuxDetReadout::AddParticleStep(
  			int	inputTrackID,
			float	inputEnergyDeposited,
			float	inputEntryX,
			float	inputEntryY,
			float	inputEntryZ,
			float	inputEntryT,
			float	inputExitX,
			float	inputExitY,
			float	inputExitZ,
			float	inputExitT,
			float	inputExitMomentumX,
			float	inputExitMomentumY,
			float	inputExitMomentumZ){

    sim::AuxDetIDE auxDetIDE;
    auxDetIDE.trackID		= inputTrackID;
    auxDetIDE.energyDeposited	= inputEnergyDeposited;
    auxDetIDE.entryX		= inputEntryX;
    auxDetIDE.entryY		= inputEntryY;
    auxDetIDE.entryZ		= inputEntryZ;
    auxDetIDE.entryT		= inputEntryT;
    auxDetIDE.exitX		= inputExitX;
    auxDetIDE.exitY		= inputExitY;
    auxDetIDE.exitZ		= inputExitZ;
    auxDetIDE.exitT		= inputExitT;
    auxDetIDE.exitMomentumX	= inputExitMomentumX;
    auxDetIDE.exitMomentumY	= inputExitMomentumY;
    auxDetIDE.exitMomentumZ	= inputExitMomentumZ;

    std::vector<sim::AuxDetIDE>::iterator IDEitr
      = std::find(fAuxDetIDEs.begin(), fAuxDetIDEs.end(), auxDetIDE);

    if(IDEitr != fAuxDetIDEs.end()){ //If trackID is already in the map, update it

      IDEitr->energyDeposited += inputEnergyDeposited;
      IDEitr->exitX            = inputExitX;
      IDEitr->exitY            = inputExitY;
      IDEitr->exitZ            = inputExitZ;
      IDEitr->exitT            = inputExitT;
      IDEitr->exitMomentumX    = inputExitMomentumX;
      IDEitr->exitMomentumY    = inputExitMomentumY;
      IDEitr->exitMomentumZ    = inputExitMomentumZ;
    }
    else{  //if trackID is not in the set yet, add it
      fAuxDetIDEs.push_back(std::move(auxDetIDE));
    }//else
  }//AddParticleStep

  //---------------------------------------------------------------------------------------
  // Never used but still have to be defined for G4
  void AuxDetReadout::DrawAll()  {}
  void AuxDetReadout::PrintAll() {}

} // namespace larg4
