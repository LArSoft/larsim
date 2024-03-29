////////////////////////////////////////////////////////////////////////
/// \file  G4BadIdeaAction.cxx
/// \brief Use Geant4's user "hooks" to maintain a list of particles generated by Geant4.
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larsim/LegacyLArG4/G4BadIdeaAction.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Track.hh"

namespace larg4 {

  //----------------------------------------------------------------------------
  // Constructor.
  G4BadIdeaAction::G4BadIdeaAction(int trkOption) : fNoIncomingMuons(trkOption)
  {
    // trkOption comes from LArG4's fSmartStacking
    // Negative values effect action in this routine. Positive values
    // effect action in LArStackingAction.
    mf::LogWarning("G4BadIdeaAction") << "instantiating the G4BadIdeaAction \n"
                                      << "This UserAction is only to be used with "
                                      << "Geant4 v4.9.4.p02 to solve a stepping bug.\n"
                                      << "If you are using a different version of G4, "
                                      << "remove this UserAction from your list in LArG4.cxx";
  }

  //----------------------------------------------------------------------------
  // Destructor.
  G4BadIdeaAction::~G4BadIdeaAction() {}

  //----------------------------------------------------------------------------
  // With every step, add to the particle's trajectory.
  void G4BadIdeaAction::SteppingAction(const G4Step* step)
  {
    ////////////////////////////////////////////////////////////////////
    ///\todo Do not copy the code below. Contact Brian Rebel, Eric Church,
    ///\todo Bill Seligman and Andrzej Szelc for reasons why not to!
    ////////////////////////////////////////////////////////////////////////

    // If the step size is such that the particle appears to be "stuck"
    // in its trajectory, give it a kick.
    const double epsilon = 5000. * std::numeric_limits<double>::epsilon();
    const double stepSize = step->GetStepLength();
    G4Track* nonConstTrack = const_cast<G4Track*>(step->GetTrack());

    if (step->GetTrack()->GetCurrentStepNumber() > 5e4 && stepSize < epsilon &&
        (step->GetTrack()->GetCurrentStepNumber() % 1000) == 1) {

      // Cast away the const-ness of the pointer to G4Step.
      // This is dangerous. Don't do this at home. We're
      // only doing this because we're desperate.
      // The need to do this is the result of a bug in Geant4 v4.9.4.p02
      // We should no longer call this code when we move beyond that version
      const double kick = 0.001;

      G4ThreeVector aValue = nonConstTrack->GetPosition();

      mf::LogWarning("G4BadIdeaAction")
        << "##### In endless loop. Kicking particle by "
        << " (+0.001,+0.001,+0.001) --- "
        << " PDG and encoding " << step->GetTrack()->GetDynamicParticle()->GetPDGcode() << " "
        << step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()
        << " current step number: " << step->GetTrack()->GetCurrentStepNumber()
        << " stepsize: " << stepSize << " x,y,z  " << aValue.x() << " " << aValue.y() << " "
        << aValue.z();

      G4ThreeVector translate(kick, kick, kick);
      aValue += translate;

      nonConstTrack->SetPosition(aValue);
    }

    if (fNoIncomingMuons < 0) {
      // This is for overlays of, say, rock muons, which we have
      G4StepPoint* thePrePoint = step->GetPreStepPoint();
      G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
      G4String thePrePVname = thePrePV->GetName();
      G4StepPoint* thePostPoint = step->GetPostStepPoint();
      G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
      G4String thePostPVname("null");
      if (thePostPV) thePostPVname = thePostPV->GetName();

      if (abs(step->GetTrack()->GetDynamicParticle()->GetPDGcode()) == 13 &&
          thePostPVname.contains("volTPCActive") && !thePrePVname.contains("volTPCActive"))
        ((G4Track*)nonConstTrack)->SetTrackStatus(fStopAndKill);
    }

    return;
  }

} // namespace LArG4
