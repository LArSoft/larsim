////////////////////////////////////////////////////////////////////////
/// \file VisualizationAction.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
/// VisualizationAction.cc
/// 19-Mar-2002 Bill Seligman
///
/// Use UserAction to implement the standard visualization control for
/// a typical Geant4 job.  Everything in this class comes from the
/// Geant4 examples; the only difference is that it's put into an
/// UserAction class.
///
/// 25-Feb-2009 WGS: Revised for FMWK/LArSoft

#include "larsim/LArG4/VisualizationAction.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "Geant4/G4Run.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4Event.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4VVisManager.hh"
#include "Geant4/G4TrajectoryContainer.hh"
#include "Geant4/G4Trajectory.hh"
#include "Geant4/G4ios.hh"

#include <cmath>
#include <cstring>

namespace larg4 {

  VisualizationAction::VisualizationAction()
  {
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    m_energyCutoff = lgp->VisualizationEnergyCut()*GeV;
    m_drawNeutrals = lgp->DrawNeutrals();

  }

  VisualizationAction::~VisualizationAction()
  {
  }

  void VisualizationAction::BeginOfRunAction( const G4Run* )
  {
    // If graphics is ON
    if(G4VVisManager::GetConcreteInstance()){
      G4cout << "Geant4 visualisation is ON" << G4endl;
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/tracking/storeTrajectory 1");
      UI->ApplyCommand("/vis/scene/notifyhandlers");
    }
  }

  void VisualizationAction::EndOfRunAction( const G4Run* )
  {
    // If graphics is ON
    if (G4VVisManager::GetConcreteInstance())
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  void VisualizationAction::BeginOfEventAction( const G4Event* )
  {}

  void VisualizationAction::EndOfEventAction(const G4Event* a_event)
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

    if (pVVisManager)
      {
	G4TrajectoryContainer* trajectoryContainer = a_event->GetTrajectoryContainer();
	G4int n_trajectories = 0;
	if (trajectoryContainer)
	  n_trajectories = trajectoryContainer->entries();

	for (G4int i=0; i < n_trajectories; i++)
	  {
	    G4VTrajectory* trajectory = (*trajectoryContainer)[i];

	    // Only draw "interesting" trajectories, according to cuts
	    // and limits set by the user.

	    G4bool drawThisTrack = false;
	    if ( trajectory->GetInitialMomentum().mag() > m_energyCutoff )
	      drawThisTrack = true;

	    if ( ! m_drawNeutrals )
	      {
		G4double charge = trajectory->GetCharge();
		// electron-, muon-, and tau-neutrino codes:
		if ( charge == 0 )
		  drawThisTrack = false;
	      }

	    if ( drawThisTrack )
	      trajectory->DrawTrajectory(50);
	  }
      }
  }

  


} // namespace LArG4
