// Adapted for LArSoft by Ben Jones, MIT, Sept 09
//
// This class is based upon the default Geant4 optical physics constructor.
// However, for optical physics in LArSoft:
//
// - Cerenkov light can be enabled or disabled as required
//
// - The default geant4 reflection process is swapped for a simpler
//     custom boundary process in LArG4
//
// Parameters relating to the physics in this constructor can be set
// in the Utilities/LArParameters service.
//
//
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

#include "larsim/LArG4/OpticalPhysics.hh"
#include "Geant4/G4Alpha.hh"
#include "Geant4/G4AntiProton.hh"
#include "Geant4/G4Cerenkov.hh"
#include "Geant4/G4Deuteron.hh"
#include "Geant4/G4Electron.hh"
#include "Geant4/G4Gamma.hh"
#include "Geant4/G4GenericIon.hh"
#include "Geant4/G4He3.hh"
#include "Geant4/G4KaonMinus.hh"
#include "Geant4/G4KaonPlus.hh"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4MuonMinus.hh"
#include "Geant4/G4MuonPlus.hh"
#include "Geant4/G4OpAbsorption.hh"
#include "Geant4/G4OpBoundaryProcess.hh"
#include "Geant4/G4OpRayleigh.hh"
#include "Geant4/G4OpWLS.hh"
#include "Geant4/G4OpticalPhoton.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4PionMinus.hh"
#include "Geant4/G4PionPlus.hh"
#include "Geant4/G4Positron.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4Proton.hh"
#include "Geant4/G4Scintillation.hh"
#include "Geant4/G4Triton.hh"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "larsim/LArG4/CustomPhysicsFactory.hh"
#include "larsim/LArG4/OpBoundaryProcessSimple.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Register optical physics in custom physics list

namespace larg4 {

  CustomPhysicsFactory<OpticalPhysics> optical_factory("Optical");

  //-----------------------------------------------------------
  OpticalPhysics::OpticalPhysics(G4int ver, const G4String& name)
    : G4VPhysicsConstructor(name), verbose(ver)
  {
    G4LossTableManager::Instance();
    mf::LogInfo("OpticalPhysics") << "OBJECT BEING CONSTRUCTED IN OPTICAL PHYSICS";
  }


  //-----------------------------------------------------------
  OpticalPhysics::~OpticalPhysics()
  {}

  //-----------------------------------------------------------
  void OpticalPhysics::ConstructParticle()
  {
    MF_LOG_DEBUG("OpticalPhysics") << "PARTICLES BEING CONSTRUCTED IN OPTICAL PHYSICS";
    // optical photon
    G4OpticalPhoton::OpticalPhotonDefinition();

    // gamma
    G4Gamma::Gamma();

    // leptons
    G4Electron::Electron();
    G4Positron::Positron();
    G4MuonPlus::MuonPlus();
    G4MuonMinus::MuonMinus();

    // mesons
    G4PionPlus::PionPlusDefinition();
    G4PionMinus::PionMinusDefinition();
    G4KaonPlus::KaonPlusDefinition();
    G4KaonMinus::KaonMinusDefinition();

    // barions
    G4Proton::Proton();
    G4AntiProton::AntiProton();

    // ions
    G4Deuteron::Deuteron();
    G4Triton::Triton();
    G4He3::He3();
    G4Alpha::Alpha();
    G4GenericIon::GenericIonDefinition();
  }

  //-----------------------------------------------------------
  void OpticalPhysics::ConstructProcess()
  {

    const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
    const detinfo::DetectorProperties* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Add standard EM Processes
    MF_LOG_DEBUG("OpticalPhysics") << "PROCESSES BEING CONSTRUCTED IN OPTICAL PHYSICS";

    fTheCerenkovProcess            = new G4Cerenkov("Cerenkov");
    fTheScintillationProcess       = new G4Scintillation("Scintillation");
    fTheAbsorptionProcess          = new G4OpAbsorption();
    fTheRayleighScatteringProcess  = new G4OpRayleigh();
    if(detp->SimpleBoundary())
      fTheBoundaryProcess          = new OpBoundaryProcessSimple();
    else
      fTheBoundaryProcess_g4       = new G4OpBoundaryProcess();
    fTheWLSProcess                 = new G4OpWLS();



    fTheCerenkovProcess->SetMaxNumPhotonsPerStep(700);
    fTheCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    fTheCerenkovProcess->SetTrackSecondariesFirst(false);

    fTheScintillationProcess->SetScintillationYieldFactor(1.);
    fTheScintillationProcess->SetTrackSecondariesFirst(false);

    // Use Birks Correction in the Scintillation process

    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    fTheScintillationProcess->AddSaturation(emSaturation);


    bool CerenkovLightEnabled = larp->CerenkovLightEnabled();

    mf::LogInfo("OpticalPhysics")<<"Cerenkov light enabled : " << CerenkovLightEnabled;
    static G4ParticleTable* fParticleTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator*  aParticleIterator;
    aParticleIterator=fParticleTable->GetIterator();

    aParticleIterator->reset();
    while( (*aParticleIterator)() ){
      G4ParticleDefinition* particle = aParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (fTheCerenkovProcess->IsApplicable(*particle)&&CerenkovLightEnabled) {
	pmanager->AddProcess(fTheCerenkovProcess);
	pmanager->SetProcessOrdering(fTheCerenkovProcess,idxPostStep);
	//	mf::LogInfo("OpticalPhysics")<<"OpticalPhysics : Cerenkov applicable : " << particleName;
      }
      if (fTheScintillationProcess->IsApplicable(*particle)) {
	pmanager->AddProcess(fTheScintillationProcess);
	pmanager->SetProcessOrderingToLast(fTheScintillationProcess, idxAtRest);
	pmanager->SetProcessOrderingToLast(fTheScintillationProcess, idxPostStep);
	//	mf::LogInfo("OpticalPhysics")<<"OpticalPhysics : Scintillation applicable : " << particleName;
      }

     if (particleName == "opticalphoton") {
       mf::LogInfo("OpticalPhysics") << " AddDiscreteProcess to OpticalPhoton ";
	pmanager->AddDiscreteProcess(fTheAbsorptionProcess);
	pmanager->AddDiscreteProcess(fTheRayleighScatteringProcess);
	if(detp->SimpleBoundary())
	  pmanager->AddDiscreteProcess(fTheBoundaryProcess);
	else
	  pmanager->AddDiscreteProcess(fTheBoundaryProcess_g4);
	pmanager->AddDiscreteProcess(fTheWLSProcess);
      }
    }

  }
}
