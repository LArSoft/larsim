////////////////////////////////////////////////////////////////////////
/// \file  PhysicsList.h
/// \brief Create the physics lists to be used by Geant4.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
//
// Don't be too confused by the names.  PhysicsList.h and
// PhysicsList.cxx define what the name "LArG4::PhysicList" means.
// However, that definition is mainly in terms of
// LArG4::ModularPhysicsList, a class that inherits from
// G4VModularPhysicsList.

#include "larsim/LArG4/PhysicsList.h"

#include "Geant4/G4VModularPhysicsList.hh"
#include "Geant4/G4ParallelWorldScoringProcess.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTable.hh"

#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4ChargeExchange.hh"
#include "Geant4/G4ChargeExchangeProcess.hh"
#include "nug4/G4Base/G4PhysListFactorySingleton.hh"

PHYSLISTREG3(larg4,PhysicsList,larg4::PhysicsList)
#define G4MT_physicsVector ((G4VMPLsubInstanceManager.offset[g4vmplInstanceID]).physicsVector)

namespace larg4 {

  // Constructor: call the G4 constructor.
  ModularPhysicsList::ModularPhysicsList()
    : G4VModularPhysicsList()
  {}

  // Destructor; C++ will automatically call the G4VModulePhysicsList
  // destructor, so we need do nothing here.
  ModularPhysicsList::~ModularPhysicsList()
  {}

  // This is the method we have to modify to use the Geant4 parallel geometries.
  void ModularPhysicsList::ConstructProcess()
  {
    // We don't need to modify G4VModularPhysicsList's
    // AddTransportation method.  Just invoke it directly.
    G4VModularPhysicsList::AddTransportation();

    ////////////////////////////////////////////////////////////////////////////////////////
    // This is the "new" code that has to be added to
    // G4VModularPhysicsList::ConstructProcess() to handle the
    // parallel geometry.

    art::ServiceHandle<sim::LArG4Parameters const> lgp;
    bool DisableWireplanes = lgp->DisableWireplanes();

    G4ParallelWorldScoringProcess* LArVoxelParallelWorldScoringProcess
      = new G4ParallelWorldScoringProcess("LArVoxelReadoutScoringProcess");

    G4ParallelWorldScoringProcess* OpDetParallelWorldScoringProcess
      = new G4ParallelWorldScoringProcess("OpDetReadoutScoringProcess");

    // Note that the name below MUST match the name used in the
    // LArVoxelReadoutGeometry or OpDetReadoutGeometry constructor.
    LArVoxelParallelWorldScoringProcess->SetParallelWorld("LArVoxelReadoutGeometry");
    OpDetParallelWorldScoringProcess->SetParallelWorld("OpDetReadoutGeometry");

    // Tell all the particles in the physics list to recognize the
    // LAr voxel parallel world.  "theParticleIterator" is defined in
    // G4VUserPhysicsList.hh.  Only photons recognise the OpDet parallel
    // world.
    static G4ParticleTable* fParticleTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator*  theParticleIterator;
    theParticleIterator=fParticleTable->GetIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();



      if(particle->GetParticleType()=="opticalphoton"){
	G4ProcessManager* pmanager = particle->GetProcessManager();
	pmanager->AddProcess(OpDetParallelWorldScoringProcess);
	pmanager->SetProcessOrderingToLast(OpDetParallelWorldScoringProcess, idxAtRest);
	pmanager->SetProcessOrdering(OpDetParallelWorldScoringProcess, idxAlongStep, 1);
	pmanager->SetProcessOrderingToLast(OpDetParallelWorldScoringProcess, idxPostStep);
      }

      // Only apply voxel processing in the LAr TPC if the particles are
      // charged and they have a significant life-time.
      else if(particle->GetPDGCharge() != 0  &&  !particle->IsShortLived()){
	G4ProcessManager* pmanager = particle->GetProcessManager();
	if(!DisableWireplanes){
	  pmanager->AddProcess(LArVoxelParallelWorldScoringProcess);
	  pmanager->SetProcessOrderingToLast(LArVoxelParallelWorldScoringProcess, idxAtRest);
	  pmanager->SetProcessOrdering(LArVoxelParallelWorldScoringProcess, idxAlongStep, 1);
	  pmanager->SetProcessOrderingToLast(LArVoxelParallelWorldScoringProcess, idxPostStep);
	}

	// Do secondary biasing under control of K0Bias switch.
	G4bool genSecondaries(false);
	G4bool cE(false);
	if (((particle->GetParticleName() == G4MuonPlus::MuonPlus()->GetParticleName() ||
	      particle->GetParticleName() == G4MuonMinus::MuonMinus()->GetParticleName()
	      )
	     &&  lgp->UseCustomPhysics()
	     ) ||
	    ((particle->GetParticleName() == G4Proton::Proton()->GetParticleName() ||
	      particle->GetParticleName() == G4Neutron::Neutron()->GetParticleName() ||
	      particle->GetParticleName() == G4KaonZeroShort::KaonZeroShort()->GetParticleName() ||
	      particle->GetParticleName() == G4KaonZeroLong::KaonZeroLong()->GetParticleName() ||
	      particle->GetParticleName() == G4Lambda::Lambda()->GetParticleName()
	      )
	     &&  lgp->UseCustomPhysics()
	     )
	    ){
	  std::vector<std::string> EnabledPhysics = lgp->EnabledPhysics();
	  for(std::vector<std::string>::const_iterator it = EnabledPhysics.begin();
	      it != EnabledPhysics.end(); it++){
	    std::string PhysicsName=(*it);
	    if (!PhysicsName.compare("SynchrotronAndGN") && lgp->K0Bias()){
	      genSecondaries = true;
	    }
	    if (!PhysicsName.compare("ChargeExchange")){
	      cE = true;
	    }
	  }// end loop over enabled physics
	}// end if using custom physics

	if (genSecondaries){
	  G4int nSecondaries(lgp->K0Bias());
	  G4int fXSBias(lgp->MNXSBias());
	  G4int xSBias(lgp->MNXBias());
	  mf::LogInfo("PhysicsList: ") << "Turning on WrappedMuNuclear for "
				       << particle->GetParticleName()
				       << "s with " << nSecondaries
				       << " appropriately weighted secondaries."
				       << " XSBias is set to "
				       << xSBias
				       << " and the cross-section is increased by a factor of "
				       << fXSBias << ".";
	  G4MuonNuclearProcess* g4MNI = new G4MuonNuclearProcess();
	  MuNuclearSplittingProcessXSecBias* munuclSplitting = new MuNuclearSplittingProcessXSecBias();
	  G4bool active(true);

	  munuclSplitting->SetNSplit(nSecondaries,xSBias,fXSBias);
	  munuclSplitting->SetIsActive(active);
	  munuclSplitting->RegisterProcess(g4MNI);
	  // Middle  -1 with no biasing.
	  // The middle +1 enforces AlongStepDoIt active for MuNuclear, I claim.
	  if (xSBias&&(fXSBias>1)) pmanager->AddProcess(munuclSplitting,-1,1,1);
	  else pmanager->AddProcess(munuclSplitting,-1,-1,1);
	}// end if generating secondaries

	if (cE){
	  mf::LogInfo("PhysicsList: ") << "Turning on ChargeExchange for "
				       << particle->GetParticleName() << "s.";
	  G4ChargeExchange* model = new G4ChargeExchange();
	  G4ChargeExchangeProcess* p = new G4ChargeExchangeProcess();

	  p->RegisterMe(model);
	  pmanager->AddDiscreteProcess(p);

	}// end if doing charge exchange
      }// end if short lived particle
    }// end loop over particles

    // End of code "added" to G4VModularPhysicsList::ConstructProcess()
    ////////////////////////////////////////////////////////////////////////////////////////

    // This code is also unchanged from
    // G4VModularPhysicsList::ConstructProcess(); it means "activate
    // the physics processes and particle combinations we've specified
    // in physicsVector."  The physicsVector is built up by the
    // pre-supplied Geant4 physics list; see PhysicsList.h for the
    // name of that list.  "physicsVector" is defined in
    // G4VModularPhysicsList.hh.
    G4PhysConstVector::iterator itr;
    for (itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
      (*itr)->ConstructProcess();
    }
  }// end ConstructProcess

} // namespace LArG4
