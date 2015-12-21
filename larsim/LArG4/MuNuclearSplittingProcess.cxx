//
// Code adapted from EventBiasing from Jane Tinslay, SLAC. 2008. Now at Cisco,
// I think.
//
// This is a variance reduction technique.
//
// The point is to create Nsplit secondaries with each muNucl process, so
// that we may, e.g., create many, otherwise-rare, K0s which
// may charge exchange into pernicious K+s which then mimic proton dk.
//
// echurch@fnal.gov, May, 2011.
//


#include "larsim/LArG4/MuNuclearSplittingProcess.h"

#include <stdexcept> // std::runtime_error

#include "Geant4/G4WrapperProcess.hh"
#include "Geant4/G4VParticleChange.hh"

#include "Geant4/G4SDManager.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4Event.hh"
#include "Geant4/G4HCofThisEvent.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4TrackStatus.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4KaonZeroLong.hh"
#include "Geant4/G4ios.hh"

#include <TMath.h>

namespace larg4 {

G4VParticleChange* MuNuclearSplittingProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{

  G4VParticleChange* particleChange = new G4VParticleChange();
  G4double weight = track.GetWeight()/fNSplit;
  std::vector<G4Track*> secondaries; // Secondary store

  // Loop over PostStepDoIt method to generate multiple secondaries.
  // The point is for each value of ii up to Nsplit we want to re-toss
  // all secondaries for the muNucl process. Each toss gives back numSec
  // secondaries, which will be a different sample of species/momenta each time.
  for (unsigned int ii=0; ii<(unsigned int)fNSplit; ii++) {
    particleChange = pRegProcess->PostStepDoIt(track, step);
    if (!particleChange)
      throw std::runtime_error("MuNuclearSplittingProcess::PostStepDoIt(): no particle change");
    G4int j(0);
    G4int numSec(particleChange->GetNumberOfSecondaries());
    for (j=0; j<numSec; j++)
      {
  G4Track* newSec = new G4Track(*(particleChange->GetSecondary(j)));
  G4String pdgstr = newSec->GetParticleDefinition()->GetParticleName();
  G4double ke = newSec->GetKineticEnergy()/GeV;
  G4int pdg = newSec->GetParticleDefinition()->GetPDGEncoding();
  if (abs(pdg)==310 ||abs(pdg)==311 || abs(pdg)==3122 || abs(pdg)==2112)
    {
      //      std::cout << "MuNuclSplProc: Creating " << pdgstr << " of Kinetic Energy " << ke << " GeV. numSec is " << j << std::endl;

      if (pdg==G4KaonZeroShort::KaonZeroShort()->GetPDGEncoding()&&G4UniformRand()<0.50)
        {
    pdg = G4KaonZeroLong::KaonZeroLong()->GetPDGEncoding();
    pdgstr = G4KaonZeroLong::KaonZeroLong()->GetParticleName();
    G4LorentzVector pK0L(newSec->GetMomentum(),TMath::Sqrt(TMath::Power(G4KaonZeroLong::KaonZeroLong()->GetPDGMass(),2)+TMath::Power(newSec->GetMomentum().mag(),2)));
    G4DynamicParticle *newK0L = new G4DynamicParticle(G4KaonZeroLong::KaonZeroLong(),pK0L);

    G4Track* newSecK0L = new G4Track(newK0L,track.GetGlobalTime(),track.GetPosition());
    secondaries.push_back(newSecK0L);
        }
      else
        {
    secondaries.push_back(newSec);
        }

      if (abs(pdg)==130 ||abs(pdg)==310 ||abs(pdg)==311 || abs(pdg)==3122)
        {
    // WrappedmuNuclear always produces K0s if it produces a K0 at all.
    // Let's make half of these K0Ls.
    std::cout << "MuNuclSplProc: Creating " << pdgstr << " of Kinetic Energy " << ke << " GeV. numSec is " << j << std::endl;
        }

    }
      }
  }

  particleChange->SetNumberOfSecondaries(secondaries.size());
  particleChange->SetSecondaryWeightByProcess(true);
  //int numSecAdd = secondaries.size();
  std::vector<G4Track*>::iterator iter = secondaries.begin(); // Add all secondaries
  while (iter != secondaries.end()) {
    G4Track* myTrack = *iter;
    G4String pdgstr = myTrack->GetParticleDefinition()->GetParticleName();
    //G4double ke = myTrack->GetKineticEnergy()/GeV;
    G4int pdg = myTrack->GetParticleDefinition()->GetPDGEncoding();
    if (abs(pdg)==130 ||abs(pdg)==310 ||abs(pdg)==311 || abs(pdg)==3122 || abs(pdg)==2112)
      //      std::cout << "MuNuclSplProc: Adding " << pdgstr << " of Kinetic Energy " << ke << " GeV." << std::endl;
    // Will ask for PdgCode and only Add K0s.
    myTrack->SetWeight(weight);
    particleChange->AddSecondary(myTrack);
    iter++;
  }
  return particleChange;
}

} // end namespace
