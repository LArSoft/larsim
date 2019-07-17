//
// This(these) is (are) a variance reduction technique(s).
//
//
// Secondary biasing code adapted from EventBiasing from Jane Tinslay, SLAC. 2008.
// (Who's now at Cisco, I think.) XSBiasing code from http://geant4.cern.ch/collaboration/workshops/workshop2002/slides/docs/flei/biasing.pdf. (And NOT from http://www.lcsim.org/software/geant4/doxygen/html/classG4HadronicProcess.html, where I *believe* the point is to bias primaries and secondaries from sub-classes of Cross-Sections in an already chosen
// broader class of XSs. This doesn't help us. We need more occurrences (albeit with small
// weights) of a specific process.
//
//
// The point is to create Nsplit secondaries with each muNucl process, so
// that we may, e.g., create many, otherwise-rare, K0s which
// may charge exchange into pernicious K+s which then mimic proton dk. And,
// additionally now to make the XSection for the muNucl process itself higher,
// keeping track of appropriate weights everywhere.
//
// echurch@fnal.gov, May, 2011.
//


#include "larsim/LArG4/MuNuclearSplittingProcessXSecBias.h"

#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4KaonZeroShort.hh"
#include "Geant4/G4LorentzVector.hh"
#include "Geant4/G4ParticleChange.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4VParticleChange.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4KaonZeroLong.hh"
#include "Geant4/G4ios.hh"

#include "cetlib_except/exception.h"

#include <TMath.h>

namespace larg4 {



  G4VParticleChange* MuNuclearSplittingProcessXSecBias::PostStepDoIt(const G4Track& track, const G4Step& step)
  {

    G4VParticleChange* pChange = pRegProcess->PostStepDoIt( track, step );
    pChange->SetVerboseLevel(0);
    pChange->SetSecondaryWeightByProcess(true);

    G4double secWeight = 1.0;
    if (fNSplit>=1) secWeight = 1.0/fNSplit;
    std::vector<G4Track*> secondaries; // Secondary store

    G4double nw=1.0;
    // Let's go bias the underlying xsection.
    if(xBiasMode==1)
      {
  G4double s = step.GetStepLength();
  G4double x = pRegProcess->GetCurrentInteractionLength();
  nw =  pChange->GetParentWeight();
  // need to undo the change made to the parent track in
  // most recent step (1/wc).
  // It seems to be possible to get here w.o. first having hit the
  // AlongStep Method that initialied wc. So protect against division y 0.
  if (wc < 1.0e-10) wc = 1.0;
  nw *= 1/wc * (1. - exp(-s/x))/(1 - exp(-eFactor*s/x)) ;
  G4cout << " MNSPXSB.PostStepDoit(): original wt = " <<  pChange->GetParentWeight() << " eF = " << eFactor << G4endl;
  G4cout << " MNSPXSB.PostStepDoit(): New weight = " << nw << " wc = " << wc << " s= " << s <<" x = " << x <<  G4endl;
  ((G4ParticleChange*)pChange)->ProposeParentWeight(nw) ;
      }
    else if (xBiasMode==2)
      {
  nw = pChange->GetParentWeight()*XBiasSecondaryWeight();
  ((G4ParticleChange*)pChange)->ProposeParentWeight(nw) ;
  //  for (G4int i = 0; i < pChange->GetNumberOfSecondaries(); i++)
  //pChange->GetSecondary(i)->SetWeight(nw);
  if(G4UniformRand()<XBiasSurvivalProbability())
    { // need to add the primary back in to balance the weight
      // the position of the extra track should be moved to the next volume, but I don't know how yet.
      // need to change number of secondary first before adding an additional track
      pChange->SetNumberOfSecondaries(pChange->GetNumberOfSecondaries()+1);
      G4ThreeVector mDirection = track.GetDynamicParticle()->GetMomentumDirection();
      G4DynamicParticle* aD = new G4DynamicParticle((track.GetDynamicParticle())->GetDefinition(),
                G4LorentzVector(mDirection,
                    (track.GetDynamicParticle())->GetKineticEnergy()+step.GetTotalEnergyDeposit()));
      G4Track* secTrk = new G4Track(aD,step.GetDeltaTime(),step.GetDeltaPosition());
      pChange->AddSecondary(secTrk);
      pChange->GetSecondary(pChange->GetNumberOfSecondaries()-1)->
        SetWeight(XBiasSurvivalProbability()*track.GetWeight());
    }
  else
    {
      nw = track.GetWeight();
      pChange->ProposeParentWeight(nw);
    }
      }

    secWeight = secWeight*nw;

    // Loop over PostStepDoIt method to generate multiple secondaries.
    // The point is for each value of ii up to Nsplit we want to re-toss
    // all secondaries for the muNucl process. Each toss gives back numSec
    // secondaries, which will be a different sample of species/momenta each time.
    G4VParticleChange* particleChange = new G4VParticleChange();

    for (unsigned int ii=0; ii<(unsigned int)fNSplit; ii++) {
      particleChange = pRegProcess->PostStepDoIt(track, step);
      if (!particleChange)
        throw std::runtime_error("MuNuclearSplittingProcessXSecBias::PostStepDoIt(): no particle change");
      G4int j(0);
      G4int numSec(particleChange->GetNumberOfSecondaries());
      // Don't change weight of last secondary. It's the just-added primary in this mode.
      for (j=0; j<numSec; j++)
  {
    G4Track* newSec = new G4Track(*(particleChange->GetSecondary(j)));
    G4String pdgstr = newSec->GetParticleDefinition()->GetParticleName();
    //G4double ke = newSec->GetKineticEnergy()/GeV;
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
      }
  }
    }

    pChange->SetNumberOfSecondaries(secondaries.size());
    pChange->SetSecondaryWeightByProcess(true);
    //pChange->ProposeTrackStatus(fStopAndKill); // End of the story for this muon.

    //int numSecAdd = secondaries.size();
    std::vector<G4Track*>::iterator iter = secondaries.begin(); // Add all secondaries
    while (iter != secondaries.end()) {
      G4Track* myTrack = *iter;
      G4String pdgstr = myTrack->GetParticleDefinition()->GetParticleName();
      G4double ke = myTrack->GetKineticEnergy()/CLHEP::GeV;
      G4int pdg = myTrack->GetParticleDefinition()->GetPDGEncoding();
      if (abs(pdg)==130 ||abs(pdg)==310 ||abs(pdg)==311 || abs(pdg)==3122 || abs(pdg)==2112)
  {
    if (pdg!=2112)
      std::cout << "MuNuclSplXSB(): Adding " << pdgstr << " of Kinetic Energy and weight " << ke << " GeV and " << secWeight << "." << std::endl;

    myTrack->SetWeight(secWeight);
    pChange->AddSecondary(myTrack);
  }
      iter++;
    }
    return pChange;

  }

G4VParticleChange* MuNuclearSplittingProcessXSecBias::AlongStepDoIt(
                                                    const G4Track& track,
                                                    const G4Step& step
                )
  {

    wc = 1.0;

    if (xBiasMode==2)
      {
  /* I don't understand this. If I new it and Initialize(track) it works, but
     it goes off and grinds through the event.
     If I try to use pRegProcess->AlongStepDoIt() it crashes
     at the GetParentWeight() just below.
  */
  G4VParticleChange* pC  = new G4VParticleChange();
  pC->Initialize(track);
  //pC = pRegProcess->AlongStepDoIt( track, step);

  G4double nw = pC->GetParentWeight()*XBiasSecondaryWeight();
  pC->ProposeParentWeight(nw) ;
  for (G4int i = 0; i < pC->GetNumberOfSecondaries(); i++)
    pC->GetSecondary(i)->SetWeight(nw);
  //
  if(G4UniformRand()<XBiasSurvivalProbability())
    {
      pC->SetNumberOfSecondaries(pC->GetNumberOfSecondaries()+1);
      G4ThreeVector mDirection = track.GetDynamicParticle()->GetMomentumDirection();
      G4DynamicParticle* aD = new G4DynamicParticle((track.GetDynamicParticle())->GetDefinition(),
                G4LorentzVector(mDirection,
                    (track.GetDynamicParticle())->GetKineticEnergy()+step.GetTotalEnergyDeposit()));
      G4Track* secTrk = new G4Track(aD,track.GetGlobalTime(),track.GetPosition()); // +step.GetDeltaTime(), ,+step.GetDeltaPosition()
      secTrk->SetGoodForTrackingFlag(true);
      pC->AddSecondary(secTrk); // Add the primary right back.
      pC->GetSecondary(pC->GetNumberOfSecondaries()-1)->SetWeight(XBiasSurvivalProbability()*track.GetWeight());
      //pC->ProposeTrackStatus(fStopAndKill); // 8-June-2011-21:36 MDT, from Finland, I (EDC) added this. Seems right. Makes it go, anyway. Almost too fast.
    }

  return pC;
      }
    else if (xBiasMode==1)
      {
  fParticleChange.Initialize(track) ;
  G4double s = step.GetStepLength();
  G4double x = pRegProcess->GetCurrentInteractionLength();
  //fParticleChange = &(pRegProcess->AlongStepDoIt( track, step));
  wc =  exp(-(1. - eFactor)*s/x);
  G4double w = wc * fParticleChange.GetParentWeight();
  fParticleChange.SetParentWeightByProcess(false);
  fParticleChange.ProposeParentWeight(w) ;
  if (w>100)
    {
      //    G4cout << " MNSPXSB.AlongStepDoit(): GPW is " << w/wc << " wc = " << wc << ", CIL = "<< x << G4endl;
    }

  return &fParticleChange;
      }

    if (xBiasMode!=1 && xBiasMode!=2)
      {
  throw cet::exception("Incorrectly set Bias Mode. ")
    << "Set XBiasMode to 0 or 1. " << xBiasMode << " not allowed!\n";
      }
  fParticleChange.Initialize(track) ;
  return &fParticleChange;
  }


 G4double MuNuclearSplittingProcessXSecBias::XBiasSurvivalProbability()
  {
   G4double result = 0;
   G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
   G4double biasedProbability = 1.-std::exp(-nLTraversed);
   G4double realProbability = 1-std::exp(-nLTraversed/eFactor);
   result = (biasedProbability-realProbability)/biasedProbability;
   return result;
  }

  G4double MuNuclearSplittingProcessXSecBias::XBiasSecondaryWeight()
  {
    G4double result = 0;
    G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
    result = 1./eFactor*std::exp(-nLTraversed/eFactor*(1-1./eFactor));
    return result;
  }


  G4double MuNuclearSplittingProcessXSecBias::GetTotalNumberOfInteractionLengthTraversed()
  {
    return theInitialNumberOfInteractionLength
      -G4VProcess::theNumberOfInteractionLengthLeft;
  }


  G4double MuNuclearSplittingProcessXSecBias::PostStepGetPhysicalInteractionLength( const G4Track& track,
                                          G4double   previousStepSize,
                                          G4ForceCondition* condition )
  {
    // I believe this is the money line for the whole thing. By shrinking the
    // interaction length in MuNuclear's process, it's more likely that
    // this process will be chosen in the competition among processes.
    if(xBiasMode)
      {
  G4double psgPIL = (1./eFactor) * pRegProcess->PostStepGetPhysicalInteractionLength(
                    track,
             previousStepSize,
                                     condition );
  return psgPIL;
      }

    return pRegProcess->PostStepGetPhysicalInteractionLength(
                                     track,
             previousStepSize,
                                     condition );

}


G4double MuNuclearSplittingProcessXSecBias::AlongStepGetPhysicalInteractionLength( const G4Track& track,
                                     G4double previousStepSize,
                                     G4double currentMinimumStep,
                                     G4double& proposedSafety,
                                     G4GPILSelection* selection )
{
  if (xBiasMode==2)
    {
      // I am pretty sure this mode does not work because MuNuclear doesn't interact
      // AlongStep. Below often (always?) returns a negative number, which I have to
      // swizzle.
      G4double intLen;
      intLen = (1./eFactor) *
      pRegProcess->AlongStepGetPhysicalInteractionLength( track,previousStepSize,
                currentMinimumStep,proposedSafety,selection);
      if (intLen<0) return currentMinimumStep;
      else return intLen;
    }
  else
    return DBL_MAX ;
}


} // end namespace
