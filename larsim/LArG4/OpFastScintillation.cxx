// Class adapted for LArSoft by Ben Jones, MIT 10/10/12
//
// This class is a physics process based on the standard Geant4
// scintillation process.
//
// It has been stripped down and adapted to form the backbone of 
// the LArG4 fast optical simulation.  Photons, instead of being
// produced and added to the geant4 particle stack, are logged
// and used to predict the visibility of this step to each PMT in
// the detector.
//
// The photonvisibilityservice looks up the visibility of the relevant
// xyz point, and if a photon is detected at a given PMT, one OnePhoton 
// object is logged in the OpDetPhotonTable
//
// At the end of the event, the OpDetPhotonTable is read out
// by LArG4, and detected photons are stored in the event.
//
// This process can be used alongside the standard Cerenkov process,
// which does step geant4 opticalphotons.  Both the fast scintillation
// table and the geant4 sensitive detectors are read out by LArG4 to
// produce a combined SimPhoton collection.
//
// Added disclaimer : This code is gross.  Thats basically because
//  it adheres to the original, gross Geant4 implementation.
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
//
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        OpFastScintillation.cc
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//              of energy deposited by particle type
//              Thanks to Zach Hartwig (Department of Nuclear
//              Science and Engineeering - MIT)
//              2010-09-22 by Peter Gumplinger
//              > scintillation rise time included, thanks to
//              > Martin Goettlich/DESY
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2004-08-05 by Peter Gumplinger
//              > changed StronglyForced back to Forced in GetMeanLifeTime
//              2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumberOfPhotons
//              2002-11-07 by Peter Gumplinger
//              > now allow for fast and slow scintillation component
//              2002-11-05 by Peter Gumplinger
//              > now use scintillation constants from G4Material
//              2002-05-09 by Peter Gumplinger
//              > use only the PostStepPoint location for the origin of
//                scintillation photons when energy is lost to the medium
//                by a neutral particle
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              2001-09-17, migration of Materials to pure STL (mma)
//              2003-06-03, V.Ivanchenko fix compilation warnings
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

//#include "Geant4/g4ios.hh"
#include "TLorentzVector.h"

#include "Geant4/globals.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4EmProcessSubType.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4Poisson.hh"
#include "Geant4/G4VPhysicalVolume.hh"


#include "larsim/LArG4/ParticleListAction.h"
#include "larsim/LArG4/IonizationAndScintillation.h"
#include "larsim/LArG4/OpFastScintillation.hh"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/LArG4/OpDetPhotonTable.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// support libraries
#include "cetlib_except/exception.h"

#include "TRandom3.h"
#include "TMath.h"
#include "boost/algorithm/string.hpp"

namespace larg4{

/////////////////////////
// Class Implementation
/////////////////////////
  
  //////////////
  // Operators
  //////////////
  
  // OpFastScintillation::operator=(const OpFastScintillation &right)
  // {
  // }
  
  /////////////////
  // Constructors
  /////////////////
  
  OpFastScintillation::OpFastScintillation(const G4String& processName, G4ProcessType type)      
  : G4VRestDiscreteProcess(processName, type)
  , bPropagate(!(art::ServiceHandle<sim::LArG4Parameters>()->NoPhotonPropagation()))
  {
        G4cout<<"BEA1 timing distribution"<<G4endl;
        SetProcessSubType(25);

        fTrackSecondariesFirst = false;
        fFiniteRiseTime = false;


	YieldFactor=1.0;
        ExcitationRatio = 1.0;
	
	const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
	
        scintillationByParticleType = larp->ScintByParticleType();

        theFastIntegralTable = NULL;
        theSlowIntegralTable = NULL;

        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        BuildThePhysicsTable();

        emSaturation = NULL;

        if (bPropagate) {
          art::ServiceHandle<phot::PhotonVisibilityService> pvs;
          if(pvs->IncludePropTime()) {
            pvs->SetDirectLightPropFunctions(functions_vuv, fd_break, fd_max, ftf1_sampling_factor);
            pvs->SetReflectedCOLightPropFunctions(functions_vis, ft0_max, ft0_break_point);
          }
        }
}

  OpFastScintillation::OpFastScintillation(const OpFastScintillation& rhs)
    :  G4VRestDiscreteProcess(rhs.GetProcessName(), rhs.GetProcessType())
  {
    theSlowIntegralTable        = rhs.GetSlowIntegralTable();
    theFastIntegralTable        = rhs.GetFastIntegralTable();
    
    fTrackSecondariesFirst      = rhs.GetTrackSecondariesFirst();
    fFiniteRiseTime             = rhs.GetFiniteRiseTime();
    YieldFactor                 = rhs.GetScintillationYieldFactor();
    ExcitationRatio             = rhs.GetScintillationExcitationRatio();
    scintillationByParticleType = rhs.GetScintillationByParticleType();
    emSaturation                = rhs.GetSaturation();

    BuildThePhysicsTable();
  }


        ////////////////
        // Destructors
        ////////////////

OpFastScintillation::~OpFastScintillation()
{
	if (theFastIntegralTable != NULL) {
           theFastIntegralTable->clearAndDestroy();
           delete theFastIntegralTable;
        }
        if (theSlowIntegralTable != NULL) {
           theSlowIntegralTable->clearAndDestroy();
           delete theSlowIntegralTable;
        }

}

        ////////////
        // Methods
        ////////////

// AtRestDoIt
// ----------
//
G4VParticleChange*
OpFastScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()

{
        return OpFastScintillation::PostStepDoIt(aTrack, aStep);
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
OpFastScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is 
// generated according to the scintillation yield formula, distributed 
// evenly along the track segment and uniformly into 4pi.

{
        aParticleChange.Initialize(aTrack);

	// Check that we are in a material with a properties table, if not
	// just return
        const G4Material* aMaterial = aTrack.GetMaterial();
        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable)
             return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

        G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
      
        G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4ThreeVector p0 = aStep.GetDeltaPosition().unit();


	///////////////////////////////////////////////////////////////////////////////////
	//   This is the old G4 way - but we do things differently - Ben J, Oct Nov 2012.
	///////////////////////////////////////////////////////////////////////////////////
	//
	//     if (MeanNumberOfPhotons > 10.)
	//      {
	//        G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
	//        NumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);
	//      }
	//     else
	//      {
	//        NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
	//      }
	//
	//
	//
	//        if (NumPhotons <= 0)
	//        {
	//  // return unchanged particle and no secondaries 
	//
	//           aParticleChange.SetNumberOfSecondaries(0);
	//
	//           return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
	//        }
	//
	//
	//       aParticleChange.SetNumberOfSecondaries(NumPhotons);
	//
	//
	//        if (fTrackSecondariesFirst) {
	//           if (aTrack.GetTrackStatus() == fAlive )
	//                  aParticleChange.ProposeTrackStatus(fSuspend);
	//        }
	//
	//
	//
	//
        ////////////////////////////////////////////////////////////////////////////////////
	//

	
	////////////////////////////////////////////////////////////////////////////////////
	//  The fast sim way - Ben J, Nov 2012
	////////////////////////////////////////////////////////////////////////////////////
	//
	//

	// We don't want to produce any trackable G4 secondaries
	aParticleChange.SetNumberOfSecondaries(0);

	
        // Retrieve the Scintillation Integral for this material  
        // new G4PhysicsOrderedFreeVector allocated to hold CII's
	

	// Some explanation for later improvements to scint yield code:
	//
	// What does G4 do here?
	//  It produces light in 2 steps, fast (scnt=1) then slow (scnt=2)
	//
	// The ratio of slow photons to fast photons is related	by the yieldratio
	//  parameter.  G4's poisson fluctuating scheme is a bit different to ours
	//  - we should check that they are equivalent.
	//
	// G4 poisson fluctuates the number of initial photons then divides them
	//  with a constant factor between fast + slow, whereas we poisson 
	//  fluctuate separateyly the fast and slow detection numbers.
	//
	
	// get the number of photons produced from the IonizationAndScintillation
	// singleton
	larg4::IonizationAndScintillation::Instance()->Reset(&aStep);
	double MeanNumberOfPhotons = larg4::IonizationAndScintillation::Instance()->NumberScintillationPhotons();
//  double stepEnergy          = larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit()/CLHEP::MeV;
        RecordPhotonsProduced(aStep, MeanNumberOfPhotons);//, stepEnergy);
	
	if (verboseLevel>0) {
	  G4cout << "\n Exiting from OpFastScintillation::DoIt -- NumberOfSecondaries = " 
		 << aParticleChange.GetNumberOfSecondaries() << G4endl;
	}
	
	
	return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


//-------------------------------------------------------------
  void OpFastScintillation::ProcessStep( const G4Step& step)
  {
    if(step.GetTotalEnergyDeposit() <= 0) return;

    OpDetPhotonTable::Instance()->AddEnergyDeposit
      (-1,
       -1,
       (double)(step.GetTotalEnergyDeposit()/CLHEP::MeV), //energy in MeV
       (float)(step.GetPreStepPoint()->GetPosition().x()/CLHEP::cm),
       (float)(step.GetPreStepPoint()->GetPosition().y()/CLHEP::cm),
       (float)(step.GetPreStepPoint()->GetPosition().z()/CLHEP::cm),
       (float)(step.GetPostStepPoint()->GetPosition().x()/CLHEP::cm),
       (float)(step.GetPostStepPoint()->GetPosition().y()/CLHEP::cm),
       (float)(step.GetPostStepPoint()->GetPosition().z()/CLHEP::cm),
       (double)(step.GetPreStepPoint()->GetGlobalTime()),
       (double)(step.GetPostStepPoint()->GetGlobalTime()),
       //step.GetTrack()->GetTrackID(),
       ParticleListAction::GetCurrentTrackID(),
       step.GetTrack()->GetParticleDefinition()->GetPDGEncoding(),
       step.GetPreStepPoint()->GetPhysicalVolume()->GetName()
       );
  }
  
//-------------------------------------------------------------

bool OpFastScintillation::RecordPhotonsProduced(const G4Step& aStep, double MeanNumberOfPhotons)//, double stepEnergy)
{
  // make sure that whatever happens afterwards, the energy deposition is stored
  art::ServiceHandle<sim::LArG4Parameters> lgp;
  if(lgp->FillSimEnergyDeposits())
    ProcessStep(aStep);
  
  
  // Get the pointer to the fast scintillation table
  OpDetPhotonTable * fst = OpDetPhotonTable::Instance();
  OpDetPhotonTable* litefst = OpDetPhotonTable::Instance();

  const G4Track * aTrack = aStep.GetTrack();

  G4StepPoint const* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint const* pPostStepPoint = aStep.GetPostStepPoint();
  
  const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
  const G4Material* aMaterial = aTrack->GetMaterial();

  G4int materialIndex = aMaterial->GetIndex();
  G4int tracknumber = aStep.GetTrack()->GetTrackID();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  //G4double      t0 = pPreStepPoint->GetGlobalTime() - fGlobalTimeOffset;
  G4double      t0 = pPreStepPoint->GetGlobalTime();
  
  
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();

  G4MaterialPropertyVector* Fast_Intensity = 
    aMaterialPropertiesTable->GetProperty("FASTCOMPONENT"); 
  G4MaterialPropertyVector* Slow_Intensity =
    aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");
  
  if (!Fast_Intensity && !Slow_Intensity )
    return 1;
  
  if(!bPropagate)
    return 0;

  // Get the visibility vector for this point
  art::ServiceHandle<phot::PhotonVisibilityService> pvs;
  size_t const NOpChannels = pvs->NOpChannels();


  G4int nscnt = 1;
  if (Fast_Intensity && Slow_Intensity) nscnt = 2;

  
  double Num = 0;
  double YieldRatio=0;

  
  if (scintillationByParticleType) {
    // The scintillation response is a function of the energy
    // deposited by particle types.
    
    // Get the definition of the current particle
    G4ParticleDefinition *pDef = aParticle->GetDefinition();
    
    // Obtain the G4MaterialPropertyVectory containing the
    // scintillation light yield as a function of the deposited
    // energy for the current particle type
    
    // Protons
    if(pDef==G4Proton::ProtonDefinition()) 
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("PROTONYIELDRATIO");

      }
    
    // Muons
    else if(pDef==G4MuonPlus::MuonPlusDefinition()||pDef==G4MuonMinus::MuonMinusDefinition())
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("MUONYIELDRATIO");
      }
    
    // Pions
    else if(pDef==G4PionPlus::PionPlusDefinition()||pDef==G4PionMinus::PionMinusDefinition())
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("PIONYIELDRATIO");
      }
    
    // Kaons
    else if(pDef==G4KaonPlus::KaonPlusDefinition()||pDef==G4KaonMinus::KaonMinusDefinition())
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("KAONYIELDRATIO");
      }
    
    // Alphas
    else if(pDef==G4Alpha::AlphaDefinition())
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("ALPHAYIELDRATIO");
      }
    
    // Electrons (must also account for shell-binding energy
    // attributed to gamma from standard PhotoElectricEffect)
    else if(pDef==G4Electron::ElectronDefinition() ||
	    pDef==G4Gamma::GammaDefinition())
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("ELECTRONYIELDRATIO");
      }
    
    // Default for particles not enumerated/listed above
    else
      {
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("ELECTRONYIELDRATIO");
      }
    
    // If the user has not specified yields for (p,d,t,a,carbon)
    // then these unspecified particles will default to the 
    // electron's scintillation yield
    if(YieldRatio==0){
      {
	
	YieldRatio = aMaterialPropertiesTable->
	  GetConstProperty("ELECTRONYIELDRATIO");
	
      }
    }
  }

  double const xyz[3] = { x0[0]/CLHEP::cm, x0[1]/CLHEP::cm, x0[2]/CLHEP::cm };
  float const* Visibilities = pvs->GetAllVisibilities(xyz);

  float const* ReflVisibilities = nullptr;
  float const* ReflT0s = nullptr;

  TF1* ParPropTimeTF1 = nullptr;

  if(pvs->StoreReflected()) {
    ReflVisibilities = pvs->GetAllVisibilities(xyz,true);
    if(pvs->StoreReflT0())
      ReflT0s = pvs->GetReflT0s(xyz); 
  }
  if(pvs->IncludeParPropTime())
  {
    ParPropTimeTF1 = pvs->GetTimingTF1(xyz);
  }
  /*
  // For Kazu to debug # photons generated using csv file, by default should be commented out 
  // but don't remove as it's useful. Separated portion of code relevant to this block
  // is labeled as "CASE-DEBUG DO NOT REMOVE THIS COMMENT"
  // (2017 May 2nd ... "kazuhiro at nevis dot columbia dot edu")
  //
  static bool first_time=false;
  if(!first_time) {
    std::cout<<"LOGMEid,pdg,mass,ke,te,de,x,y,z,t,dr,mean_npe,gen_npe,det_npe"<<std::endl;
    first_time=true;
  }

  std::cout<<"LOGME"
           <<aStep.GetTrack()->GetTrackID()<<","
           <<aParticle->GetDefinition()->GetPDGEncoding()<<","
           <<aParticle->GetDefinition()->GetPDGMass()/CLHEP::MeV<<","
           <<pPreStepPoint->GetKineticEnergy()<<","
           <<pPreStepPoint->GetTotalEnergy()<<","
           <<aStep.GetTotalEnergyDeposit()<<","
           <<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<","<<t0<<","
           <<aStep.GetDeltaPosition().mag()<<","
           <<MeanNumberOfPhotons<<","<<std::flush;

  double gen_photon_ctr=0;
  double det_photon_ctr=0;
  */  
  for (G4int scnt = 1; scnt <= nscnt; scnt++) {
    
    G4double ScintillationTime = 0.*CLHEP::ns;
    G4double ScintillationRiseTime = 0.*CLHEP::ns;
    G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;
    
    if (scnt == 1) {
      if (nscnt == 1) {
	if(Fast_Intensity){
	  ScintillationTime   = aMaterialPropertiesTable->
	    GetConstProperty("FASTTIMECONSTANT");
	  if (fFiniteRiseTime) {
	    ScintillationRiseTime = aMaterialPropertiesTable->
	      GetConstProperty("FASTSCINTILLATIONRISETIME");
	  }
	  ScintillationIntegral =
	    (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
	}
	if(Slow_Intensity){
	  ScintillationTime   = aMaterialPropertiesTable->
	    GetConstProperty("SLOWTIMECONSTANT");
	  if (fFiniteRiseTime) {
	    ScintillationRiseTime = aMaterialPropertiesTable->
	      GetConstProperty("SLOWSCINTILLATIONRISETIME");
	  }
	  ScintillationIntegral =
	    (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
	}
      }//endif nscnt=1
      else {
	if(YieldRatio==0) 
	  YieldRatio = aMaterialPropertiesTable->
	    GetConstProperty("YIELDRATIO");
	
	
	if ( ExcitationRatio == 1.0 ) {
	  Num = std::min(YieldRatio,1.0)*MeanNumberOfPhotons;
	}
	else {
	  Num = std::min(ExcitationRatio,1.0)*MeanNumberOfPhotons;
	}
	ScintillationTime   = aMaterialPropertiesTable->
		  GetConstProperty("FASTTIMECONSTANT");
	if (fFiniteRiseTime) {
	  ScintillationRiseTime = aMaterialPropertiesTable->
	    GetConstProperty("FASTSCINTILLATIONRISETIME");
	}
	ScintillationIntegral =
	  (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
      }//endif nscnt!=1
    }//end scnt=1
    
    else {
      Num = MeanNumberOfPhotons - Num;
      ScintillationTime   =   aMaterialPropertiesTable->
	GetConstProperty("SLOWTIMECONSTANT");
      if (fFiniteRiseTime) {
                    ScintillationRiseTime = aMaterialPropertiesTable->
		      GetConstProperty("SLOWSCINTILLATIONRISETIME");
      }
      ScintillationIntegral =
	(G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
    }
    
    if (!ScintillationIntegral) continue;

    //gen_photon_ctr += Num; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
    
    // Max Scintillation Integral
    
    //            G4double CIImax = ScintillationIntegral->GetMaxValue();
    
    
    //std::cout << "++++++++++++" << Num << "++++++++++" << std::endl;
    

    // here we go: now if visibilities are invalid, we are in trouble
    //if (!Visibilities && (NOpChannels > 0)) {
    //  throw cet::exception("OpFastScintillator")
    //    << "Photon library does not cover point ( " << xyz[0] << ", "
    //    << xyz[1] << ", " << xyz[2] << " ) cm.\n";
    //}
    
    if(!Visibilities){
    }else{
      std::map<int, int> DetectedNum;

      std::map<int, int> ReflDetectedNum;

      for(size_t OpDet=0; OpDet!=NOpChannels; OpDet++)
      {
    		G4int DetThisPMT = G4int(G4Poisson(Visibilities[OpDet] * Num));

    		if(DetThisPMT>0) 
        {
		      DetectedNum[OpDet]=DetThisPMT;
		      //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
		      //   //   it->second<<" " << Num << " " << DetThisPMT;  

		      //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
        }
		if(pvs->StoreReflected()) {
		  G4int ReflDetThisPMT = G4int(G4Poisson(ReflVisibilities[OpDet] * Num));
		  if(ReflDetThisPMT>0)
		    {
		      ReflDetectedNum[OpDet]=ReflDetThisPMT;
		    }
                }

      }
	    // Now we run through each PMT figuring out num of detected photons
	
      if(lgp->UseLitePhotons())
      {
        std::map<int, std::map<int, int>> StepPhotonTable;
        // And then add these to the total collection for the event     
        for(std::map<int,int>::const_iterator itdetphot = DetectedNum.begin();
                itdetphot!=DetectedNum.end(); ++itdetphot)
        {
          std::map<int, int>  StepPhotons;
          for (G4int i = 0; i < itdetphot->second; ++i)
          {
            G4double deltaTime = aStep.GetStepLength() /
                ((pPreStepPoint->GetVelocity()+ pPostStepPoint->GetVelocity())/2.);


            if (ScintillationRiseTime==0.0) {
                deltaTime = deltaTime -
                    ScintillationTime * std::log( G4UniformRand() );
            } else {
                deltaTime = deltaTime +
                    sample_time(ScintillationRiseTime, ScintillationTime);
            }
            //Warning: TF1::GetNdim()==1 will tell us if the TF1 is really defined or it is the default one.
            //This will fix a segfault when using timing and interpolation.
	    if(pvs->IncludeParPropTime()&&ParPropTimeTF1&&ParPropTimeTF1[itdetphot->first].GetNdim()==1) {
              deltaTime += ParPropTimeTF1[itdetphot->first].GetRandom();
	    }

            G4double aSecondaryTime = t0 + deltaTime;
            float Time = aSecondaryTime;
            int ticks = static_cast<int>(Time);
            StepPhotons[ticks]++;
          }
         StepPhotonTable[itdetphot->first] = StepPhotons;
         //Iterate over Step Photon Table to add photons to OpDetBacktrackerRecords.

         sim::OpDetBacktrackerRecord tmpOpDetBTRecord(itdetphot->first);
         //int thisG4TrackID = (aStep.GetTrack())->GetTrackID();
         int thisG4TrackID = ParticleListAction::GetCurrentTrackID();
         CLHEP::Hep3Vector prePoint  = (aStep.GetPreStepPoint())->GetPosition();
         CLHEP::Hep3Vector postPoint = (aStep.GetPostStepPoint())->GetPosition();
         //Note the use of xO (letter O) instead of x0. This is to differentiate the positions here with the earlier declared double* x0
         double xO = ( ( (prePoint.getX() + postPoint.getX() ) / 2.0) / CLHEP::cm );
         double yO = ( ( (prePoint.getY() + postPoint.getY() ) / 2.0) / CLHEP::cm );
         double zO = ( ( (prePoint.getZ() + postPoint.getZ() ) / 2.0) / CLHEP::cm );
         double const xyzPos[3] = {xO,yO,zO};
//         double energy  = ( aStep.GetTotalEnergyDeposit() / CLHEP::MeV );
         double energy  = larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit()/CLHEP::MeV;
         //Loop over StepPhotons to get number of photons detected at each time for this channel and G4Step.
         for(std::map<int,int>::iterator stepPhotonsIt = StepPhotons.begin(); stepPhotonsIt != StepPhotons.end(); ++stepPhotonsIt)
         {
           int photonTime = stepPhotonsIt->first;
           int numPhotons = stepPhotonsIt->second;
           tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, photonTime, numPhotons, xyzPos, energy);
         }
         //Add OpDetBackTrackerRecord. (opdetphotonTABLE->instance().addOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord BTRrecord)
         litefst->AddOpDetBacktrackerRecord(tmpOpDetBTRecord);
        }
        litefst->AddPhoton(&StepPhotonTable);
      }
      else
       {
	 // We need to know the positions of the different optical detectors                                                      
	 art::ServiceHandle<geo::Geometry> geo;
	 // And then add these to the total collection for the event	    
	 for(std::map<int,int>::const_iterator itdetphot = DetectedNum.begin();
	     itdetphot!=DetectedNum.end(); ++itdetphot)
	   {
	     G4double propagation_time = 0;
	     std::vector<double> arrival_time_dist_vuv;
	     if(pvs->IncludePropTime()) {
	       // Get VUV photons arrival time distribution from the parametrization 
	       double OpDetCenter[3];
	       geo->OpDetGeoFromOpDet(itdetphot->first).GetCenter(OpDetCenter);
	       G4ThreeVector OpDetPoint(OpDetCenter[0]*CLHEP::cm,OpDetCenter[1]*CLHEP::cm,OpDetCenter[2]*CLHEP::cm);
	       double distance_in_cm = (x0 - OpDetPoint).mag()/CLHEP::cm; // this must be in CENTIMETERS! 
	       arrival_time_dist_vuv = GetVUVTime(distance_in_cm, itdetphot->second);//in ns
	       if(!(arrival_time_dist_vuv.size()>0)) {
		 //G4cout<<"Number of VUV detected photons = "<<itdetphot->second<<" ... too small, => No photons simulated!"<<G4endl;
		 continue;
	       }
	     }
	     
	        for (G4int i = 0; i < itdetphot->second; ++i) 
          {
            G4double deltaTime = aStep.GetStepLength() /
                ((pPreStepPoint->GetVelocity()+
                  pPostStepPoint->GetVelocity())/2.);
		
            if (ScintillationRiseTime==0.0) 
            {
                deltaTime = deltaTime - 
                    ScintillationTime * std::log( G4UniformRand() );
            } 
            else 
            {
              deltaTime = deltaTime +
                  sample_time(ScintillationRiseTime, ScintillationTime);
            }		
		
            G4double aSecondaryTime = t0 + deltaTime;
	    if(pvs->IncludePropTime())
	      propagation_time = arrival_time_dist_vuv.at(i)*CLHEP::ns;
	
            // The sim photon in this case stores its production point and time
            TVector3 PhotonPosition(x0[0],x0[1],x0[2]);
		
            // We don't know anything about the momentum dir, so set it to be Z		
            float Energy = 9.7*CLHEP::eV;
            float Time = aSecondaryTime;
	    if(pvs->IncludePropTime())
	      Time += propagation_time;
	    
            // Make a photon object for the collection
	    sim::OnePhoton PhotToAdd;
            PhotToAdd.InitialPosition  = PhotonPosition;
            PhotToAdd.Energy           = Energy;
            PhotToAdd.Time             = Time;
            PhotToAdd.SetInSD          = false;
	    PhotToAdd.MotherTrackID    = tracknumber;

            fst->AddPhoton(itdetphot->first, std::move(PhotToAdd));
          }
        }
	 // repeat the same for the reflected light, in case it has been sored
	 if(pvs->StoreReflected())
	   {
	     for(std::map<int,int>::const_iterator itdetphot = ReflDetectedNum.begin();
		 itdetphot!=ReflDetectedNum.end(); ++itdetphot)
	       {
		 G4double propagation_time = 0;
		 std::vector<double> arrival_time_dist_vis;
		 if(pvs->IncludePropTime()) {
		   // Get Visible photons arrival time distribution from the parametrization
		   double OpDetCenter[3];
		   geo->OpDetGeoFromOpDet(itdetphot->first).GetCenter(OpDetCenter);
		   G4ThreeVector OpDetPoint(OpDetCenter[0]*CLHEP::cm,OpDetCenter[1]*CLHEP::cm,OpDetCenter[2]*CLHEP::cm);
		   // currently unused:
		 //  double distance_in_cm = (x0 - OpDetPoint).mag()/CLHEP::cm; // this must be in CENTIMETERS!
		   double t0_vis_lib = ReflT0s[itdetphot->first];		   
		   arrival_time_dist_vis = GetVisibleTimeOnlyCathode(t0_vis_lib, itdetphot->second);
		   if(!(arrival_time_dist_vis.size()>0)) {
		     //G4cout<<"Number of Visible detected photons = "<<itdetphot->second<<" ... too small, => No photons simulated!"<<G4endl; 
		     continue;
		   }
		 }
		 for (G4int i = 0; i < itdetphot->second; ++i)
		   {
		     G4double deltaTime = aStep.GetStepLength() /
		       ((pPreStepPoint->GetVelocity()+
			 pPostStepPoint->GetVelocity())/2.);

		     if (ScintillationRiseTime==0.0)
		       {
                deltaTime = deltaTime -
		  ScintillationTime * std::log( G4UniformRand() );
		       }
		     else
		       {
              deltaTime = deltaTime +
		sample_time(ScintillationRiseTime, ScintillationTime);
		       }
		     G4double aSecondaryTime = t0 + deltaTime;
		     if(pvs->IncludePropTime())		   
		       propagation_time = arrival_time_dist_vis.at(i)*CLHEP::ns;

		     TVector3 PhotonPosition(x0[0],x0[1],x0[2]);

		     // We don't know anything about the momentum dir, so set it to be Z                                                                     
		     std::map<double,double> tpbemission=lar::providerFrom<detinfo::LArPropertiesService>()->TpbEm();
		     const int nbins=tpbemission.size();

		     double * parent=new double[nbins];

		     int ii=0;
		     for( std::map<double, double>::iterator iter = tpbemission.begin(); iter != tpbemission.end(); ++iter)
		       { parent[ii++]=(*iter).second;
		       }

		     CLHEP::RandGeneral rgen0(parent,nbins);

		     double energy_reemited = rgen0.fire()*((*(--tpbemission.end())).first-(*tpbemission.begin()).first)+(*tpbemission.begin()).first;

		     float Energy = energy_reemited*CLHEP::eV;
		     float Time = aSecondaryTime;
		     if(pvs->IncludePropTime())
		       Time += propagation_time;

		     // Make a photon object for the collection  
		     sim::OnePhoton PhotToAdd;
		     PhotToAdd.InitialPosition  = PhotonPosition;
		     PhotToAdd.Energy           = Energy;
		     PhotToAdd.Time             = Time;
		     PhotToAdd.MotherTrackID    = tracknumber;
		     PhotToAdd.SetInSD          = false;
		     fst->AddPhoton(itdetphot->first, std::move(PhotToAdd));
		     delete [] parent;
		   }     
	       }
	   }
      }
    }
  }

  //std::cout<<gen_photon_ctr<<","<<det_photon_ctr<<std::endl; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
  return 0;
}


// BuildThePhysicsTable for the scintillation process
// --------------------------------------------------
//

void OpFastScintillation::BuildThePhysicsTable()
{
        if (theFastIntegralTable && theSlowIntegralTable) return;

        const G4MaterialTable* theMaterialTable = 
                               G4Material::GetMaterialTable();
        G4int numOfMaterials = G4Material::GetNumberOfMaterials();

        // create new physics table
	
        if(!theFastIntegralTable)theFastIntegralTable = new G4PhysicsTable(numOfMaterials);
        if(!theSlowIntegralTable)theSlowIntegralTable = new G4PhysicsTable(numOfMaterials);

        // loop for materials

        for (G4int i=0 ; i < numOfMaterials; i++)
        {
                G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
					new G4PhysicsOrderedFreeVector();
                G4PhysicsOrderedFreeVector* bPhysicsOrderedFreeVector =
                                        new G4PhysicsOrderedFreeVector();

                // Retrieve vector of scintillation wavelength intensity for
                // the material from the material's optical properties table.

                G4Material* aMaterial = (*theMaterialTable)[i];

                G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                aMaterial->GetMaterialPropertiesTable();

                if (aMaterialPropertiesTable) {

                   G4MaterialPropertyVector* theFastLightVector = 
                   aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");

                   if (theFastLightVector) {

                      // Retrieve the first intensity point in vector
                      // of (photon energy, intensity) pairs 

                      G4double currentIN = (*theFastLightVector)[0];

                      if (currentIN >= 0.0) {

                         // Create first (photon energy, Scintillation 
                         // Integral pair  

                         G4double currentPM = theFastLightVector->Energy(0);

                         G4double currentCII = 0.0;

                         aPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon energy, intensity)
                         // pairs stored for this material  

                         for (size_t i = 1;
                              i < theFastLightVector->GetVectorLength();
                              i++)
                         {
                                currentPM = theFastLightVector->Energy(i);
                                currentIN = (*theFastLightVector)[i];

                                currentCII = 0.5 * (prevIN + currentIN);

                                currentCII = prevCII +
                                             (currentPM - prevPM) * currentCII;

                                aPhysicsOrderedFreeVector->
                                    InsertValues(currentPM, currentCII);

                                prevPM  = currentPM;
                                prevCII = currentCII;
                                prevIN  = currentIN;
                         }

                      }
                   }

                   G4MaterialPropertyVector* theSlowLightVector =
                   aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

                   if (theSlowLightVector) {

                      // Retrieve the first intensity point in vector
                      // of (photon energy, intensity) pairs

                      G4double currentIN = (*theSlowLightVector)[0];

                      if (currentIN >= 0.0) {

                         // Create first (photon energy, Scintillation
                         // Integral pair

                         G4double currentPM = theSlowLightVector->Energy(0);

                         G4double currentCII = 0.0;

                         bPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon energy, intensity)
                         // pairs stored for this material

                         for (size_t i = 1;
                              i < theSlowLightVector->GetVectorLength();
                              i++)
                         {
                                currentPM = theSlowLightVector->Energy(i);
                                currentIN = (*theSlowLightVector)[i];

                                currentCII = 0.5 * (prevIN + currentIN);

                                currentCII = prevCII +
                                             (currentPM - prevPM) * currentCII;

                                bPhysicsOrderedFreeVector->
                                    InsertValues(currentPM, currentCII);

                                prevPM  = currentPM;
                                prevCII = currentCII;
                                prevIN  = currentIN;
                         }

                      }
                   }
                }

        // The scintillation integral(s) for a given material
        // will be inserted in the table(s) according to the
        // position of the material in the material table.

        theFastIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
        theSlowIntegralTable->insertAt(i,bPhysicsOrderedFreeVector);

        }
}

// Called by the user to set the scintillation yield as a function
// of energy deposited by particle type

void OpFastScintillation::SetScintillationByParticleType(const G4bool scintType)
{
        if (emSaturation) {
           G4Exception("OpFastScintillation::SetScintillationByParticleType", "Scint02",
                       JustWarning, "Redefinition: Birks Saturation is replaced by ScintillationByParticleType!");
           RemoveSaturation();
        }
        scintillationByParticleType = scintType;
}

// GetMeanFreePath
// ---------------
//

G4double OpFastScintillation::GetMeanFreePath(const G4Track&,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;

        return DBL_MAX;

}

// GetMeanLifeTime
// ---------------
//

G4double OpFastScintillation::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
        *condition = Forced;

        return DBL_MAX;

}

G4double OpFastScintillation::sample_time(G4double tau1, G4double tau2)
{
// tau1: rise time and tau2: decay time

        while(1) {
          // two random numbers
          G4double ran1 = G4UniformRand();
          G4double ran2 = G4UniformRand();
          //
          // exponential distribution as envelope function: very efficient
          //
          G4double d = (tau1+tau2)/tau2;
          // make sure the envelope function is 
          // always larger than the bi-exponential
          G4double t = -1.0*tau2*std::log(1-ran1);
          G4double g = d*single_exp(t,tau2);
          if (ran2 <= bi_exp(t,tau1,tau2)/g) return t;
        }
        return -1.0;
}

//                         ======TIMING PARAMETRIZATION=====           //

// Parametrization of the VUV light timing (result from direct transport + Rayleigh scattering ONLY)
// using a landau + expo function.The function below returns the arrival time distribution given the
// distance IN CENTIMETERS between the scintillation/ionization point and the optical detectotr.
std::vector<double> OpFastScintillation::GetVUVTime(double distance, int number_photons) {

  //-----Distances in cm and times in ns-----// 
  //gRandom->SetSeed(0);
  std::vector<double> arrival_time_distrb;
  arrival_time_distrb.clear();

  // Parametrization data: 

  if(distance < 10 || distance > fd_max) {
    G4cout<<"WARNING: Parametrization of Direct Light not fully reliable"<<G4endl;
    G4cout<<"Too close/far to the PMT  -> set 0 VUV photons(?)!!!!!!"<<G4endl;
    return arrival_time_distrb;
  }
  //signals (remember this is transportation) no longer than 1us  
  const double signal_t_range = 1000.;
  const double vuv_vgroup = 10.13;//cm/ns     
  double t_direct = distance/vuv_vgroup;
  // Defining the two functions (Landau + Exponential) describing the timing vs distance  
  double pars_landau[3] = {functions_vuv[1]->Eval(distance), functions_vuv[2]->Eval(distance),
			   pow(10.,functions_vuv[0]->Eval(distance))};
  if(distance > fd_break) {
    pars_landau[0]=functions_vuv[6]->Eval(distance);
    pars_landau[1]=functions_vuv[2]->Eval(fd_break);
    pars_landau[2]=pow(10.,functions_vuv[5]->Eval(distance));
  }
  TF1 *flandau = new TF1("flandau","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
  flandau->SetParameters(pars_landau);
  double pars_expo[2] = {functions_vuv[3]->Eval(distance), functions_vuv[4]->Eval(distance)};
  if(distance > (fd_break - 50.)) {
    pars_expo[0] = functions_vuv[7]->Eval(distance);
    pars_expo[1] = functions_vuv[4]->Eval(fd_break - 50.);
  }
  TF1 *fexpo = new TF1("fexpo","expo",0, signal_t_range/2);
  fexpo->SetParameters(pars_expo);
  //this is to find the intersection point between the two functions:      
  TF1 *fint = new TF1("fint",finter_d,flandau->GetMaximumX(),3*t_direct,5);
  double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
  fint->SetParameters(parsInt);
  double t_int = fint->GetMinimumX();
  double minVal = fint->Eval(t_int);
  //the functions must intersect!!!     
  if(minVal>0.015)
    G4cout<<"WARNING: Parametrization of Direct Light discontinuous (landau + expo)!!!!!!"<<G4endl;

  TF1 *fVUVTiming =  new TF1("fTiming",LandauPlusExpoFinal,0,signal_t_range,6);
  double parsfinal[6] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
  fVUVTiming->SetParameters(parsfinal);
  // Set the number of points used to sample the function                         

  int f_sampling = 1000;
  if(distance < 50)
    f_sampling *= ftf1_sampling_factor;
  fVUVTiming->SetNpx(f_sampling);

  for(int i=0; i<number_photons; i++)
    arrival_time_distrb.push_back(fVUVTiming->GetRandom());

  //deleting ...                                                                                     
  delete flandau;
  delete fexpo;
  delete fint;
  delete fVUVTiming;
  G4cout<<"BEAMAUS timing distribution hecha"<<G4endl;
  return arrival_time_distrb;
}

// Parametrization of the Visible light timing (result from direct transport + Rayleigh scattering ONLY) 
// using a landau + exponential function. The function below returns the arrival time distribution given the 
// time of the first visible photon in the PMT. The light generated has been reflected by the cathode ONLY.
std::vector<double> OpFastScintillation::GetVisibleTimeOnlyCathode(double t0, int number_photons){
  //-----Distances in cm and times in ns-----//  
  //gRandom->SetSeed(0);  
               
  std::vector<double> arrival_time_distrb;
  arrival_time_distrb.clear();
  // Parametrization data:

  if(t0 < 8 || t0 > ft0_max) {
    G4cout<<"WARNING: Parametrization of Cathode-Only reflected Light not fully reliable"<<G4endl;
    G4cout<<"Too close/far to the PMT  -> set 0 Visible photons(?)!!!!!!"<<G4endl;
    return arrival_time_distrb;
  }
  //signals (remember this is transportation) no longer than 1us     
  const double signal_t_range = 1000.;
  double pars_landau[3] = {functions_vis[1]->Eval(t0), functions_vis[2]->Eval(t0),
                           pow(10.,functions_vis[0]->Eval(t0))};
  double pars_expo[2] = {functions_vis[3]->Eval(t0), functions_vis[4]->Eval(t0)};
  if(t0 > ft0_break_point) {
    pars_landau[0] = -0.798934 + 1.06216*t0;
    pars_landau[1] = functions_vis[2]->Eval(ft0_break_point);
    pars_landau[2] = pow(10.,functions_vis[0]->Eval(ft0_break_point));
    pars_expo[0] = functions_vis[3]->Eval(ft0_break_point);
    pars_expo[1] = functions_vis[4]->Eval(ft0_break_point);
  }

  // Defining the two functions (Landau + Exponential) describing the timing vs t0                  
  TF1 *flandau = new TF1("flandau","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
  flandau->SetParameters(pars_landau);
  TF1 *fexpo = new TF1("fexpo","expo",0, signal_t_range/2);
  fexpo->SetParameters(pars_expo);
  //this is to find the intersection point between the two functions:                                                          
  TF1 *fint = new TF1("fint",finter_d,flandau->GetMaximumX(),2*t0,5);
  double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
  fint->SetParameters(parsInt);
  double t_int = fint->GetMinimumX();
  double minVal = fint->Eval(t_int);
  //the functions must intersect!!!                                             
  if(minVal>0.015)
    G4cout<<"WARNING: Parametrization of Direct Light discontinuous (landau + expo)!!!!!!"<<G4endl;

  TF1 *fVisTiming =  new TF1("fTiming",LandauPlusExpoFinal,0,signal_t_range,6);
  double parsfinal[6] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
  fVisTiming->SetParameters(parsfinal);
  // Set the number of points used to sample the function                                  

  int f_sampling = 1000;
  if(t0 < 20)
    f_sampling *= ftf1_sampling_factor;
  fVisTiming->SetNpx(f_sampling);

  for(int i=0; i<number_photons; i++)
    arrival_time_distrb.push_back(fVisTiming->GetRandom());

  //deleting ...                                    

  delete flandau;
  delete fexpo;
  delete fint;
  delete fVisTiming;
  G4cout<<"Timing distribution BEAMAUS!"<<G4endl;
  return arrival_time_distrb;
}

double finter_d(double *x, double *par) {

  double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
  double y2 = TMath::Exp(par[3]+x[0]*par[4]);

  return TMath::Abs(y1 - y2);
}
double LandauPlusExpoFinal(double *x, double *par)
{
  // par0 = joining point
  // par1 = Landau MPV 
  // par2 = Landau widt
  // par3 = normalization 
  // par4 = Expo cte  
  // par5 = Expo tau 
  double y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
  double y2 = TMath::Exp(par[4]+x[0]*par[5]);
  if(x[0] > par[0]) y1 = 0.;
  if(x[0] < par[0]) y2 = 0.;

  return (y1 + y2);
}


double finter_r(double *x, double *par) {

  double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
  double y2 = par[5]*TMath::Landau(x[0],par[3],par[4]);

  return TMath::Abs(y1 - y2);
}



}
