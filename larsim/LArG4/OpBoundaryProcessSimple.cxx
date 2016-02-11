// Adapted for LArSoft by Ben Jones, MIT, March 2010
//
// This class invokes a simplified model of optical reflections at
// boundaries between different materials.  The relevant reflectivities
// are read in by an implementation of the MaterialPropertiesLoader.
//
// The required parameters are total reflectance and radio of diffuse
// to specular reflectance.  Each photon crossing a boundary with a
// defined reflectance is randomly either reflected or absorbed and killed
// according to the supplied probability.
//
// Every reflected photon with a defined diffuse reflection fraction
// is then randomly either diffusely or specularly reflected according
// to the supplied probability.  All materials with no defined
// reflectance are assumed to be black and
//  absorb all incident photons.
//
// This physics process is loaded in the OpticalPhysics physics constructor.
//
// This class is based on the G4OpBoundaryProcess class in Geant4


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
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Implementation
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_OPBOUNDARYPROCESSSIMPLE_CXX
#define LARG4_OPBOUNDARYPROCESSSIMPLE_CXX 1

#include "TH1.h"
#include "TTree.h"

#include "Geant4/G4ios.hh"
#include "Geant4/G4OpProcessSubType.hh"
#include "Geant4/G4GeometryTolerance.hh"

#include "larsim/LArG4/OpBoundaryProcessSimple.hh"
#include "larsim/Simulation/LArG4Parameters.h"

//#define G4DEBUG_OPTICAL



namespace larg4 {

  //Constructor

  OpBoundaryProcessSimple::OpBoundaryProcessSimple(const G4String& processName,
               G4ProcessType type)
    : G4VDiscreteProcess(processName, type)
  {
    if ( 1 > 0) {
      G4cout << GetProcessName() << " is created " << G4endl;
    }


    SetProcessSubType(fOpBoundary);

    fTheStatus = Undefined;

    fCarTolerance = G4GeometryTolerance::GetInstance()
      ->GetSurfaceTolerance();

  }

  // OpBoundaryProcessSimple::OpBoundaryProcessSimple(const OpBoundaryProcessSimple &right)
  // {
  // }


  //Destructor

  OpBoundaryProcessSimple::~OpBoundaryProcessSimple()
  {
  }




  //Action to take after making each step of an optical photon - described in file header.

  // PostStepDoIt
  G4VParticleChange*
  OpBoundaryProcessSimple::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
  {

    // Note - these have to be loaded here, since this object
    // is constructed before LArG4, and hence before the
    // LArG4 parameters are read from config

    fTheStatus = Undefined;
    aParticleChange.Initialize(aTrack);


    G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();


    if (pPostStepPoint->GetStepStatus() != fGeomBoundary){
      fTheStatus = NotAtBoundary;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    if (aTrack.GetStepLength()<=fCarTolerance/2){
      fTheStatus = StepTooSmall;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    G4Material* Material1 = pPreStepPoint  -> GetMaterial();
    G4Material* Material2 = pPostStepPoint -> GetMaterial();


    if(Material1 != Material2)
      {

  art::ServiceHandle<sim::LArG4Parameters> lgp;

  fVerbosity = lgp->OpVerbosity();



  G4double thePhotonMomentum;

  G4ThreeVector OldMomentum;
  G4ThreeVector OldPolarization;

  G4ThreeVector NewMomentum;
  G4ThreeVector NewPolarization;

  G4ThreeVector theGlobalNormal;
  G4ThreeVector theFacetNormal;

  std::string Material1Name = Material1->GetName();
  std::string Material2Name = Material2->GetName();

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();


  thePhotonMomentum = aParticle->GetTotalMomentum();
  OldMomentum       = aParticle->GetMomentumDirection();
  OldPolarization   = aParticle->GetPolarization();



  if (fVerbosity>9) std::cout<<"OpBoundaryProcessSimple Debug: Photon " <<aTrack.GetTrackID() <<" momentum "<<thePhotonMomentum << " Material1 " << Material1Name.c_str() << " Material 2 " << Material2Name.c_str() << std::endl;

  G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();

  G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();

  G4ThreeVector theLocalPoint = theNavigator->
    GetGlobalToLocalTransform().
    TransformPoint(theGlobalPoint);

  G4ThreeVector theLocalNormal;   // Normal points back into volume

  G4bool valid;
  theLocalNormal = theNavigator->GetLocalExitNormal(&valid);

  if (valid) {
    theLocalNormal = -theLocalNormal;
  }
  else {
    G4cerr << " OpBoundaryProcessSimple/PostStepDoIt(): "
     << " The Navigator reports that it returned an invalid normal"
     << G4endl;
  }

  theGlobalNormal = theNavigator->GetLocalToGlobalTransform().
    TransformAxis(theLocalNormal);

  if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4DEBUG_OPTICAL
    G4cerr << " OpBoundaryProcessSimple/PostStepDoIt(): "
     << " theGlobalNormal points the wrong direction "
     << G4endl;
#endif
    theGlobalNormal = -theGlobalNormal;
  }

  G4MaterialPropertiesTable* aMaterialPropertiesTable;
  G4MaterialPropertyVector* Reflectance;
  G4MaterialPropertyVector* DiffuseReflectanceFraction;


  aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();


  if (aMaterialPropertiesTable) {

   std::stringstream PropertyName;
   PropertyName<<"REFLECTANCE_"<<Material2Name.c_str();
   Reflectance =         aMaterialPropertiesTable->GetProperty(PropertyName.str().c_str());

   PropertyName.str("");
   PropertyName<<"DIFFUSE_REFLECTANCE_FRACTION_"<<Material2Name.c_str();
   DiffuseReflectanceFraction=aMaterialPropertiesTable->GetProperty(PropertyName.str().c_str());

    if (Reflectance)
      {
        double theReflectance = Reflectance->Value(thePhotonMomentum);
        if( G4BooleanRand(theReflectance))
    {
      if(DiffuseReflectanceFraction)
        {
          double theDiffuseReflectanceFraction = DiffuseReflectanceFraction->Value(thePhotonMomentum);
          if(G4BooleanRand(theDiffuseReflectanceFraction))
      {
        fTheStatus=SimpleDiffuse;
      }
          else
      {
        fTheStatus=SimpleSpecular;
      }
        }
      else
        {
          fTheStatus=SimpleSpecular;
        }
    }
        else
    {
        fTheStatus=SimpleAbsorbed;
    }

      }
    else
      {
        fTheStatus=SimpleAbsorbedNoRefl;
      }
  }




  // Take action according to track status

  if ( fTheStatus == SimpleDiffuse ) {
    NewMomentum = G4LambertianRand(theGlobalNormal);
    theFacetNormal = (NewMomentum - OldMomentum).unit();
  }

  else if (fTheStatus==SimpleSpecular) {
    theFacetNormal = theGlobalNormal;
    G4double PdotN = OldMomentum * theFacetNormal;
    NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
  }

  else if ((fTheStatus==SimpleAbsorbed) || (fTheStatus==SimpleAbsorbedNoRefl)) {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  else if ((fTheStatus==NotAtBoundary)||(fTheStatus==StepTooSmall))
    {
      NewMomentum=OldMomentum;
    }
  else aParticleChange.ProposeTrackStatus(fStopAndKill);

  NewMomentum = NewMomentum.unit();
  aParticleChange.ProposeMomentumDirection(NewMomentum);
      }
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);

  }



  // Compulsary method for G4DiscreteProcesses.  Serves no actual function here.

  // GetMeanFreePath
  G4double OpBoundaryProcessSimple::GetMeanFreePath(const G4Track& ,
                G4double ,
                G4ForceCondition* condition)
  {
    *condition = Forced;

    return DBL_MAX;
  }
}


#endif


