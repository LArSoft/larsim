// Class ageodapted for LArSoft by Ben Jones, MIT 10/10/12
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
#include "TFormula.h"
//#include "Math/SpecFuncMathMore.h"
#include "TSystem.h"
#include <cmath>

#include "boost/algorithm/string.hpp"

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

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
      // Loading the position of each optical channel, neccessary for the parametrizatiuons of Nhits and prop-time 
      static art::ServiceHandle<geo::Geometry> geo;
   
      for(size_t i = 0; i != pvs->NOpChannels(); i++)
	{
	  double OpDetCenter_i[3];
	  std::vector<double> OpDetCenter_v;
	  geo->OpDetGeoFromOpDet(i).GetCenter(OpDetCenter_i);
	  OpDetCenter_v.assign(OpDetCenter_i, OpDetCenter_i +3);
	  fOpDetCenter.push_back(OpDetCenter_v);
	}
      
      if(pvs->IncludePropTime()) {
	std::cout << "Using parameterisation of timings." << std::endl;
        //OLD VUV time parapetrization (to be removed soon)
        pvs->SetDirectLightPropFunctions(functions_vuv, fd_break, fd_max, ftf1_sampling_factor);
        pvs->SetReflectedCOLightPropFunctions(functions_vis, ft0_max, ft0_break_point);
	//New VUV time parapetrization
	pvs->LoadTimingsForVUVPar(fparameters, fstep_size, fmax_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance);

	// create vector of empty TF1s that will be replaces with the parameterisations that are generated as they are required
	// default TF1() constructor gives function with 0 dimensions, can then check numDim to qucikly see if a parameterisation has been generated  
	int num_params = (fmax_d - 25) / fstep_size;  // for d < 25cm, no parameterisaton, a delta function is used instead
	std::vector<TF1> VUV_timing_temp(num_params,TF1());
	VUV_timing = VUV_timing_temp;
    
	// initialise vectors to contain range parameterisations sampled to in each case
	// when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling is regenerated, this is the slow part!
	std::vector<double> VUV_empty(num_params, 0);
	VUV_max = VUV_empty;
	VUV_min = VUV_empty;

        // VIS time parameterisation
        if (pvs->StoreReflected()) {
	  // load parameters
	  pvs->LoadTimingsForVISPar(fdistances_refl, fcut_off_pars, ftau_pars, fvis_vmean, fn_LAr_vis, fn_LAr_vuv, fplane_depth);	 	  
	}

      }
      if(pvs->UseNhitsModel()) {
	std::cout << "Using semi-analytic model for number of hits:" << std::endl;	
	
	// LAr absorption length in cm
	std::map<double, double> abs_length_spectrum = lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
	std::vector<double> x_v, y_v;
	for(auto elem : abs_length_spectrum) {
	  x_v.push_back(elem.first);
	  y_v.push_back(elem.second);
	}
	fL_abs_vuv =  interpolate(x_v, y_v, 9.7, false);

	// Load Gaisser-Hillas corrections for VUV semi-analytic hits 
	std::cout<<"Loading the GH corrections"<<std::endl;
	pvs->LoadGHForVUVCorrection(fGHvuvpars, fzdimension, fydimension, fradius, foptical_detector_type);
        fdelta_angulo = 10.; // angle bin size

	// initialise gaisser hillas functions for VUV Rayleigh scattering correction
	double pars_ini[4] = {0,0,0,0};
	std::cout << "Initalising Gaisser-Hillas parameters" << std::endl;
	for(int bin = 0; bin < 9; bin++) {
	  GHvuv[bin] =  new TF1("GH",GaisserHillas,0.,2000,4);      
	  for(int j=0; j < 4; j++) {
	    // loads parameter set read in from fcl
	    pars_ini[j] = fGHvuvpars[j][bin];
	  }
	  GHvuv[bin]->SetParameters(pars_ini);
	}
	
	if(pvs->StoreReflected()) {
	  // Load corrections for VIS semi-anlytic hits
	  std::cout << "Loading vis corrections"<<std::endl;
	  pvs->LoadParsForVISCorrection(fvispars, fplane_depth, fcathode_zdimension, fcathode_ydimension, fcathode_centre, fzdimension,fydimension,fradius, foptical_detector_type);
	 std::cout << "fzdimension = " << fzdimension << ", fydimension = " << fydimension << ", fcathode_ydimension = " << fcathode_ydimension << ", fcathode_zdimension = " << fcathode_zdimension << std::endl;
	  // initialise vis correction functions
	  double pars_ini_vis[6] = {0,0,0,0,0,0};
	  std::cout << "Initialising visible correction parameters" << std::endl;
	  for (int bin = 0; bin < 9; bin++) {
	    VIS_pol[bin] = new TF1 ("pol", "pol5", 0, 2000);
	    for (int j = 0; j < 6; j++){
	      // loads paramter set read in from fcl
	      pars_ini_vis[j] = fvispars[j][bin];
	    }
	    VIS_pol[bin]->SetParameters(pars_ini_vis);
	  }	
	}
      }
    }
    tpbemission=lar::providerFrom<detinfo::LArPropertiesService>()->TpbEm();
    const int nbins=tpbemission.size();
    double * parent=new double[nbins];
    int ii=0;
    for( std::map<double, double>::iterator iter = tpbemission.begin(); iter != tpbemission.end(); ++iter)
    {
      parent[ii++]=(*iter).second;
    }
    rgen0 = new CLHEP::RandGeneral(parent,nbins);
    
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
    // The ratio of slow photons to fast photons is related        by the yieldratio
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

    const G4Track * aTrack = aStep.GetTrack();

    G4StepPoint const* pPreStepPoint  = aStep.GetPreStepPoint();
    // unused G4StepPoint const* pPostStepPoint = aStep.GetPostStepPoint();
  
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
    auto const& Visibilities = pvs->GetAllVisibilities(xyz);

    phot::MappedCounts_t ReflVisibilities;

    
    // Store timing information in the object for use in propagation_time method
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
    
      if(!Visibilities && !pvs->UseNhitsModel()){
      }else{
        std::map<int, int> DetectedNum;

        std::map<int, int> ReflDetectedNum;

        for(size_t OpDet=0; OpDet!=NOpChannels; OpDet++)
        {
          G4int DetThisPMT = 0.;
	  if(Visibilities && !pvs->UseNhitsModel()){
	    DetThisPMT = G4int(G4Poisson(Visibilities[OpDet] * Num));
	  }
	  else {
	    TVector3 ScintPoint( xyz[0], xyz[1], xyz[2] );
	    TVector3 OpDetPoint(fOpDetCenter.at(OpDet)[0], fOpDetCenter.at(OpDet)[1], fOpDetCenter.at(OpDet)[2]); 
	    DetThisPMT = VUVHits(Num, ScintPoint, OpDetPoint, foptical_detector_type);
         
	  }
	  
          if(DetThisPMT>0) 
          {
            DetectedNum[OpDet]=DetThisPMT;

            //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
            //   //   it->second<<" " << Num << " " << DetThisPMT;  

            //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
          }
          if(pvs->StoreReflected()) {
	    G4int ReflDetThisPMT = 0;
	    if (!pvs->UseNhitsModel()){
	      ReflDetThisPMT = G4int(G4Poisson(ReflVisibilities[OpDet] * Num));
	    }
	    else {
	      TVector3 ScintPoint( xyz[0], xyz[1], xyz[2] );
	      TVector3 OpDetPoint(fOpDetCenter.at(OpDet)[0], fOpDetCenter.at(OpDet)[1], fOpDetCenter.at(OpDet)[2]); 
	      ReflDetThisPMT = VISHits(Num, ScintPoint, OpDetPoint, foptical_detector_type);
	    }	
            
            if(ReflDetThisPMT>0)
            {
	      ReflDetectedNum[OpDet]=ReflDetThisPMT;
	   
            }
          }

        }

        
        // Now we run through each PMT figuring out num of detected photons
        for (int Reflected = 0; Reflected <= 1; Reflected++) {
          // Only do the reflected loop if we have reflected visibilities
          if (Reflected && !pvs->StoreReflected())
            continue;
          
          std::map<int,int>::const_iterator itstart;
          std::map<int,int>::const_iterator itend;
          if (Reflected) {
            itstart = ReflDetectedNum.begin();
            itend   = ReflDetectedNum.end();
          }
          else{
            itstart = DetectedNum.begin();
            itend   = DetectedNum.end();
          }
          
          for(std::map<int,int>::const_iterator itdetphot=itstart; itdetphot!=itend; ++itdetphot)
          {
            int OpChannel = itdetphot->first;
            int NPhotons  = itdetphot->second;

            // Set up the OpDetBTR information
            sim::OpDetBacktrackerRecord tmpOpDetBTRecord(OpChannel);
            int thisG4TrackID = ParticleListAction::GetCurrentTrackID();
            double xyzPos[3];
            average_position(aStep, xyzPos);
            double Edeposited  = larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit()/CLHEP::MeV;

            // Get the transport time distribution
	    std::vector<double> arrival_time_dist = propagation_time(x0, OpChannel, NPhotons, Reflected);
            // Loop through the photons
            for (G4int i = 0; i < NPhotons; ++i)
            {
	      //std::cout<<"VUV time correction: "<<arrival_time_dist[i]<<std::endl;
              G4double Time = t0
                + scint_time(aStep, ScintillationTime, ScintillationRiseTime) 
		+ arrival_time_dist[i]*CLHEP::ns;

              // Always store the BTR
              tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, Time, 1, xyzPos, Edeposited);

              // Store as lite photon or as OnePhoton
              if(lgp->UseLitePhotons())
              {
                fst->AddLitePhoton(OpChannel, static_cast<int>(Time), 1, Reflected);
              }
              else {
                // The sim photon in this case stores its production point and time
                TVector3 PhotonPosition( x0[0], x0[1], x0[2] );

                float PhotonEnergy = 0;
                if (Reflected)  PhotonEnergy = reemission_energy()*CLHEP::eV;
                else            PhotonEnergy = 9.7*CLHEP::eV;
            
                // Make a photon object for the collection
                sim::OnePhoton PhotToAdd;
                PhotToAdd.InitialPosition  = PhotonPosition;
                PhotToAdd.Energy           = PhotonEnergy;
                PhotToAdd.Time             = Time;
                PhotToAdd.SetInSD          = false;
                PhotToAdd.MotherTrackID    = tracknumber;

                fst->AddPhoton(OpChannel, std::move(PhotToAdd), Reflected);
              }
            }
            fst->AddOpDetBacktrackerRecord(tmpOpDetBTRecord, Reflected);
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

  G4double OpFastScintillation::scint_time(const G4Step& aStep,
                                           G4double ScintillationTime,
                                           G4double ScintillationRiseTime) const
  {
    G4StepPoint const* pPreStepPoint  = aStep.GetPreStepPoint();
    G4StepPoint const* pPostStepPoint = aStep.GetPostStepPoint();
    G4double avgVelocity = (pPreStepPoint->GetVelocity() + pPostStepPoint->GetVelocity())/2.;
    
    G4double deltaTime = aStep.GetStepLength() / avgVelocity;

    if (ScintillationRiseTime==0.0) {
      deltaTime = deltaTime -
        ScintillationTime * std::log( G4UniformRand() );
    } else {
      deltaTime = deltaTime +
        sample_time(ScintillationRiseTime, ScintillationTime);
    }
    return deltaTime;
  }


  std::vector<double> OpFastScintillation::propagation_time(G4ThreeVector x0, int OpChannel, int NPhotons, bool Reflected) //const
  {

    static art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    // Initialize vector of the right length with all 0's
    std::vector<double> arrival_time_dist(NPhotons, 0);


    if (pvs->IncludeParPropTime() && pvs->IncludePropTime()) {
      throw cet::exception("OpFastScintillation") << "Cannot have both propagation time models simultaneously.";
    }

    else if (pvs->IncludeParPropTime() && !(ParPropTimeTF1  && (ParPropTimeTF1[OpChannel].GetNdim()==1)) )
    {
      //Warning: TF1::GetNdim()==1 will tell us if the TF1 is really defined or it is the default one.
      //This will fix a segfault when using timing and interpolation.
      G4cout << "WARNING: Requested parameterized timing, but no function found. Not applying propagation time." << G4endl;
    }
    
    else if (pvs->IncludeParPropTime()) {
      if (Reflected)
        throw cet::exception("OpFastScintillation") << "No parameterized propagation time for reflected light";
        
      for (int i = 0; i < NPhotons; i++) {
        arrival_time_dist[i] = ParPropTimeTF1[OpChannel].GetRandom();
      }
    }

    else if (pvs->IncludePropTime()) {
      // Get VUV photons arrival time distribution from the parametrization 
      G4ThreeVector OpDetPoint(fOpDetCenter.at(OpChannel)[0]*CLHEP::cm,fOpDetCenter.at(OpChannel)[1]*CLHEP::cm,fOpDetCenter.at(OpChannel)[2]*CLHEP::cm);
      
      if (!Reflected) {
        double distance_in_cm = (x0 - OpDetPoint).mag()/CLHEP::cm; // this must be in CENTIMETERS! 
        arrival_time_dist = getVUVTime(distance_in_cm, NPhotons); // in ns

      }
      else {
	TVector3 ScintPoint( x0[0]/CLHEP::cm, x0[1]/CLHEP::cm, x0[2]/CLHEP::cm ); // in cm
	TVector3 OpDetPoint_tv3(fOpDetCenter.at(OpChannel)[0], fOpDetCenter.at(OpChannel)[1], fOpDetCenter.at(OpChannel)[2]); // in cm
        arrival_time_dist = getVISTime(ScintPoint, OpDetPoint_tv3, NPhotons); // in ns
      }
    }
  
    return arrival_time_dist;
  }
 


  
  G4double OpFastScintillation::sample_time(G4double tau1, G4double tau2) const
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


  double OpFastScintillation::reemission_energy() const
  {
    return rgen0->fire()*((*(--tpbemission.end())).first-(*tpbemission.begin()).first)+(*tpbemission.begin()).first;
  }


  void OpFastScintillation::average_position(G4Step const& aStep, double *xyzPos) const
  {
    CLHEP::Hep3Vector prePoint  = (aStep.GetPreStepPoint())->GetPosition();
    CLHEP::Hep3Vector postPoint = (aStep.GetPostStepPoint())->GetPosition();
    xyzPos[0] = ( ( (prePoint.getX() + postPoint.getX() ) / 2.0) / CLHEP::cm );
    xyzPos[1] = ( ( (prePoint.getY() + postPoint.getY() ) / 2.0) / CLHEP::cm );
    xyzPos[2] = ( ( (prePoint.getZ() + postPoint.getZ() ) / 2.0) / CLHEP::cm );
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
    //G4cout<<"BEAMAUS timing distribution hecha"<<G4endl;
    return arrival_time_distrb;
  }

// Parametrization of the Visible light timing (result from direct transport + Rayleigh scattering ONLY) 
// using a landau + exponential function. The function below returns the arrival time distribution given the 
// time of the first visible photon in the PMT. The light generated has been reflected by the cathode ONLY.
  std::vector<double> OpFastScintillation::GetVisibleTimeOnlyCathode(double t0, int number_photons) {
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
    //G4cout<<"Timing distribution BEAMAUS!"<<G4endl;
    return arrival_time_distrb;
  }

  // New Parametrization code 
  // parameterisation generation function
  void OpFastScintillation::generateparam(int index) {
    // get distance 
    double distance_in_cm = (index * fstep_size) + 25;
    
    // time range
    const double signal_t_range = 5000.;
   
    // parameterisation TF1    
    TF1* fVUVTiming;
      
    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm/fvuv_vgroup_mean;
    double t_direct_min = distance_in_cm/fvuv_vgroup_max;
      
    // Defining the model function(s) describing the photon transportation timing vs distance 
    // Getting the landau parameters from the time parametrization
    double* pars_landau = interpolate(fparameters[0], fparameters[2], fparameters[3], fparameters[1], distance_in_cm, true);
    // Deciding which time model to use (depends on the distance)
    // defining useful times for the VUV arrival time shapes
    if(distance_in_cm >= finflexion_point_distance) {
      double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
      // Set model: Landau 
      fVUVTiming =  new TF1("fVUVTiming",model_far,0,signal_t_range,4);
      fVUVTiming->SetParameters(pars_far);
    }
    else {
      // Set model: Landau + Exponential 
      fVUVTiming =  new TF1("fVUVTiming",model_close,0,signal_t_range,7); 
      // Exponential parameters
      double pars_expo[2];   
      // Getting the exponential parameters from the time parametrization
      pars_expo[1] = interpolate(fparameters[4], fparameters[5], distance_in_cm, true);
      //For simplicity, not considering the small dependency with the offset angle in pars_expo[0]
      //Using the value for the [30,60deg] range. fparameters[6] and fparameters[8] are the values
      //for [0,30deg] range and [60,90deg] range respectively
      pars_expo[0] = fparameters[7].at(0) + fparameters[7].at(1)*distance_in_cm;
      pars_expo[0] *= pars_landau[2];
      pars_expo[0] = log(pars_expo[0]);
      // this is to find the intersection point between the two functions:
      TF1* fint = new TF1("fint",finter_d,pars_landau[0],4*t_direct_mean,5);
      double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
      fint->SetParameters(parsInt);
      double t_int = fint->GetMinimumX();
      double minVal = fint->Eval(t_int);
      // the functions must intersect - output warning if they don't
      if(minVal>0.015) {
	std::cout<<"WARNING: Parametrization of VUV light discontinuous for distance = " << distance_in_cm << std::endl;
	std::cout<<"WARNING: This shouldn't be happening " << std::endl;
      }
      delete fint;   
      double parsfinal[7] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1], t_direct_min};
      fVUVTiming->SetParameters(parsfinal);    
    }

    // set the number of points used to sample parameterisation
    // for shorter distances, peak is sharper so more sensitive sampling required - values could be optimised, but since these are only generate once difference is not significant
    int f_sampling;
    if (distance_in_cm < 50) { f_sampling = 10000; }
    else if (distance_in_cm < 100){ f_sampling = 5000; }
    else { f_sampling = 1000; }
    fVUVTiming->SetNpx(f_sampling);    

    // calculate max and min distance relevant to sample parameterisation
    // max 
    const int nq_max=1;
    double xq_max[nq_max];
    double yq_max[nq_max];
    xq_max[0] = 0.99;   // include 99%, 95% cuts out a lot of tail and time difference is negligible extending this
    fVUVTiming->GetQuantiles(nq_max,yq_max,xq_max);
    double max = yq_max[0];
    // min 
    double min = t_direct_min;

    // generate the sampling
    // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
    // all subsequent calls check if it has been generated previously and are ~100+ times quicker
    //double arrival_time = fVUVTiming->GetRandom(min,max);
    // add timing to the vector of timings and range to vectors of ranges
    VUV_timing[index] = *fVUVTiming;
    VUV_max[index] = max;
    VUV_min[index] = min;

    delete fVUVTiming;
  }

  // VUV arrival times calculation function
  std::vector<double> OpFastScintillation::getVUVTime(double distance, int number_photons) {
    
    // pre-allocate memory
    std::vector<double> arrival_time_distrb;
    arrival_time_distrb.clear();
    arrival_time_distrb.reserve(number_photons);
    
    // distance < 25cm
    if (distance < 25) {
      // times are fixed shift i.e. direct path only
      double t_prop_correction = distance/fvuv_vgroup_mean;
      for (int i = 0; i < number_photons; i++){
	arrival_time_distrb.push_back(t_prop_correction);
      }
    }
    // distance >= 25cm
    else {
      // determine nearest parameterisation in discretisation
      int index = std::round((distance - 25) / fstep_size);
      // check whether required parameterisation has been generated, generating if not
      if (VUV_timing[index].GetNdim() == 0) {
	generateparam(index);
      }
      // randomly sample parameterisation for each photon
      for (int i = 0; i < number_photons; i++){
	arrival_time_distrb.push_back(VUV_timing[index].GetRandom(VUV_min[index],VUV_max[index])); 
      }  
    }
    return arrival_time_distrb;

  }

  // VIS arrival times calculation functions
  std::vector<double> OpFastScintillation::getVISTime(TVector3 ScintPoint, TVector3 OpDetPoint, int number_photons) {
    // *************************************************************************************************
    //     Calculation of earliest arrival times and corresponding unsmeared distribution
    // *************************************************************************************************
   
    // set plane_depth for correct TPC:
    double plane_depth;
    if (ScintPoint[0] < 0) {
      plane_depth = -fplane_depth;
    }
    else {
      plane_depth = fplane_depth;
    }
 
    // calculate point of reflection for shortest path accounting for difference in refractive indicies    
    // vectors for storing results
    TVector3 image(0,0,0);
    TVector3 bounce_point(0,0,0);
                                 
    // distance to wall    
    TVector3 v_to_wall(plane_depth-ScintPoint[0],0,0);
    
    // hotspot is point on wall where TPB is activated most intensely by the scintillation
    TVector3 hotspot(plane_depth,ScintPoint[1],ScintPoint[2]);
                                                     
    // define "image" by reflecting over plane
    image = hotspot + v_to_wall*(fn_LAr_vis/fn_LAr_vuv);
                                                                 
    // find point of intersection with plane j of ray from the PMT to the image
    TVector3 tempvec = (OpDetPoint-image).Unit();
    double tempnorm= ((image-hotspot).Mag())/std::abs(tempvec[0]);
    bounce_point = image + tempvec*tempnorm;
    
    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point-ScintPoint).Mag();
    double Visdist = (OpDetPoint-bounce_point).Mag();

    // calculate times taken by each part
    std::vector<double> VUVTimes  = getVUVTime(VUVdist, number_photons);
    std::vector<double> ReflTimes(number_photons,Visdist/fvis_vmean);
   
    // sum parts to get total transport times times                    
    std::vector<double> transport_time_vis(number_photons,0);
    for (int i=0; i<number_photons; i++) {
      transport_time_vis[i] = VUVTimes[i] + ReflTimes[i];
    }
   
    // *************************************************************************************************
    //      Smearing of arrival time distribution
    // *************************************************************************************************
           
    // calculate fastest time possible
    // vis part
    double vis_time = Visdist/fvis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < 25){
      vuv_time = VUVdist/fvuv_vgroup_mean;
    }
    else {
      // find index of required parameterisation
      int index = std::round((VUVdist - 25) / fstep_size);
      // find shortest time
      vuv_time = VUV_min[index];
    }
    // sum
    double fastest_time = vis_time + vuv_time;
    
    // calculate angle alpha between scintillation point and reflection point
    double cosine_alpha = sqrt(pow(ScintPoint[0] - bounce_point[0],2)) / VUVdist;
    double alpha = acos(cosine_alpha)*180./CLHEP::pi;
    
    // determine smearing parameters using interpolation of generated points: 
    // 1). tau = exponential smearing factor, varies with distance and angle
    // 2). cutoff = largest smeared time allowed, preventing excessively large times caused by exponential
    // distance to cathode
    double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
    // angular bin
    unsigned int alpha_bin = alpha / 10;
    if (alpha_bin >= ftau_pars.size()) {
      alpha_bin = ftau_pars.size() - 1;      // default to the largest available bin if alpha larger than parameterised region; i.e. last bin effectively [last bin start value, 90] deg bin
    }
   // cut-off and tau
    double cutoff = interpolate( fdistances_refl, fcut_off_pars[alpha_bin], distance_cathode_plane, true );
    double tau = interpolate( fdistances_refl, ftau_pars[alpha_bin], distance_cathode_plane, true );
   
   // fail-safe if tau extrapolate goes wrong, drops below zero since last distance close to zero [did not occur in testing, but possible]
    if (tau < 0){
      tau = 0;
    } 
	
    // apply smearing:
    for (int i=0; i < number_photons; i++){
      double arrival_time = transport_time_vis[i];
      double arrival_time_smeared;
      // if time is already greater than cutoff or minimum smeared time would be greater than cutoff, do not apply smearing
      if (arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(1.0))-1) >= cutoff) {
        arrival_time_smeared = arrival_time;
      }
      // otherwise smear
      else {
        int counter = 0;
        // loop until time generated is within cutoff limit
        // most are within single attempt, very few take more than two
        do {
          // don't attempt smearings too many times for cases near cutoff (very few cases, not smearing these makes negigible difference)
          if (counter >= 10){
            arrival_time_smeared = arrival_time; // don't smear
            break;
          }
          else {
            // generate random number in appropriate range            
            double x = gRandom->Uniform(0.5,1.0);
      	    // apply the exponential smearing
      	    arrival_time_smeared = arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(x))-1);
	  }
	  // increment counter
	  counter++;
	}  while (arrival_time_smeared > cutoff);
      }
      transport_time_vis[i] = arrival_time_smeared;
    }
    return transport_time_vis;
  }

  // VUV semi-analytic hits calculation
  int OpFastScintillation::VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type) {
    // check optical channel is in same TPC as scintillation light, if not return 0 hits
    // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero, otherwise return 0 hits 
    if (((ScintPoint[0] < 0) != (OpDetPoint[0] < 0)) && std::abs(OpDetPoint[0]) > 10){	
      return 0;	
    }    
 
    // distance and angle between ScintPoint and OpDetPoint
    double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
    double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
    double theta = acos(cosine)*180./CLHEP::pi;
    
    // calculate solid angle:
    double solid_angle = 0;
    double d, h;
    // Arapucas
    if (optical_detector_type == 0) {
      // set Arapuca geometry struct for solid angle function
      acc detPoint; 
      detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
      detPoint.w = fydimension; detPoint.h = fzdimension; // width and height in cm of arapuca active window

      // get scintillation point coordinates relative to arapuca window centre
      TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;  

      // calculate solid angle
      solid_angle = Rectangle_SolidAngle(detPoint, ScintPoint_rel);

    }
    // PMTs
    else if (optical_detector_type == 1) {
      // offset in z-y plane
      d = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
      // drift distance (in x)
      h =  sqrt(pow(ScintPoint[0] - OpDetPoint[0],2));
      // Solid angle of a disk
      solid_angle = Disk_SolidAngle(d, h, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" <<std:: endl;
    }  

    // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
    double hits_geo = exp(-1.*distance/fL_abs_vuv) * (solid_angle / (4*CLHEP::pi)) * Nphotons_created;

    // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
    // offset angle bin
    int j = (theta/fdelta_angulo);
    double hits_rec = gRandom->Poisson( GHvuv[j]->Eval(distance)*hits_geo/cosine );

    // round to integer value, cannot have non-integer number of hits
    int hits_vuv = std::round(hits_rec);

    return hits_vuv;
  }

  // VIS hits semi-analytic model calculation 
  int OpFastScintillation::VISHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type) {
     // check optical channel is in same TPC as scintillation light, if not return 0 hits
     // temporary method working for SBND, DUNE 1x2x6; to be replaced to work in full DUNE geometry
     // check x coordinate has same sign or is close to zero, otherwise return 0 hits
     if (((ScintPoint[0] < 0) != (OpDetPoint[0] < 0)) && std::abs(OpDetPoint[0]) > 10){	
       return 0;	
     }   

    // set plane_depth for correct TPC:
     double plane_depth;
     if (ScintPoint[0] < 0) {
       plane_depth = -fplane_depth;
     }
     else {
       plane_depth = fplane_depth;
     }
 
     // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:
     
     // set cathode plane struct for solid angle function
     acc cathode_plane; 
     cathode_plane.ax = plane_depth; cathode_plane.ay = fcathode_centre[1]; cathode_plane.az = fcathode_centre[2];       	// centre coordinates of cathode plane
     cathode_plane.w = fcathode_ydimension; cathode_plane.h = fcathode_zdimension;                        				// width and height in cm
     
     // get scintpoint coords relative to centre of cathode plane
     TVector3 cathodeCentrePoint(plane_depth,fcathode_centre[1],fcathode_centre[2]);
     TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint; 
     
     // calculate solid angle of cathode from the scintillation point
     double solid_angle_cathode = Rectangle_SolidAngle(cathode_plane, ScintPoint_relative);
     
     // calculate distance and angle between ScintPoint and hotspot
     // vast majority of hits in hotspot region directly infront of scintpoint,therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
     double distance_cathode = std::abs(plane_depth - ScintPoint[0]);
     double cosine_cathode = 1;
     double theta_cathode = 0;
     
     // calculate hits on cathode plane via geometric acceptance
     double cathode_hits_geo = exp(-1.*distance_cathode/fL_abs_vuv) * (solid_angle_cathode / (4.*CLHEP::pi)) * Nphotons_created;
     
     // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
     // offset angle bin
     int j = (theta_cathode/fdelta_angulo);
     double cathode_hits_rec = GHvuv[j]->Eval(distance_cathode)*cathode_hits_geo/cosine_cathode;
      
    
     // 2). calculate number of these hits which reach the optical detector from the hotspot via solid angle 

     // hotspot coordinates         
     TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);
     
     // get hotspot coordinates relative to detpoint
     TVector3 emission_relative = hotspot - OpDetPoint;
   
     // calculate solid angle of optical channel
     double solid_angle_detector = 0;
     // rectangular aperture
     if (optical_detector_type == 0) {
      	// set rectangular aperture geometry struct for solid angle function
     	acc detPoint; 
     	detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  	// centre coordinates of optical detector
     	detPoint.w = fydimension; detPoint.h = fzdimension;	 						// width and height in cm of optical detector active window [rectangular aperture]
	// calculate solid angle
	solid_angle_detector = Rectangle_SolidAngle(detPoint, emission_relative);
     }
     // disk aperture    
     else if (optical_detector_type == 1) {
	// offset in z-y plane
      	double d = sqrt(pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
      	// drift distance (in x)
      	double h =  sqrt(pow(hotspot[0] - OpDetPoint[0],2)); 
	// calculate solid angle
	solid_angle_detector = Disk_SolidAngle(d, h, fradius);
     }
     else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" <<std::endl;
     }
    
     // calculate number of hits via geometeric acceptance
     double hits_geo = (solid_angle_detector / (2.*CLHEP::pi)) * cathode_hits_rec;	// 2*pi rather than 4*pi due to presence of reflective foils (vm2000)
  
     // calculate distances and angles for application of corrections
     // distance to hotspot
     double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
     // distance from hotspot to optical detector
     double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
     //  angle between hotspot and optical detector
     double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
     double theta_vis = acos(cosine_vis)*180./CLHEP::pi; 
     
     // apply correction 
     int k = (theta_vis/fdelta_angulo);
     double hits_rec = gRandom->Poisson(VIS_pol[k]->Eval(distance_vuv)*hits_geo/cosine_vis);
    
     // round final result, number of hits integer
     int hits_vis = std::round(hits_rec);
                                  
     return hits_vis; 
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

  double model_close(double *x, double *par)
  {
    // par0 = joining point
    // par1 = Landau MPV
    // par2 = Landau width
    // par3 = normalization
    // par4 = Expo cte
    // par5 = Expo tau
    // par6 = t_min
  
    double y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
    double y2 = TMath::Exp(par[4]+x[0]*par[5]);
    if(x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
    if(x[0] < par[0]) y2 = 0.;
  
    return (y1 + y2);
  }

  double model_far(double *x, double *par)
  {
    // par1 = Landau MPV
    // par2 = Landau width
    // par3 = normalization
    // par0 = t_min

    double y = par[3]*TMath::Landau(x[0],par[1],par[2]);
    if(x[0] <= par[0]) y = 0.;
  
    return y;
  }
  double GaisserHillas(double *x,double *par) {
    //This is the Gaisser-Hillas function
    double X_mu_0=par[3];
    double Normalization=par[0];
    double Diff=par[1]-X_mu_0;
    double Term=pow((*x-X_mu_0)/Diff,Diff/par[2]);
    double Exponential=TMath::Exp((par[1]-*x)/par[2]);
  
    return (Normalization*Term*Exponential);
  }

  //======================================================================

  //   Returns interpolated value at x from parallel arrays ( xData, yData )
  //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
  //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)

  double interpolate( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
  {
    int size = xData.size();
    int i = 0;                                          // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                         // special case: beyond right end
      {
	i = size - 2;
      }
    else
      {
	while ( x > xData[i+1] ) i++;
      }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
    if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
      {
	if ( x < xL ) yR = yL;
	if ( x > xR ) yL = yR;
      }
    double dydx = ( yR - yL ) / ( xR - xL );            // gradient
    return yL + dydx * ( x - xL );                      // linear interpolation
  }

  double* interpolate( std::vector<double> &xData, std::vector<double> &yData1, std::vector<double> &yData2,
		       std::vector<double> &yData3, double x, bool extrapolate)
  {
    int size = xData.size();
    int i = 0;                                          // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                         // special case: beyond right end
      {
	i = size - 2;
      }
    else
      {
	while ( x > xData[i+1] ) i++;
      }
    double xL = xData[i], xR = xData[i+1];// points on either side (unless beyond ends)
    double yL1 = yData1[i], yR1 = yData1[i+1], yL2 = yData2[i], yR2 = yData2[i+1], yL3 = yData3[i], yR3 = yData3[i+1]; 
  
    if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
      {
	if ( x < xL ) {yR1 = yL1; yR2 = yL2; yR3 = yL3;}
	if ( x > xR ) {yL1 = yR1; yL2 = yR2; yL3 = yR3;}
      }
    double dydx1 = ( yR1 - yL1 ) / ( xR - xL );            // gradient
    double dydx2 = ( yR2 - yL2 ) / ( xR - xL );
    double dydx3 = ( yR3 - yL3 ) / ( xR - xL );

    double *yy = new double[3]; 
    yy[0] = yL1 + dydx1 * ( x - xL );// linear interpolations
    yy[1] = yL2 + dydx2 * ( x - xL );
    yy[2] = yL3 + dydx3 * ( x - xL );
  
    return yy;                      
  }

  // solid angle of circular aperture
  double Disk_SolidAngle(double* x, double *p) {
    const double d = x[0];
    const double h = x[1];
    const double b = p[0];
    if(b <= 0. || d < 0. || h <= 0.) return 0.; 
    const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
    if(d == 0) {
      return 2.*TMath::Pi()*(1.-aa);
    }
    const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
    const double cc = 4*b*d/((b+d)*(b+d));
 
    if(TMath::Abs(boost::math::ellint_1(bb) - bb) < 1e-10 && TMath::Abs(boost::math::ellint_3(cc,bb) - cc) <1e-10) {
      throw(std::runtime_error("Problem loading ELLIPTIC INTEGRALS running Disk_SolidAngle!"));
    }
    if(d < b) {
     return 2.*TMath::Pi() - 2.*aa*(boost::math::ellint_1(bb) + TMath::Sqrt(1.-cc)*boost::math::ellint_3(bb,cc));
    }
    if(d == b) {
      return TMath::Pi() - 2.*aa*boost::math::ellint_1(bb);
    }
    if(d > b) {
      return 2.*aa*(TMath::Sqrt(1.-cc)*boost::math::ellint_3(bb,cc) - boost::math::ellint_1(bb));
    }
    return 0.;
  }

  double Disk_SolidAngle(double d, double h, double b) {
    double x[2] = { d, h };
    double p[1] = { b };

    return Disk_SolidAngle(x,p);
    }
 
  // solid angle of rectanglular aperture 

  double Rectangle_SolidAngle(double a, double b, double d){

    double aa = a/(2.0*d);
    double bb = b/(2.0*d);
    double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
    return 4*std::acos(std::sqrt(aux));

  }
  
  double Rectangle_SolidAngle(acc& out, TVector3 v){

    //v is the position of the track segment with respect to 
    //the center position of the arapuca window 
 
    // arapuca plane fixed in x direction

    if( v.Y()==0.0 && v.Z()==0.0){
      return Rectangle_SolidAngle(out.w,out.h,v.X());
    }
  
    if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
      double A, B, a, b, d;
      A = std::abs(v.Y())-out.w/2.0;
      B = std::abs(v.Z())-out.h/2.0;
      a = out.w;
      b = out.h;
      d = std::abs(v.X());
      double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
      return to_return;
    }
  
    if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
      double A, B, a, b, d;
      A = -std::abs(v.Y())+out.w/2.0;
      B = -std::abs(v.Z())+out.h/2.0;
      a = out.w;
      b = out.h;
      d = std::abs(v.X());
      double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
      return to_return;
    }

    if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
      double A, B, a, b, d;
      A = std::abs(v.Y())-out.w/2.0;
      B = -std::abs(v.Z())+out.h/2.0;
      a = out.w;
      b = out.h;
      d = std::abs(v.X());
      double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
      return to_return;
    }

    if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
      double A, B, a, b, d;
      A = -std::abs(v.Y())+out.w/2.0;
      B = std::abs(v.Z())-out.h/2.0;
      a = out.w;
      b = out.h;
      d = std::abs(v.X());
      double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
      return to_return;
    }
    // error message if none of these cases, i.e. something has gone wrong!
    std::cout << "Warning: invalid solid angle call." << std::endl;
    return 0.0;
    }

}
