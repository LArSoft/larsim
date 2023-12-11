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
// Updated:     
//              2022-11-11 By Jiaoyang Li <jiaoyang.li@ed.ac.uk>
//              > Implemented the semi-analytical model to optical fast
//              simulation. https://doi.org/10.1140/epjc/s10052-021-09119-3
//              2010-10-20 Allow the scintillation yield to be a function
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
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"         // geo::vect::fillCoords()

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// support libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h" // To-do: not sure whether is applicable in that stage. 

#include "TRandom3.h"
#include "TMath.h"
#include "boost/algorithm/string.hpp"

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

// Define a new policy *not* internally promoting RealType to double:
typedef boost::math::policies::policy<
  // boost::math::policies::digits10<8>,
  boost::math::policies::promote_double<false>>
  noLDoublePromote;
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
    , fActiveVolumes{extractActiveVolumes(*(lar::providerFrom<geo::Geometry>()))}
    , bPropagate(!(art::ServiceHandle<sim::LArG4Parameters>()->NoPhotonPropagation()))
    , fUseNhitsModel(bPropagate ? art::ServiceHandle<phot::PhotonVisibilityService const>().get()->UseNhitsModel() : false)
    , fOnlyActiveVolume(usesSemiAnalyticModel())
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
      geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

      if (usesSemiAnalyticModel()) {  
        if (geom.Ncryostats() > 1U) {
          if (fOnlyOneCryostat) {
            mf::LogWarning("OpFastScintillation")
              << std::string(80, '=') << "\nA detector with " << geom.Ncryostats()
              << " cryostats is configured"
              << " , and semi-analytic model is requested for scintillation photon propagation."
              << " THIS CONFIGURATION IS NOT SUPPORTED and it is open to bugs"
              << " (e.g. scintillation may be detected only in cryostat #0)."
              << "\nThis would be normally a fatal error, but it has been forcibly overridden."
              << "\n"
              << std::string(80, '=');
          }
          else {
            throw art::Exception(art::errors::Configuration)
              << "Photon propagation via semi-analytic model is not supported yet"
              << " on detectors with more than one cryostat.";
          }
        }
      } // if

      // Loading the position of each optical channel, neccessary for the parametrizatiuons of Nhits and prop-time
      
      geo::Point_t const Cathode_centre{geom.TPC(0, 0).GetCathodeCenter().X(),
                                        fActiveVolumes[0].CenterY(),
                                        fActiveVolumes[0].CenterZ()};
      mf::LogTrace("OpFastScintillation") << "Cathode_centre: " << Cathode_centre << " cm";

      for (size_t const i : util::counter(pvs->NOpChannels())) {
        geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpDet(i);
        fOpDetCenter.push_back(opDet.GetCenter());
          
        fOpDetType.push_back(2); // Disk PMTs
        //std::cout<<"Radio: "<<geom.OpDetGeoFromOpDet(i).RMax()<<std::endl;
        fOpDetLength.push_back(-1);
        fOpDetHeight.push_back(-1);

      }

      if(pvs->IncludePropTime()) {
        pvs->SetDirectLightPropFunctions(functions_vuv, fd_break, fd_max, ftf1_sampling_factor);
        pvs->SetReflectedCOLightPropFunctions(functions_vis, ft0_max, ft0_break_point);
      }

      if (usesSemiAnalyticModel()) {
        mf::LogVerbatim("OpFastScintillation")
          << "OpFastScintillation: using semi-analytic model for number of hits";

        // LAr absorption length in cm
        std::map<double, double> abs_length_spectrum =
          lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
        std::vector<double> x_v, y_v;
        for (auto elem : abs_length_spectrum) {
          x_v.push_back(elem.first);
          y_v.push_back(elem.second);
        }
        fL_abs_vuv =
          std::round(interpolate(x_v, y_v, 9.7, false)); // 9.7 eV: peak of VUV emission spectrum

        // Load Gaisser-Hillas corrections for VUV semi-analytic hits
        std::cout << "Loading the GH corrections" << std::endl;
        pvs->LoadVUVSemiAnalyticProperties(fIsFlatPDCorr, fdelta_angulo_vuv, fradius);
        if (fIsFlatPDCorr) {
          pvs->LoadGHFlat(fGHvuvpars_flat, fborder_corr_angulo_flat, fborder_corr_flat);
        }else {
            throw cet::exception("OpFastScintillation")
            << "Both isFlatPDCorr and isDomePDCorr parameters are false, at least one type of "
               "parameterisation is required for the semi-analytic light simulation."
            << "\n";
        }

        // cathode center coordinates required for corrections
        fcathode_centre = geom.TPC(0, 0).GetCathodeCenter();
        fcathode_centre[1] = fActiveVolumes[0].CenterY();
        fcathode_centre[2] = fActiveVolumes[0].CenterZ(); // to get full cathode dimension rather than just single tpc
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
    geo::Point_t const ScintPoint = {xyz[0], xyz[1], xyz[2]};
    if (fOnlyActiveVolume && !isScintInActiveVolume(ScintPoint)) return 0;

    float const* Visibilities = pvs->GetAllVisibilities(xyz);

    float const* ReflVisibilities = nullptr;

    
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
      if (!Visibilities && !usesSemiAnalyticModel()) continue; // skip if not using photon library or semi-analytical model
     
      std::map<int, int> DetectedNum;

      std::map<int, int> ReflDetectedNum;

      if (Visibilities && !usesSemiAnalyticModel()) { // photon library
        for(size_t OpDet=0; OpDet!=NOpChannels; OpDet++)
        {
          G4int DetThisPMT = G4int(G4Poisson(Visibilities[OpDet] * Num));

          if(DetThisPMT>0) // direct light
          {
            DetectedNum[OpDet]=DetThisPMT; 
            //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
            //   //   it->second<<" " << Num << " " << DetThisPMT;  

            //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
          }
          if(pvs->StoreReflected()) { // reflected light
            G4int ReflDetThisPMT = G4int(G4Poisson(ReflVisibilities[OpDet] * Num));
            if(ReflDetThisPMT>0)
            {
              ReflDetectedNum[OpDet]=ReflDetThisPMT;
            }
          }
        }
      }
      else { // semi-analytical model for direct light.
        detectedDirectHits(DetectedNum, Num, ScintPoint);
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


  std::vector<double> OpFastScintillation::propagation_time(G4ThreeVector x0, int OpChannel, int NPhotons, bool Reflected) const
  {
    static art::ServiceHandle<geo::Geometry> geo;
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
      double OpDetCenter[3];
      geo->OpDetGeoFromOpDet(OpChannel).GetCenter(OpDetCenter);
      G4ThreeVector OpDetPoint(OpDetCenter[0]*CLHEP::cm,OpDetCenter[1]*CLHEP::cm,OpDetCenter[2]*CLHEP::cm);
      if (!Reflected) {
        double distance_in_cm = (x0 - OpDetPoint).mag()/CLHEP::cm; // this must be in CENTIMETERS! 
        arrival_time_dist = GetVUVTime(distance_in_cm, NPhotons);//in ns
      }
      else {
        double t0_vis_lib = ReflT0s[OpChannel];                   
        arrival_time_dist = GetVisibleTimeOnlyCathode(t0_vis_lib, NPhotons);
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
  std::vector<double> OpFastScintillation::GetVUVTime(double distance, int number_photons) const {

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
  std::vector<double> OpFastScintillation::GetVisibleTimeOnlyCathode(double t0, int number_photons) const{
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

  // ---------------------------------------------------------------------------
  bool
  OpFastScintillation::usesSemiAnalyticModel() const
  {
    return fUseNhitsModel;
  } // OpFastScintillation::usesSemiAnalyticModel()

  // ---------------------------------------------------------------------------
  void
  OpFastScintillation::detectedDirectHits(std::map<int, int>& DetectedNum,
                                          const double Num,
                                          geo::Point_t const& ScintPoint) const
  {
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;
    for (size_t const OpDet : util::counter(pvs->NOpChannels())) {
      //if (!isOpDetInSameTPC(ScintPoint, fOpDetCenter.at(OpDet))) continue;

      // set detector struct for solid angle function
      const OpFastScintillation::OpticalDetector op{
        fOpDetHeight.at(OpDet), fOpDetLength.at(OpDet),
        fOpDetCenter.at(OpDet), fOpDetType.at(OpDet)};
      const int DetThis = VUVHits(Num, ScintPoint, op);
      if (DetThis > 0) {
        DetectedNum[OpDet] = DetThis;
        //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
        //   //   it->second<<" " << Num << " " << DetThisPMT;
        //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
      }
    }
  }
  // VUV semi-analytic hits calculation
  int
  OpFastScintillation::VUVHits(const double Nphotons_created,
                               geo::Point_t const& ScintPoint_v,
                               OpticalDetector const& opDet) const
  {
    // the interface has been converted into geo::Point_t, the implementation not yet
    std::array<double, 3U> ScintPoint;
    std::array<double, 3U> OpDetPoint;
    geo::vect::fillCoords(ScintPoint, ScintPoint_v);
    geo::vect::fillCoords(OpDetPoint, opDet.OpDetPoint);

    // distance and angle between ScintPoint and OpDetPoint
    double distance = dist(&ScintPoint[0], &OpDetPoint[0], 3);
    double cosine = std::abs(ScintPoint[0] - OpDetPoint[0]) / distance;
    // double theta = std::acos(cosine) * 180. / CLHEP::pi;
    double theta = fast_acos(cosine) * 180. / CLHEP::pi;

    // calculate solid angle:
    double solid_angle = 0;
    // Arapucas/Bars (rectangle)
    if (opDet.type == 2) {
      // offset in z-y plane
      double d = dist(&ScintPoint[1], &OpDetPoint[1], 2);
      // drift distance (in x)
      double h = std::abs(ScintPoint[0] - OpDetPoint[0]);
      // Solid angle of a disk
      solid_angle = Disk_SolidAngle(d, h, fradius);
    }else {
      std::cout << "Error: Invalid optical detector type. 2 = disk" << std::endl;
    }

    // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
    double hits_geo = std::exp(-1. * distance / fL_abs_vuv) * (solid_angle / (4 * CLHEP::pi)) * Nphotons_created;

    // apply Gaisser-Hillas correction for Rayleigh scattering distance
    // and angular dependence offset angle bin
    const size_t j = (theta / fdelta_angulo_vuv);

    // determine GH parameters, accounting for border effects
    // radial distance from centre of detector (Y-Z)
    double r = dist(ScintPoint, fcathode_centre, 2, 1);
    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0; double s2 = 0; double s3 = 0;
    // flat PDs
    if (opDet.type == 2 && fIsFlatPDCorr){
      pars_ini[0] = fGHvuvpars_flat[0][j];
      pars_ini[1] = fGHvuvpars_flat[1][j];
      pars_ini[2] = fGHvuvpars_flat[2][j];
      pars_ini[3] = fGHvuvpars_flat[3][j];
      s1 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[0], theta, true);
      s2 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[1], theta, true);
      s3 = interpolate( fborder_corr_angulo_flat, fborder_corr_flat[2], theta, true);
    }
    else std::cout << "Error: Invalid optical detector type. 2 = disk. Or corrections for chosen optical detector type missing." << std::endl;

    // add border correction
    pars_ini[0] = pars_ini[0] + s1 * r;
    pars_ini[1] = pars_ini[1] + s2 * r;
    pars_ini[2] = pars_ini[2] + s3 * r;
    pars_ini[3] = pars_ini[3];

    // calculate correction
    double GH_correction = Gaisser_Hillas(distance, pars_ini);

    // calculate number photons
    double hits_rec = GH_correction * hits_geo / cosine;
    int hits_vuv = std::round(G4Poisson(hits_rec));
    return hits_vuv;
  }

  bool
  OpFastScintillation::isScintInActiveVolume(geo::Point_t const& ScintPoint)
  {
    //semi-analytic approach only works in the active volume
    return fActiveVolumes[0].ContainsPosition(ScintPoint);
  }

  G4double
  OpFastScintillation::Gaisser_Hillas(const double x, const double* par) const
  {
    double X_mu_0 = par[3];
    double Normalization = par[0];
    double Diff = par[1] - X_mu_0;
    double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
    double Exponential = std::exp((par[1] - x) / par[2]);
    return (Normalization * Term * Exponential);
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

  //======================================================================
  //   Returns interpolated value at x from parallel arrays ( xData, yData )
  //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
  //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
  double
  OpFastScintillation::interpolate(const std::vector<double>& xData,
                                   const std::vector<double>& yData,
                                   double x,
                                   bool extrapolate,
                                   size_t i) const
  {
    if (i == 0) {
      size_t size = xData.size();
      if (x >= xData[size - 2]) { // special case: beyond right end
        i = size - 2;
      }
      else {
        while (x > xData[i + 1])
          i++;
      }
    }
    double xL = xData[i];
    double xR = xData[i + 1];
    double yL = yData[i];
    double yR = yData[i + 1]; // points on either side (unless beyond ends)
    if (!extrapolate) {       // if beyond ends of array and not extrapolating
      if (x < xL) return yL;
      if (x > xR) return yL;
    }
    const double dydx = (yR - yL) / (xR - xL); // gradient
    return yL + dydx * (x - xL);               // linear interpolation
  }
  void
  OpFastScintillation::interpolate3(std::array<double, 3>& inter,
                                    const std::vector<double>& xData,
                                    const std::vector<double>& yData1,
                                    const std::vector<double>& yData2,
                                    const std::vector<double>& yData3,
                                    double x,
                                    bool extrapolate) const
  {
    size_t size = xData.size();
    size_t i = 0;               // find left end of interval for interpolation
    if (x >= xData[size - 2]) { // special case: beyond right end
      i = size - 2;
    }
    else {
      while (x > xData[i + 1])
        i++;
    }
    double xL = xData[i];
    double xR = xData[i + 1]; // points on either side (unless beyond ends)
    double yL1 = yData1[i];
    double yR1 = yData1[i + 1];
    double yL2 = yData2[i];
    double yR2 = yData2[i + 1];
    double yL3 = yData3[i];
    double yR3 = yData3[i + 1];

    if (!extrapolate) { // if beyond ends of array and not extrapolating
      if (x < xL) {
        inter[0] = yL1;
        inter[1] = yL2;
        inter[2] = yL3;
        return;
      }
      if (x > xR) {
        inter[0] = yL1;
        inter[1] = yL2;
        inter[2] = yL3;
        return;
      }
    }
    const double m = (x - xL) / (xR - xL);
    inter[0] = m * (yR1 - yL1) + yL1;
    inter[1] = m * (yR2 - yL2) + yL2;
    inter[2] = m * (yR3 - yL3) + yL3;
  }

  //......................................................................
  // solid angle of circular aperture
  double OpFastScintillation::Disk_SolidAngle(const double d, const double h, const double b) const
  {
    if (b <= 0. || d < 0. || h <= 0.) return 0.;
    const double leg2 = (b + d) * (b + d);
    const double aa = std::sqrt(h * h / (h * h + leg2));
    if (isApproximatelyZero(d)) { return 2. * CLHEP::pi * (1. - aa); }
    double bb = 2. * std::sqrt(b * d / (h * h + leg2));
    double cc = 4. * b * d / leg2;

    if (isDefinitelyGreaterThan(d, b)) {
      try {
        return 2. * aa *
               (std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()) -
                boost::math::ellint_1(bb, noLDoublePromote()));
      }
      catch (std::domain_error& e) {
        if (isApproximatelyEqual(d, b, 1e-9)) {
          mf::LogWarning("SemiAnalyticalModel")
            << "Elliptic Integral in Disk_SolidAngle() given parameters "
               "outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("SemiAnalyticalModel")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters "
               "outside domain.\n"
            << "\nbb: " << bb << "\ncc: " << cc << "Exception message: " << e.what();
          return 0.;
        }
      }
    }
    if (isDefinitelyLessThan(d, b)) {
      try {
        return 2. * CLHEP::pi -
               2. * aa *
                 (boost::math::ellint_1(bb, noLDoublePromote()) +
                  std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()));
      }
      catch (std::domain_error& e) {
        if (isApproximatelyEqual(d, b, 1e-9)) {
          mf::LogWarning("SemiAnalyticalModel")
            << "Elliptic Integral in Disk_SolidAngle() given parameters "
               "outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("SemiAnalyticalModel")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters "
               "outside domain.\n"
            << "\nbb: " << bb << "\ncc: " << cc << "Exception message: " << e.what();
          return 0.;
        }
      }
    }
    if (isApproximatelyEqual(d, b)) {
      return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
    }
    return 0.;
  }

  // ---------------------------------------------------------------------------
  std::vector<geo::BoxBoundedGeo>
  OpFastScintillation::extractActiveVolumes(geo::GeometryCore const& geom)
  {
    std::vector<geo::BoxBoundedGeo> activeVolumes;
    activeVolumes.reserve(geom.Ncryostats());

    for (geo::CryostatGeo const& cryo : geom.IterateCryostats()) {

      // can't use it default-constructed since it would always include origin
      geo::BoxBoundedGeo box{cryo.TPC(0).ActiveBoundingBox()};

      //for (geo::TPCGeo const& TPC : cryo.IterateTPCs())
      box.ExtendToInclude(cryo.TPC(0).ActiveBoundingBox());

      activeVolumes.push_back(std::move(box));

    } // for cryostats

    return activeVolumes;
  } // OpFastScintillation::extractActiveVolumes()
  // ---------------------------------------------------------------------------
  double fast_acos(double x)
  {
    double negate = double(x < 0);
    x = std::abs(x);
    x -= double(x > 1.0) * (x - 1.0); // <- equivalent to min(1.0,x), but faster
    double ret = -0.0187293;
    ret = ret * x;
    ret = ret + 0.0742610;
    ret = ret * x;
    ret = ret - 0.2121144;
    ret = ret * x;
    ret = ret + 1.5707288;
    ret = ret * std::sqrt(1.0 - x);
    ret = ret - 2. * negate * ret;
    return negate * 3.14159265358979 + ret;
  }
}
