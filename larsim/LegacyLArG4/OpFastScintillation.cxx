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

#include "Geant4/globals.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4Poisson.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4Alpha.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4Electron.hh"
#include "Geant4/G4ExceptionSeverity.hh"
#include "Geant4/G4Gamma.hh"
#include "Geant4/G4KaonMinus.hh"
#include "Geant4/G4KaonPlus.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4MaterialPropertyVector.hh"
#include "Geant4/G4MaterialTable.hh"
#include "Geant4/G4MuonMinus.hh"
#include "Geant4/G4MuonPlus.hh"
#include "Geant4/G4ParticleChange.hh"
#include "Geant4/G4PhysicsVector.hh"
#include "Geant4/G4PionMinus.hh"
#include "Geant4/G4PionPlus.hh"
#include "Geant4/G4Proton.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4ios.hh"

#include "larsim/LegacyLArG4/ParticleListAction.h"
#include "larsim/LegacyLArG4/IonizationAndScintillation.h"
#include "larsim/LegacyLArG4/OpFastScintillation.hh"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/LegacyLArG4/OpDetPhotonTable.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// support libraries
#include "cetlib_except/exception.h"

#include "TRandom3.h"
#include "TMath.h"
#include <cmath>
#include <limits>

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

// Define a new policy *not* internally promoting RealType to double:
typedef boost::math::policies::policy<
  // boost::math::policies::digits10<8>,
  boost::math::policies::promote_double<false>
  > noLDoublePromote;

namespace larg4 {

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
    , bPropagate(!(art::ServiceHandle<sim::LArG4Parameters const>()->NoPhotonPropagation()))
  {
    SetProcessSubType(25); // TODO: unhardcode
    fTrackSecondariesFirst = false;
    fFiniteRiseTime = false;
    YieldFactor = 1.0;
    ExcitationRatio = 1.0;

    const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();

    scintillationByParticleType = larp->ScintByParticleType();

    theFastIntegralTable = NULL;
    theSlowIntegralTable = NULL;

    if (verboseLevel > 0) {
      G4cout << GetProcessName() << " is created " << G4endl;
    }

    BuildThePhysicsTable();
    emSaturation = NULL;

    if (bPropagate) {
      art::ServiceHandle<phot::PhotonVisibilityService const> pvs;

      // Loading the position of each optical channel, neccessary for the parametrizatiuons of Nhits and prop-time
      static art::ServiceHandle<geo::Geometry const> geo;

      // Find boundary of active volume
      fminx = 1e9;
      fmaxx = -1e9;
      fminy = 1e9;
      fmaxy = -1e9;
      fminz = 1e9;
      fmaxz = -1e9;
      for (size_t i = 0; i < geo->NTPC(); ++i) {
        const geo::TPCGeo &tpc = geo->TPC(i);
        if (fminx > tpc.MinX()) fminx = tpc.MinX();
        if (fmaxx < tpc.MaxX()) fmaxx = tpc.MaxX();
        if (fminy > tpc.MinY()) fminy = tpc.MinY();
        if (fmaxy < tpc.MaxY()) fmaxy = tpc.MaxY();
        if (fminz > tpc.MinZ()) fminz = tpc.MinZ();
        if (fmaxz < tpc.MaxZ()) fmaxz = tpc.MaxZ();
      }
      std::cout << "Active volume boundaries:" << std::endl;
      std::cout << "minx: " << fminx << "  maxx: " << fmaxx << std::endl;
      std::cout << "miny: " << fminy << "  maxy: " << fmaxy << std::endl;
      std::cout << "minz: " << fminz << "  maxz: " << fmaxz << std::endl;

      TVector3 Cathode_centre(geo->TPC(0, 0).GetCathodeCenter().X(),
                              (fminy + fmaxy) / 2, (fminz + fmaxz) / 2);
      std::cout << "Cathode_centre: " << Cathode_centre.X()
                << "  " << Cathode_centre.Y() << "  " << Cathode_centre.Z() << std::endl;

      // std::cout << "\nInitialize acos_arr with " << acos_bins+1
      //           << " hence with a resolution of " << 1./acos_bins << std::endl;
      // for(size_t i=0; i<=acos_bins; ++i){
      //   acos_arr[i] = std::acos(i/double(acos_bins));
      // }

      for(size_t i = 0; i != pvs->NOpChannels(); i++) {
        double OpDetCenter_i[3];
        std::vector<double> OpDetCenter_v;
        geo->OpDetGeoFromOpDet(i).GetCenter(OpDetCenter_i);
        OpDetCenter_v.assign(OpDetCenter_i, OpDetCenter_i + 3);
        fOpDetCenter.push_back(OpDetCenter_v);
        int type_i = -1;
        if(strcmp(geo->OpDetGeoFromOpDet(i).Shape()->IsA()->GetName(), "TGeoBBox") == 0) {
          type_i = 0;//Arapucas
          fOpDetLength.push_back(geo->OpDetGeoFromOpDet(i).Length());
          fOpDetHeight.push_back(geo->OpDetGeoFromOpDet(i).Height());
        }
        else {
          type_i = 1;//PMTs
          //    std::cout<<"Radio: "<<geo->OpDetGeoFromOpDet(i).RMax()<<std::endl;
          fOpDetLength.push_back(-1);
          fOpDetHeight.push_back(-1);
        }
        fOpDetType.push_back(type_i);
        // std::cout <<"OpChannel: "<<i<<"  Optical_Detector_Type: "<< type_i <<"  APERTURE_height: "
                  // <<geo->OpDetGeoFromOpDet(i).Height()<<"  APERTURE_width: "<<geo->OpDetGeoFromOpDet(i).Length()<< std::endl;
      }

      if(pvs->IncludePropTime()) {
        std::cout << "Using parameterisation of timings." << std::endl;
        //OLD VUV time parapetrization (to be removed soon)
        //pvs->SetDirectLightPropFunctions(functions_vuv, fd_break, fd_max, ftf1_sampling_factor);
        //pvs->SetReflectedCOLightPropFunctions(functions_vis, ft0_max, ft0_break_point);
        //New VUV time parapetrization
        pvs->LoadTimingsForVUVPar(fparameters, fstep_size, fmax_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance);

        // create vector of empty TF1s that will be replaces with the parameterisations that are generated as they are required
        // default TF1() constructor gives function with 0 dimensions, can then check numDim to qucikly see if a parameterisation has been generated
        const size_t num_params = (fmax_d - 25) / fstep_size;  // for d < 25cm, no parameterisaton, a delta function is used instead // TODO: unhardcode
        std::vector<TF1> VUV_timing_temp(num_params, TF1());
        VUV_timing = VUV_timing_temp;

        // initialise vectors to contain range parameterisations sampled to in each case
        // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling is regenerated, this is the slow part!
        std::vector<double> VUV_empty(num_params, 0);
        VUV_max = VUV_empty;
        VUV_min = VUV_empty;

        // VIS time parameterisation
        if (pvs->StoreReflected()) {
          // load parameters
          pvs->LoadTimingsForVISPar(fdistances_refl, fcut_off_pars, ftau_pars,
                                    fvis_vmean, fn_LAr_vis, fn_LAr_vuv);
        }
      }
      if(pvs->UseNhitsModel()) {
        std::cout << "Using semi-analytic model for number of hits:" << std::endl;
        fUseNhitsModel = true;
        // LAr absorption length in cm
        std::map<double, double> abs_length_spectrum = lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
        std::vector<double> x_v, y_v;
        for(auto elem : abs_length_spectrum) {
          x_v.push_back(elem.first);
          y_v.push_back(elem.second);
        }
        fL_abs_vuv =  interpolate(x_v, y_v, 9.7, false);// TODO: unhardcode //TODO: unsafe casting double to int

        // Load Gaisser-Hillas corrections for VUV semi-analytic hits
        std::cout << "Loading the GH corrections" << std::endl;
        pvs->LoadGHForVUVCorrection(fGHvuvpars, fborder_corr, fradius);
        fdelta_angulo = 10.; // angle bin size
        //Needed for Nhits-model border corrections (in cm)
        fYactive_corner = (fmaxy - fminy) / 2;
        fZactive_corner = (fmaxz - fminz) / 2;

        fYcathode = Cathode_centre.Y();
        fZcathode = Cathode_centre.Z();
        fReference_to_corner = std::sqrt(fYactive_corner*fYactive_corner +
                                         fZactive_corner*fZactive_corner);

        std::cout << "For border corrections: "
                  << fborder_corr[0] << "  " << fborder_corr[1] << std::endl;
        std::cout << "Photocathode-plane centre (z,y) = (" << fZcathode << ", " << fYcathode << ") "
                  << "and corner (z, y) = (" << fZactive_corner << ", " << fYactive_corner << ")" << std::endl;
        std::cout << "Reference_to_corner: " << fReference_to_corner << std::endl;

        if(pvs->StoreReflected()) {
          // Load corrections for VIS semi-anlytic hits
          std::cout << "Loading vis corrections" << std::endl;
          pvs->LoadParsForVISCorrection(fvispars, fradius);
          fStoreReflected = true;

          if (pvs->ApplyVISBorderCorrection()) {
            // load border corrections
            std::cout << "Loading vis border corrections" << std::endl;
            pvs->LoadParsForVISBorderCorrection(fvis_border_distances_x,
                                                fvis_border_distances_r, fvis_border_correction);
            fApplyVisBorderCorrection = true;
            fVisBorderCorrectionType = pvs->VISBorderCorrectionType();
          }
          else fApplyVisBorderCorrection = false;

          // cathode dimensions required for corrections
          fcathode_centre = geo->TPC(0, 0).GetCathodeCenter();
          fcathode_centre[1] = (fmaxy + fminy) / 2; fcathode_centre[2] = (fmaxz + fminz) / 2; // to get full cathode dimension rather than just single tpc
          fcathode_ydimension = fmaxy - fminy;
          fcathode_zdimension = fmaxz - fminz;
          // set cathode plane struct for solid angle function
          cathode_plane.h = fcathode_ydimension; cathode_plane.w = fcathode_zdimension;
          fplane_depth = std::abs(fcathode_centre[0]);
        }
        else fStoreReflected = false;
      }
      else fUseNhitsModel = false;
    }
    tpbemission = lar::providerFrom<detinfo::LArPropertiesService>()->TpbEm();
    const size_t nbins = tpbemission.size();
    double * parent = new double[nbins];
    size_t ii = 0;
    for(auto iter = tpbemission.begin(); iter != tpbemission.end(); ++iter) {
      parent[ii++] = (*iter).second;
    }
    rgen0 = new CLHEP::RandGeneral(parent, nbins);
    delete [] parent;
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
    // double stepEnergy          = larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit()/CLHEP::MeV;
    RecordPhotonsProduced(aStep, MeanNumberOfPhotons);//, stepEnergy);
    if (verboseLevel > 0) {
      G4cout << "\n Exiting from OpFastScintillation::DoIt -- NumberOfSecondaries = "
             << aParticleChange.GetNumberOfSecondaries() << G4endl;
    }
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }


  void OpFastScintillation::ProcessStep( const G4Step& step)
  {
    if(step.GetTotalEnergyDeposit() <= 0) return;

    OpDetPhotonTable::Instance()->AddEnergyDeposit
    (-1,
     -1,
     1.0,  //scintillation yield
     (double)(step.GetTotalEnergyDeposit() / CLHEP::MeV), //energy in MeV
     (float)(step.GetPreStepPoint()->GetPosition().x() / CLHEP::cm),
     (float)(step.GetPreStepPoint()->GetPosition().y() / CLHEP::cm),
     (float)(step.GetPreStepPoint()->GetPosition().z() / CLHEP::cm),
     (float)(step.GetPostStepPoint()->GetPosition().x() / CLHEP::cm),
     (float)(step.GetPostStepPoint()->GetPosition().y() / CLHEP::cm),
     (float)(step.GetPostStepPoint()->GetPosition().z() / CLHEP::cm),
     (double)(step.GetPreStepPoint()->GetGlobalTime()),
     (double)(step.GetPostStepPoint()->GetGlobalTime()),
     //step.GetTrack()->GetTrackID(),
     ParticleListAction::GetCurrentTrackID(),
     step.GetTrack()->GetParticleDefinition()->GetPDGEncoding(),
     step.GetPreStepPoint()->GetPhysicalVolume()->GetName()
    );
  }


  bool OpFastScintillation::RecordPhotonsProduced(const G4Step& aStep, double MeanNumberOfPhotons)//, double stepEnergy)
  {
    // make sure that whatever happens afterwards, the energy deposition is stored
    art::ServiceHandle<sim::LArG4Parameters const> lgp;
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
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    size_t const NOpChannels = pvs->NOpChannels();

    G4int nscnt = 1;
    if (Fast_Intensity && Slow_Intensity) nscnt = 2;

    double Num = 0;
    double YieldRatio = 0;

    if (scintillationByParticleType) {
      // The scintillation response is a function of the energy
      // deposited by particle types.

      // Get the definition of the current particle
      G4ParticleDefinition *pDef = aParticle->GetDefinition();

      // Obtain the G4MaterialPropertyVectory containing the
      // scintillation light yield as a function of the deposited
      // energy for the current particle type

      // Protons
      if(pDef == G4Proton::ProtonDefinition()) {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("PROTONYIELDRATIO");
      }
      // Muons
      else if(pDef == G4MuonPlus::MuonPlusDefinition() ||
              pDef == G4MuonMinus::MuonMinusDefinition()) {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("MUONYIELDRATIO");
      }
      // Pions
      else if(pDef == G4PionPlus::PionPlusDefinition() ||
              pDef == G4PionMinus::PionMinusDefinition()) {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("PIONYIELDRATIO");
      }
      // Kaons
      else if(pDef == G4KaonPlus::KaonPlusDefinition() ||
              pDef == G4KaonMinus::KaonMinusDefinition()) {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("KAONYIELDRATIO");
      }
      // Alphas
      else if(pDef == G4Alpha::AlphaDefinition()) {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("ALPHAYIELDRATIO");
      }
      // Electrons (must also account for shell-binding energy
      // attributed to gamma from standard PhotoElectricEffect)
      else if(pDef == G4Electron::ElectronDefinition() ||
              pDef == G4Gamma::GammaDefinition()) {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("ELECTRONYIELDRATIO");
      }
      // Default for particles not enumerated/listed above
      else {
        YieldRatio = aMaterialPropertiesTable->
                     GetConstProperty("ELECTRONYIELDRATIO");
      }
      // If the user has not specified yields for (p,d,t,a,carbon)
      // then these unspecified particles will default to the
      // electron's scintillation yield
      if(YieldRatio == 0) {
        YieldRatio = aMaterialPropertiesTable->
          GetConstProperty("ELECTRONYIELDRATIO");
      }
    }

    double const xyz[3] = { x0[0]/CLHEP::cm, x0[1]/CLHEP::cm, x0[2]/CLHEP::cm };
    auto const& Visibilities = pvs->GetAllVisibilities(xyz);

    phot::MappedCounts_t ReflVisibilities;

    // Store timing information in the object for use in propagationTime method
    if(pvs->StoreReflected()) {
      ReflVisibilities = pvs->GetAllVisibilities(xyz, true);
      if(pvs->StoreReflT0())
        ReflT0s = pvs->GetReflT0s(xyz);
    }
    if(pvs->IncludeParPropTime()) {
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
          if(Fast_Intensity) {
            ScintillationTime   = aMaterialPropertiesTable->
                                  GetConstProperty("FASTTIMECONSTANT");
            if (fFiniteRiseTime) {
              ScintillationRiseTime = aMaterialPropertiesTable->
                                      GetConstProperty("FASTSCINTILLATIONRISETIME");
            }
            ScintillationIntegral =
              (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
          }
          if(Slow_Intensity) {
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
          if(YieldRatio == 0)
            YieldRatio = aMaterialPropertiesTable->
                         GetConstProperty("YIELDRATIO");

          if(ExcitationRatio == 1.0) {
            Num = std::min(YieldRatio, 1.0) * MeanNumberOfPhotons;
          }
          else{
            Num = std::min(ExcitationRatio, 1.0) * MeanNumberOfPhotons;
          }
          ScintillationTime = aMaterialPropertiesTable->
                              GetConstProperty("FASTTIMECONSTANT");
          if(fFiniteRiseTime) {
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
        if(fFiniteRiseTime) {
          ScintillationRiseTime = aMaterialPropertiesTable->
                                  GetConstProperty("SLOWSCINTILLATIONRISETIME");
        }
        ScintillationIntegral =
          (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
      }

      if(!ScintillationIntegral) continue;
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

      if(!Visibilities && !pvs->UseNhitsModel()) {
      }
      else {
        std::map<size_t, int> DetectedNum;
        std::map<size_t, int> ReflDetectedNum;

        for(size_t OpDet = 0; OpDet != NOpChannels; ++OpDet) {
          G4int DetThisPMT = 0.;
          if(Visibilities && !pvs->UseNhitsModel()) {
            DetThisPMT = G4int(G4Poisson(Visibilities[OpDet] * Num));
          }
          else {
            TVector3 ScintPoint( xyz[0], xyz[1], xyz[2] );
            TVector3 OpDetPoint(fOpDetCenter.at(OpDet)[0],
                                fOpDetCenter.at(OpDet)[1],
                                fOpDetCenter.at(OpDet)[2]);
            fydimension = fOpDetHeight.at(OpDet);
            fzdimension = fOpDetLength.at(OpDet);
            // set detector struct for solid angle function
            detPoint.h = fydimension; detPoint.w = fzdimension;
            DetThisPMT = VUVHits(Num, ScintPoint, OpDetPoint, fOpDetType.at(OpDet));
          }

          if(DetThisPMT > 0) {
            DetectedNum[OpDet] = DetThisPMT;
            //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
            //   //   it->second<<" " << Num << " " << DetThisPMT;
            //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
          }
          if(pvs->StoreReflected()) {
            G4int ReflDetThisPMT = 0;
            if (!pvs->UseNhitsModel()) {
              ReflDetThisPMT = G4int(G4Poisson(ReflVisibilities[OpDet] * Num));
            }
            else {
              TVector3 ScintPoint( xyz[0], xyz[1], xyz[2] );
              TVector3 OpDetPoint(fOpDetCenter.at(OpDet)[0],
                                  fOpDetCenter.at(OpDet)[1],
                                  fOpDetCenter.at(OpDet)[2]);
              ReflDetThisPMT = VISHits(Num, ScintPoint, OpDetPoint, fOpDetType.at(OpDet));
            }
            if(ReflDetThisPMT > 0) {
              ReflDetectedNum[OpDet] = ReflDetThisPMT;
            }
          }
        }

        std::vector<double> arrival_time_dist;
        // Now we run through each PMT figuring out num of detected photons
        for (size_t Reflected = 0; Reflected <= 1; ++Reflected) {
          // Only do the reflected loop if we have reflected visibilities
          if (Reflected && !pvs->StoreReflected()) continue;

          std::map<size_t, int>::const_iterator itstart;
          std::map<size_t, int>::const_iterator itend;
          if (Reflected) {
            itstart = ReflDetectedNum.begin();
            itend   = ReflDetectedNum.end();
          }
          else {
            itstart = DetectedNum.begin();
            itend   = DetectedNum.end();
          }
          for(auto itdetphot = itstart; itdetphot != itend; ++itdetphot) {
            const size_t OpChannel = itdetphot->first;
            const int NPhotons  = itdetphot->second;

            // Set up the OpDetBTR information
            sim::OpDetBacktrackerRecord tmpOpDetBTRecord(OpChannel);
            int thisG4TrackID = ParticleListAction::GetCurrentTrackID();
            double xyzPos[3];
            average_position(aStep, xyzPos);
            double Edeposited  = 0;
            if(scintillationByParticleType) {
              //We use this when it is the only sensical information. It may be of limited use to end users.
              Edeposited = aStep.GetTotalEnergyDeposit();
            }
            else if(emSaturation) {
              //If Birk Coefficient used, log VisibleEnergies.
              Edeposited = larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit() / CLHEP::MeV;
            }
            else {
              //We use this when it is the only sensical information. It may be of limited use to end users.
              Edeposited = aStep.GetTotalEnergyDeposit();
            }

            // Get the transport time distribution
            arrival_time_dist.resize(NPhotons);
            propagationTime(arrival_time_dist, x0, OpChannel, Reflected);

            //We need to split the energy up by the number of photons so that we never try to write a 0 energy.
            Edeposited = Edeposited / double(NPhotons);

            // Loop through the photons
            for (G4int i = 0; i < NPhotons; ++i) {
              //std::cout<<"VUV time correction: "<<arrival_time_dist[i]<<std::endl;
              G4double Time = t0
                              + scint_time(aStep, ScintillationTime, ScintillationRiseTime)
                              + arrival_time_dist[i] * CLHEP::ns;

              // Always store the BTR
              tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, Time, 1, xyzPos, Edeposited);

              // Store as lite photon or as OnePhoton
              if(lgp->UseLitePhotons()) {
                fst->AddLitePhoton(OpChannel, static_cast<int>(Time), 1, Reflected);
              }
              else {
                // The sim photon in this case stores its production point and time
                TVector3 PhotonPosition( x0[0], x0[1], x0[2] );

                float PhotonEnergy = 0;
                if (Reflected)  PhotonEnergy = reemission_energy() * CLHEP::eV;
                else            PhotonEnergy = 9.7 * CLHEP::eV;// TODO: unhardcode

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
    for (G4int i = 0 ; i < numOfMaterials; i++) {
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
            InsertValues(currentPM, currentCII);

            // Set previous values to current ones prior to loop
            G4double prevPM  = currentPM;
            G4double prevCII = currentCII;
            G4double prevIN  = currentIN;

            // loop over all (photon energy, intensity)
            // pairs stored for this material
            for (size_t i = 1;
                 i < theFastLightVector->GetVectorLength();
                 i++) {
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
            InsertValues(currentPM, currentCII);

            // Set previous values to current ones prior to loop
            G4double prevPM  = currentPM;
            G4double prevCII = currentCII;
            G4double prevIN  = currentIN;

            // loop over all (photon energy, intensity)
            // pairs stored for this material
            for (size_t i = 1;
                 i < theSlowLightVector->GetVectorLength();
                 i++) {
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
      theFastIntegralTable->insertAt(i, aPhysicsOrderedFreeVector);
      theSlowIntegralTable->insertAt(i, bPhysicsOrderedFreeVector);
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


  G4double OpFastScintillation::GetMeanFreePath(const G4Track&,
      G4double,
      G4ForceCondition* condition)
  {
    *condition = StronglyForced;
    return DBL_MAX;
  }


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
    G4double avgVelocity = (pPreStepPoint->GetVelocity() +
                            pPostStepPoint->GetVelocity()) / 2.;
    G4double deltaTime = aStep.GetStepLength() / avgVelocity;
    if (ScintillationRiseTime == 0.0) {
      deltaTime = deltaTime -
                  ScintillationTime * std::log( G4UniformRand() );
    }
    else {
      deltaTime = deltaTime +
                  sample_time(ScintillationRiseTime, ScintillationTime);
    }
    return deltaTime;
  }


  void OpFastScintillation::propagationTime(std::vector<double>& arrival_time_dist, G4ThreeVector x0,
                                             const size_t OpChannel, bool Reflected) //const
  {
    static art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    if (pvs->IncludeParPropTime() && pvs->IncludePropTime()) {
      throw cet::exception("OpFastScintillation") << "Cannot have both propagation time models simultaneously.";
    }
    else if (pvs->IncludeParPropTime() && !(ParPropTimeTF1  && (ParPropTimeTF1[OpChannel].GetNdim() == 1)) ) {
      //Warning: TF1::GetNdim()==1 will tell us if the TF1 is really defined or it is the default one.
      //This will fix a segfault when using timing and interpolation.
      G4cout << "WARNING: Requested parameterized timing, but no function found. Not applying propagation time." << G4endl;
    }
    else if (pvs->IncludeParPropTime()) {
      if (Reflected)
        throw cet::exception("OpFastScintillation") << "No parameterized propagation time for reflected light";
      for(size_t i = 0; i < arrival_time_dist.size(); ++i) {
        arrival_time_dist[i] = ParPropTimeTF1[OpChannel].GetRandom();
      }
    }
    else if (pvs->IncludePropTime()) {
      // Get VUV photons arrival time distribution from the parametrization
      const G4ThreeVector OpDetPoint(fOpDetCenter.at(OpChannel)[0]*CLHEP::cm,
                                     fOpDetCenter.at(OpChannel)[1]*CLHEP::cm,
                                     fOpDetCenter.at(OpChannel)[2]*CLHEP::cm);
      if (!Reflected) {
        double distance_in_cm = (x0 - OpDetPoint).mag() / CLHEP::cm; // this must be in CENTIMETERS!
        getVUVTimes(arrival_time_dist, distance_in_cm); // in ns
      }
      else {
        TVector3 ScintPoint( x0[0]/CLHEP::cm, x0[1]/CLHEP::cm, x0[2]/CLHEP::cm ); // in cm
        TVector3 OpDetPoint_tv3(fOpDetCenter.at(OpChannel)[0],
                                fOpDetCenter.at(OpChannel)[1],
                                fOpDetCenter.at(OpChannel)[2]); // in cm
        getVISTimes(arrival_time_dist, ScintPoint, OpDetPoint_tv3); // in ns
      }
    }
  }


  G4double OpFastScintillation::sample_time(const G4double tau1, const G4double tau2) const
  {
    // tau1: rise time and tau2: decay time
    while(1) {
      // two random numbers
      G4double ran1 = G4UniformRand();
      G4double ran2 = G4UniformRand();
      //
      // exponential distribution as envelope function: very efficient
      //
      G4double d = (tau1 + tau2) / tau2;
      // make sure the envelope function is
      // always larger than the bi-exponential
      G4double t = -1.0 * tau2 * std::log(1 - ran1);
      G4double g = d * single_exp(t, tau2);
      if (ran2 <= bi_exp(t, tau1, tau2) / g) return t;
    }
    return -1.0;
  }


  double OpFastScintillation::reemission_energy() const
  {
    return rgen0->fire() * ((*(--tpbemission.end())).first - (*tpbemission.begin()).first) +
      (*tpbemission.begin()).first;
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
  /*
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
                               std::pow(10.,functions_vuv[0]->Eval(distance))};
      if(distance > fd_break) {
        pars_landau[0]=functions_vuv[6]->Eval(distance);
        pars_landau[1]=functions_vuv[2]->Eval(fd_break);
        pars_landau[2]=std::pow(10.,functions_vuv[5]->Eval(distance));
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
      TF1 *fint = new TF1("fint","finter_d",flandau->GetMaximumX(),3*t_direct,5);
      double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
      fint->SetParameters(parsInt);
      double t_int = fint->GetMinimumX();
      double minVal = fint->Eval(t_int);
      //the functions must intersect!!!
      if(minVal>0.015)
        G4cout<<"WARNING: Parametrization of Direct Light discontinuous (landau + expo)!!!!!!"<<G4endl;

      TF1 *fVUVTiming =  new TF1("fTiming","LandauPlusExpoFinal",0,signal_t_range,6);
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
                               std::pow(10.,functions_vis[0]->Eval(t0))};
      double pars_expo[2] = {functions_vis[3]->Eval(t0), functions_vis[4]->Eval(t0)};
      if(t0 > ft0_break_point) {
        pars_landau[0] = -0.798934 + 1.06216*t0;
        pars_landau[1] = functions_vis[2]->Eval(ft0_break_point);
        pars_landau[2] = std::pow(10.,functions_vis[0]->Eval(ft0_break_point));
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
    }*/


  // New Parametrization code
  // parameterisation generation function
  void OpFastScintillation::generateParam(const size_t index)
  {
    // get distance
    double distance_in_cm = (index * fstep_size) + 25;// TODO: unhardcode

    // time range
    const double signal_t_range = 5000.;// TODO: unhardcode

    // parameterisation TF1
    TF1 fVUVTiming;

    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
    double t_direct_min = distance_in_cm / fvuv_vgroup_max;

    // Defining the model function(s) describing the photon transportation timing vs distance
    // Getting the landau parameters from the time parametrization
    double pars_landau[3];
    interpolate(pars_landau, fparameters[0], fparameters[2], fparameters[3],
                fparameters[1], distance_in_cm, true);
    // Deciding which time model to use (depends on the distance)
    // defining useful times for the VUV arrival time shapes
    if(distance_in_cm >= finflexion_point_distance) {
      double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1],
                            pars_landau[2]};
      // Set model: Landau
      fVUVTiming = TF1("fVUVTiming", model_far, 0, signal_t_range, 4);
      fVUVTiming.SetParameters(pars_far);
    }
    else {
      // Set model: Landau + Exponential
      fVUVTiming = TF1("fVUVTiming", model_close, 0, signal_t_range, 7);
      // Exponential parameters
      double pars_expo[2];
      // Getting the exponential parameters from the time parametrization
      pars_expo[1] = interpolate(fparameters[4], fparameters[5],
                                 distance_in_cm, true);
      //For simplicity, not considering the small dependency with the offset angle in pars_expo[0]
      //Using the value for the [30,60deg] range. fparameters[6] and fparameters[8] are the values
      //for [0,30deg] range and [60,90deg] range respectively
      pars_expo[0] = fparameters[7].at(0) + fparameters[7].at(1) * distance_in_cm;
      pars_expo[0] *= pars_landau[2];
      pars_expo[0] = std::log(pars_expo[0]);
      // this is to find the intersection point between the two functions:
      TF1 fint = TF1("fint", finter_d, pars_landau[0], 4 * t_direct_mean, 5);
      double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2],
                           pars_expo[0], pars_expo[1]};
      fint.SetParameters(parsInt);
      double t_int = fint.GetMinimumX();
      double minVal = fint.Eval(t_int);
      // the functions must intersect - output warning if they don't
      if(minVal > 0.015) {
        std::cout << "WARNING: Parametrization of VUV light discontinuous for distance = "
                  << distance_in_cm << std::endl;
        std::cout << "WARNING: This shouldn't be happening " << std::endl;
      }
      double parsfinal[7] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2],
                             pars_expo[0], pars_expo[1], t_direct_min};
      fVUVTiming.SetParameters(parsfinal);
    }

    // set the number of points used to sample parameterisation
    // for shorter distances, peak is sharper so more sensitive sampling required
    int fsampling;// TODO: unhardcode
    if (distance_in_cm < 50) {fsampling = 10000;}
    else if (distance_in_cm < 100) {fsampling = 5000;}
    else {fsampling = 1000;}
    fVUVTiming.SetNpx(fsampling);

    // calculate max and min distance relevant to sample parameterisation
    // max
    // TODO: array instead of pointer? why this?
    const size_t nq_max = 1;
    double xq_max[nq_max];
    double yq_max[nq_max];
    xq_max[0] = 0.99;   // include 99%
    fVUVTiming.GetQuantiles(nq_max, yq_max, xq_max);
    double max = yq_max[0];
    // min
    double min = t_direct_min;

    // generate the sampling
    // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
    // all subsequent calls check if it has been generated previously and are ~100+ times quicker
    // add timing to the vector of timings and range to vectors of ranges
    VUV_timing[index] = fVUVTiming;
    VUV_max[index] = max;
    VUV_min[index] = min;
  }


  // VUV arrival times calculation function
  void OpFastScintillation::getVUVTimes(std::vector<double>& arrivalTimes, double distance)
  {
    if (distance < 25) {// TODO: unhardcode
      // times are fixed shift i.e. direct path only
      double t_prop_correction = distance / fvuv_vgroup_mean;
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = t_prop_correction;
      }
    }
    else { // distance >= 25cm
      // determine nearest parameterisation in discretisation
      int index = std::round((distance - 25) / fstep_size);// TODO: unhardcode
      // check whether required parameterisation has been generated, generating if not
      if (VUV_timing[index].GetNdim() == 0) {
        generateParam(index);
      }
      // randomly sample parameterisation for each photon
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = VUV_timing[index].GetRandom(VUV_min[index], VUV_max[index]);
      }
    }
  }


  // VIS arrival times calculation functions
  void OpFastScintillation::getVISTimes(std::vector<double>& arrivalTimes, TVector3 ScintPoint,
                                        TVector3 OpDetPoint)
  {
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

    // calculate point of reflection for shortest path accounting for difference
    // in refractive indices vectors for storing results
    TVector3 image(0, 0, 0);
    TVector3 bounce_point(0, 0, 0);

    // distance to wall
    TVector3 v_to_wall(plane_depth - ScintPoint[0], 0, 0);

    // hotspot is point on wall where TPB is activated most intensely by the scintillation
    TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);

    // define "image" by reflecting over plane
    image = hotspot + v_to_wall * (fn_LAr_vis / fn_LAr_vuv);

    // find point of intersection with plane j of ray from the PMT to the image
    TVector3 tempvec = (OpDetPoint - image).Unit();
    double tempnorm = ((image - hotspot).Mag()) / std::abs(tempvec[0]);
    bounce_point = image + tempvec * tempnorm;

    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point - ScintPoint).Mag();
    double Visdist = (OpDetPoint - bounce_point).Mag();

    // calculate times taken by each part
    getVUVTimes(arrivalTimes, VUVdist);

    // sum parts to get total transport times times
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] += Visdist/fvis_vmean;
    }

    // *************************************************************************************************
    //      Smearing of arrival time distribution
    // *************************************************************************************************
    // calculate fastest time possible
    // vis part
    double vis_time = Visdist / fvis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < 25) {// TODO: unhardcode
      vuv_time = VUVdist / fvuv_vgroup_mean;
    }
    else {
      // find index of required parameterisation
      const size_t index = std::round((VUVdist - 25) / fstep_size);// TODO: unhardcode
      // find shortest time
      vuv_time = VUV_min[index];
    }
    // sum
    double fastest_time = vis_time + vuv_time;

    // calculate angle alpha between scintillation point and reflection point
    double cosine_alpha = std::abs(ScintPoint[0] - bounce_point[0]) / VUVdist;
    // double alpha = std::acos(cosine_alpha) * 180. / CLHEP::pi;
    double alpha = fast_acos(cosine_alpha) * 180. / CLHEP::pi;

    // determine smearing parameters using interpolation of generated points:
    // 1). tau = exponential smearing factor, varies with distance and angle
    // 2). cutoff = largest smeared time allowed, preventing excessively large
    //     times caused by exponential distance to cathode
    double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
    // angular bin
    size_t alpha_bin = alpha / 10;// TODO: unhardcode
    if (alpha_bin >= ftau_pars.size()) {
      alpha_bin = ftau_pars.size() - 1;      // default to the largest available bin if alpha larger than parameterised region; i.e. last bin effectively [last bin start value, 90] deg bin
    }
    // cut-off and tau
    double cutoff = interpolate(fdistances_refl, fcut_off_pars[alpha_bin],
                                 distance_cathode_plane, true);
    double tau = interpolate(fdistances_refl, ftau_pars[alpha_bin],
                             distance_cathode_plane, true);
    if (tau < 0) {tau = 0;} // if tau extrapolate goes wrong

    // apply smearing:
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      double arrival_time_smeared;
      // if time is already greater than cutoff, do not apply smearing
      if (arrivalTimes[i]  >= cutoff) {
        continue;
      }
      // otherwise smear
      else {
        unsigned int counter = 0;
        // loop until time generated is within cutoff limit
        // most are within single attempt, very few take more than two
        do {
          // don't attempt smearings too many times
          if (counter >= 10) {// TODO: unhardcode
            arrival_time_smeared = arrivalTimes[i]; // don't smear
            break;
          }
          else {
            // generate random number in appropriate range
            double x = gRandom->Uniform(0.5, 1.0);// TODO: unhardcode
            // apply the exponential smearing
            arrival_time_smeared = arrivalTimes[i] +
              (arrivalTimes[i] - fastest_time)*(std::pow(x,-tau) - 1);
          }
          counter++;
        }
        while (arrival_time_smeared > cutoff);
      }
      arrivalTimes[i] = arrival_time_smeared;
    }
  }


  // VUV semi-analytic hits calculation
  int OpFastScintillation::VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type)
  {
    // check optical channel is in same TPC as scintillation light, if not return 0 hits
    // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero, otherwise return 0 hits
    if (((ScintPoint[0] < 0) != (OpDetPoint[0] < 0)) && std::abs(OpDetPoint[0]) > 10) {// TODO: unhardcode
      return 0;
    }
    //semi-analytic approach only works in the active volume
    if((ScintPoint[0] < fminx) || (ScintPoint[0] > fmaxx) ||
       (ScintPoint[1] < fminy) || (ScintPoint[1] > fmaxy) ||
       (ScintPoint[2] < fminz) || (ScintPoint[2] > fmaxz) ||
       (std::abs(ScintPoint[0]) <= fplane_depth)) {
      return 0;
    }

    // distance and angle between ScintPoint and OpDetPoint
    double distance = dist(&ScintPoint[0], &OpDetPoint[0], 3);
    double cosine = std::abs(ScintPoint[0] - OpDetPoint[0]) / distance;
    // double theta = std::acos(cosine) * 180. / CLHEP::pi;
    double theta = fast_acos(cosine) * 180. / CLHEP::pi;

    // calculate solid angle:
    double solid_angle = 0;
    // Arapucas
    if (optical_detector_type == 0) {
      // get scintillation point coordinates relative to arapuca window centre
      std::array<double, 3> ScintPoint_rel = {std::abs(ScintPoint.X() - OpDetPoint.X()),
                                              std::abs(ScintPoint.Y() - OpDetPoint.Y()),
                                              std::abs(ScintPoint.Z() - OpDetPoint.Z())};
      // calculate solid angle
      solid_angle = Rectangle_SolidAngle(detPoint, ScintPoint_rel);
    }
    // PMTs
    else if (optical_detector_type == 1) {
      // offset in z-y plane
      double d = dist(&ScintPoint[1], &OpDetPoint[1], 2);
      // drift distance (in x)
      double h =  std::abs(ScintPoint[0] - OpDetPoint[0]);
      // Solid angle of a disk
      solid_angle = Disk_SolidAngle(d, h, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" << std:: endl;
    }

    // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
    double hits_geo = std::exp(-1.*distance / fL_abs_vuv) * (solid_angle / (4 * CLHEP::pi)) * Nphotons_created;

    // apply Gaisser-Hillas correction for Rayleigh scattering distance
    // and angular dependence offset angle bin
    const size_t j = (theta / fdelta_angulo);// TODO:: std::round?

    //Accounting for border effects
    double z_to_corner = std::abs(ScintPoint[2] - fZactive_corner) - fZactive_corner;
    double y_to_corner = std::abs(ScintPoint[1]) - fYactive_corner;
    double distance_to_corner = std::sqrt(y_to_corner * y_to_corner +
                                          z_to_corner * z_to_corner); // in the ph-cathode plane
    double pars_ini_[4] = {fGHvuvpars[0][j] + fborder_corr[0] * (distance_to_corner - fReference_to_corner),
                           fGHvuvpars[1][j] + fborder_corr[1] * (distance_to_corner - fReference_to_corner),
                           fGHvuvpars[2][j],
                           fGHvuvpars[3][j]
                          };
    double GH_correction = Gaisser_Hillas(distance, pars_ini_);
    double hits_rec = gRandom->Poisson( GH_correction * hits_geo / cosine );
    // round to integer value, cannot have non-integer number of hits
    int hits_vuv = std::round(hits_rec);

    return hits_vuv;
  }


  // VIS hits semi-analytic model calculation
  int OpFastScintillation::VISHits(int Nphotons_created, TVector3 ScintPoint,
                                   TVector3 OpDetPoint, int optical_detector_type)
  {
    // check optical channel is in same TPC as scintillation light, if not return 0 hits
    // temporary method working for SBND, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero, otherwise return 0 hits
    if (((ScintPoint[0] < 0) != (OpDetPoint[0] < 0)) && std::abs(OpDetPoint[0]) > 10) {// TODO: unhardcode
      return 0;
    }

    //semi-analytic approach only works in the active volume
    if((ScintPoint[0] < fminx) || (ScintPoint[0] > fmaxx) ||
       (ScintPoint[1] < fminy) || (ScintPoint[1] > fmaxy) ||
       (ScintPoint[2] < fminz) || (ScintPoint[2] > fmaxz) ||
       (std::abs(ScintPoint[0]) <= fplane_depth)) {
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
    // get scintpoint coords relative to centre of cathode plane
    std::array<double, 3> ScintPoint_relative = {std::abs(ScintPoint.X() - plane_depth),
                                                 std::abs(ScintPoint.Y() - fcathode_centre[1]),
                                                 std::abs(ScintPoint.Z() - fcathode_centre[2])};
    // calculate solid angle of cathode from the scintillation point
    double solid_angle_cathode = Rectangle_SolidAngle(cathode_plane, ScintPoint_relative);
    // calculate distance and angle between ScintPoint and hotspot
    // vast majority of hits in hotspot region directly infront of scintpoint,
    // therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
    double distance_cathode = std::abs(plane_depth - ScintPoint[0]);
    double cosine_cathode = 1;
    double theta_cathode = 0;
    // calculate hits on cathode plane via geometric acceptance
    double cathode_hits_geo = std::exp(-1.*distance_cathode / fL_abs_vuv) *
      (solid_angle_cathode / (4.*CLHEP::pi)) * Nphotons_created;
    // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
    // offset angle bin
    const size_t j = (theta_cathode / fdelta_angulo);// TODO:: std::round?
    double  pars_ini_[4] = {fGHvuvpars[0][j],
                            fGHvuvpars[1][j],
                            fGHvuvpars[2][j],
                            fGHvuvpars[3][j]};
    double GH_correction = Gaisser_Hillas(distance_cathode, pars_ini_);
    double cathode_hits_rec = GH_correction * cathode_hits_geo / cosine_cathode;

    // 2). calculate number of these hits which reach the optical detector from the hotspot via solid angle
    // hotspot coordinates
    TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);

    // calculate solid angle of optical channel
    double solid_angle_detector = 0;
    // rectangular aperture
    if (optical_detector_type == 0) {
      // get hotspot coordinates relative to detpoint
      std::array<double, 3> emission_relative = {std::abs(hotspot.X() - OpDetPoint.X()),
                                                 std::abs(hotspot.Y() - OpDetPoint.Y()),
                                                 std::abs(hotspot.Z() - OpDetPoint.Z())};
      // calculate solid angle
      solid_angle_detector = Rectangle_SolidAngle(detPoint, emission_relative);
    }
    // disk aperture
    else if (optical_detector_type == 1) {
      // offset in z-y plane
      double d = dist(&hotspot[1], &OpDetPoint[1], 2);
      // drift distance (in x)
      double h =  std::abs(hotspot[0] - OpDetPoint[0]);
      // calculate solid angle
      solid_angle_detector = Disk_SolidAngle(d, h, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" << std::endl;
    }

    // calculate number of hits via geometeric acceptance
    double hits_geo = (solid_angle_detector / (2.*CLHEP::pi)) * cathode_hits_rec; // 2*pi due to presence of reflective foils

    // calculate distances and angles for application of corrections
    // distance to hotspot
    double distance_vuv = dist(&ScintPoint[0], &hotspot[0], 3);
    // distance from hotspot to optical detector
    double distance_vis = dist(&hotspot[0], &OpDetPoint[0], 3);
     //  angle between hotspot and optical detector
    double cosine_vis = std::abs(hotspot[0] - OpDetPoint[0]) / distance_vis;
    // double theta_vis = std::acos(cosine_vis) * 180. / CLHEP::pi;
    double theta_vis = fast_acos(cosine_vis) * 180. / CLHEP::pi;
    const size_t k = (theta_vis / fdelta_angulo);// TODO:: std::round?

    // apply geometric correction
    double pars_ini_vis[6] = {fvispars[0][k], fvispars[1][k], fvispars[2][k],
                              fvispars[3][k], fvispars[4][k], fvispars[5][k]};
    double geo_correction = Pol_5(distance_vuv, pars_ini_vis);
    double hits_rec = gRandom->Poisson(geo_correction * hits_geo / cosine_vis);

    // apply border correction
    int hits_vis = 0;
    if (fApplyVisBorderCorrection) {
      // calculate distance for interpolation depending on model
      double r = 0;
      if (fVisBorderCorrectionType == "Radial") {
        r = dist(&ScintPoint[1], &fcathode_centre[1], 2);
      }
      else if (fVisBorderCorrectionType == "Vertical") {
        r = std::abs(ScintPoint[1]);
      }
      else {
        std::cout << "Invalid border correction type - defaulting to using central value" << std::endl;
      }
      // interpolate in x for each r bin
      const size_t nbins_r = fvis_border_correction[k].size();
      std::vector<double> interp_vals(nbins_r, 0.0);
      for (size_t i = 0; i < nbins_r; ++i) {
        interp_vals[i] = interpolate(fvis_border_distances_x, fvis_border_correction[k][i],
                                     std::abs(ScintPoint[0]), false);
      }
      // interpolate in r
      double border_correction = interpolate(fvis_border_distances_r, interp_vals, r, false);
      // apply border correction
      double hits_rec_borders = border_correction * hits_rec / cosine_vis;

      // round final result
      hits_vis = std::round(hits_rec_borders);
    }
    else {
      // round final result
      hits_vis = std::round(hits_rec);
    }

    return hits_vis;
  }


  constexpr G4double OpFastScintillation::single_exp(const G4double t, const G4double tau2) const
  {
    return std::exp(-1.0 * t / tau2) / tau2;
  }


  constexpr G4double OpFastScintillation::bi_exp(const G4double t, const G4double tau1,
                                                 const G4double tau2) const
  {// TODO: what's up with this? ... / tau2 / tau2 ...
    return std::exp(-1.0 * t / tau2) *
      (1 - std::exp(-1.0 * t / tau1)) / tau2 / tau2 * (tau1 + tau2);
  }


  constexpr G4double OpFastScintillation::Gaisser_Hillas(const double x,
                                                         const double *par)
  {
    double X_mu_0 = par[3];
    double Normalization = par[0];
    double Diff = par[1] - X_mu_0;
    double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
    double Exponential = std::exp((par[1] - x) / par[2]);
    return (Normalization * Term * Exponential);
  }


  constexpr double OpFastScintillation::Pol_5(const double x, double *par)
  {
    // 5th order polynomial function
    double xpow = 1.;
    for(unsigned i=1; i<=5; ++i){
      xpow *= x;
      par[0] += par[i] * xpow ;
    }
    return par[0];
  }


  double finter_d(double *x, double *par)
  {
    double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
    double y2 = TMath::Exp(par[3] + x[0] * par[4]);
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
    double y1 = par[3] * TMath::Landau(x[0], par[1], par[2]);
    double y2 = TMath::Exp(par[4] + x[0] * par[5]);
    if(x[0] > par[0]) y1 = 0.;
    if(x[0] < par[0]) y2 = 0.;
    return (y1 + y2);
  }


  double finter_r(double *x, double *par)
  {
    double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
    double y2 = par[5] * TMath::Landau(x[0], par[3], par[4]);
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
    double y1 = par[3] * TMath::Landau(x[0], par[1], par[2]);
    double y2 = TMath::Exp(par[4] + x[0] * par[5]);
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
    double y = par[3] * TMath::Landau(x[0], par[1], par[2]);
    if(x[0] <= par[0]) y = 0.;
    return y;
  }


  //======================================================================
  //   Returns interpolated value at x from parallel arrays ( xData, yData )
  //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
  //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
  double interpolate(std::vector<double> &xData, std::vector<double> &yData,
                     double x, bool extrapolate)
  {
    size_t size = xData.size();
    size_t i = 0; // find left end of interval for interpolation
    if (x >= xData[size-2]) { // special case: beyond right end
      i = size - 2;
    }
    else {
      while (x > xData[i+1]) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
    if (!extrapolate) { // if beyond ends of array and not extrapolating
      if (x < xL) yR = yL;
      if (x > xR) yL = yR;
    }
    double dydx = (yR - yL) / (xR - xL); // gradient
    return yL + dydx * ( x - xL ); // linear interpolation
  }


  void interpolate(double inter[], std::vector<double> &xData,
                   std::vector<double> &yData1, std::vector<double> &yData2,
                   std::vector<double> &yData3, double x, bool extrapolate)
  {
    size_t size = xData.size();
    size_t i = 0; // find left end of interval for interpolation
    if (x >= xData[size-2]) { // special case: beyond right end
      i = size - 2;
    }
    else {
      while (x > xData[i+1]) i++;
    }
    double xL = xData[i], xR = xData[i+1]; // points on either side (unless beyond ends)
    double yL1 = yData1[i], yR1 = yData1[i+1];
    double yL2 = yData2[i], yR2 = yData2[i+1];
    double yL3 = yData3[i], yR3 = yData3[i+1];

    if (!extrapolate) { // if beyond ends of array and not extrapolating
      if (x < xL) {
        yR1 = yL1;
        yR2 = yL2;
        yR3 = yL3;
      }
      if (x > xR) {
        yL1 = yR1;
        yL2 = yR2;
        yL3 = yR3;
      }
    }
    inter[0] = (x - xL) * (yR1 - yL1) / (xR - xL);
    inter[1] = (x - xL) * ( yR2 - yL2 ) / ( xR - xL );
    inter[2] = (x - xL) * ( yR3 - yL3 ) / ( xR - xL );

    inter[0] += yL1;
    inter[1] += yL2;
    inter[2] += yL3;
  }


  // solid angle of circular aperture
  // TODO: allow greater tolerance in comparisons, by default its using:
  // std::numeric_limits<double>::epsilon(): 2.22045e-16
  // that's an unrealistic small number, better setting
  // constexpr double tolerance = 0.0000001; // 1 nm
  constexpr double OpFastScintillation::Disk_SolidAngle(const double d, const double h,
                                                        const double b)
  {
    if(b <= 0. || d < 0. || h <= 0.) return 0.;
    const double leg2 = (b + d) * (b + d);
    const double aa = std::sqrt(h * h / (h * h + leg2));
    if(isApproximatelyZero(d)) {
      return 2. * CLHEP::pi * (1. - aa);
    }
    double bb = 2.*std::sqrt(b * d / (h * h + leg2));
    double cc = 4. * b * d / leg2;

    if(isDefinitelyGreaterThan(d,b)) {
      return 2.*aa*(std::sqrt(1.-cc)*boost::math::ellint_3(bb,cc,noLDoublePromote()) -
                    boost::math::ellint_1(bb,noLDoublePromote()));
    }
    if(isDefinitelyLessThan(d,b)) {
      return 2.* CLHEP::pi -
        2.*aa*(boost::math::ellint_1(bb, noLDoublePromote()) +
               std::sqrt(1.-cc)*boost::math::ellint_3(bb,cc,noLDoublePromote()));
    }
    if(isApproximatelyEqual(d,b)) {
      return CLHEP::pi - 2.*aa*boost::math::ellint_1(bb,noLDoublePromote());
    }
    return 0.;
  }


  // solid angle of rectangular aperture
  // TODO: what's up with all of that times 2, divided by 2?
  constexpr double OpFastScintillation::Rectangle_SolidAngle(const double a, const double b,
                                                             const double d)
  {
    double aa = a / (2.0 * d);
    double bb = b / (2.0 * d);
    double aux = (1 + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
    // return 4 * std::acos(std::sqrt(aux));
    return 4 * fast_acos(std::sqrt(aux));
  }


  // TODO: allow greater tolerance in comparisons, see note above on Disk_SolidAngle()
  constexpr double OpFastScintillation::Rectangle_SolidAngle(const dims o,
                                                             const std::array<double, 3> v)
  {
    // v is the position of the track segment with respect to
    // the center position of the arapuca window

    // arapuca plane fixed in x direction
    if(isApproximatelyZero(v[1]) &&
       isApproximatelyZero(v[2])) {
      return Rectangle_SolidAngle(o.h, o.w, v[0]);
    }
    if(isDefinitelyGreaterThan(v[1], o.h/2.0) &&
       isDefinitelyGreaterThan(v[2], o.w/2.0)) {
      double A = v[1] - o.h / 2.0;
      double B = v[2] - o.w / 2.0;
      double to_return = (Rectangle_SolidAngle(2 * (A + o.h), 2 * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2 * A, 2 * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2 * (A + o.h), 2 * B, v[0]) +
                          Rectangle_SolidAngle(2 * A, 2 * B, v[0])) / 4.0;
      return to_return;
    }
    if((v[1] <= o.h / 2.0) &&
       (v[2] <= o.w / 2.0)) {
      double A = -v[1] + o.h / 2.0;
      double B = -v[2] + o.w / 2.0;
      double to_return = (Rectangle_SolidAngle(2 * (o.h - A), 2 * (o.w - B), v[0]) +
                          Rectangle_SolidAngle(2 * A, 2 * (o.w - B), v[0]) +
                          Rectangle_SolidAngle(2 * (o.h - A), 2 * B, v[0]) +
                          Rectangle_SolidAngle(2 * A, 2 * B, v[0])) / 4.0;
      return to_return;
    }
    if(isDefinitelyGreaterThan(v[1], o.h/2.0) &&
       (v[2] <= o.w / 2.0)) {
      double A = v[1] - o.h / 2.0;
      double B = -v[2] + o.w / 2.0;
      double to_return = (Rectangle_SolidAngle(2 * (A + o.h), 2 * (o.w - B), v[0]) -
                          Rectangle_SolidAngle(2 * A, 2 * (o.w - B), v[0]) +
                          Rectangle_SolidAngle(2 * (A + o.h), 2 * B, v[0]) -
                          Rectangle_SolidAngle(2 * A, 2 * B, v[0])) / 4.0;
      return to_return;
    }
    if((v[1] <= o.h / 2.0) &&
       isDefinitelyGreaterThan(v[2], o.w/2.0)) {
      double A = -v[1] + o.h / 2.0;
      double B = v[2] - o.w / 2.0;
      double to_return = (Rectangle_SolidAngle(2 * (o.h - A), 2 * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2 * (o.h - A), 2 * B, v[0]) +
                          Rectangle_SolidAngle(2 * A, 2 * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2 * A, 2 * B, v[0])) / 4.0;
      return to_return;
    }
    // error message if none of these cases, i.e. something has gone wrong!
    // std::cout << "Warning: invalid solid angle call." << std::endl;
    return 0.0;
  }


  constexpr double acos_table(const double x)
  {
    if(x < 0. || x > 1.){
      std::cout << "Range out of bounds in acos_table, only defined in [0, 1]" << std::endl;
      exit(0);
    }
    return acos_arr[std::round(acos_bins*x)];
  }


  constexpr double fast_acos(double x) {
    double negate = double(x < 0);
    x = std::abs(x);
    x -= double(x>1.0)*(x-1.0); // <- equivalent to min(1.0,x), but faster
    double ret = -0.0187293;
    ret = ret * x;
    ret = ret + 0.0742610;
    ret = ret * x;
    ret = ret - 0.2121144;
    ret = ret * x;
    ret = ret + 1.5707288;
    ret = ret * std::sqrt(1.0-x);
    ret = ret - 2 * negate * ret;
    return negate * 3.14159265358979 + ret;
  }

}
