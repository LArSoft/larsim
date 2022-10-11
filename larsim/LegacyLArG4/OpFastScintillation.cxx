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
#include "Geant4/G4Poisson.hh"
#include "Geant4/G4Proton.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4ios.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/globals.hh"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"         // geo::vect::fillCoords()
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h" // geo::vect::toTVector3()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/LegacyLArG4/IonizationAndScintillation.h"
#include "larsim/LegacyLArG4/OpDetPhotonTable.h"
#include "larsim/LegacyLArG4/OpFastScintillation.hh"
#include "larsim/LegacyLArG4/ParticleListAction.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/LArG4Parameters.h"

// support libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TGeoSphere.h"
#include "TMath.h"
#include "TRandom3.h"

#include <cassert>
#include <cmath>
#include <limits>

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

// Define a new policy *not* internally promoting RealType to double:
typedef boost::math::policies::policy<
  // boost::math::policies::digits10<8>,
  boost::math::policies::promote_double<false>>
  noLDoublePromote;

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
    , fActiveVolumes{extractActiveVolumes(*(lar::providerFrom<geo::Geometry>()))}
    , bPropagate(!(art::ServiceHandle<sim::LArG4Parameters const>()->NoPhotonPropagation()))
    , fPVS(bPropagate ? art::ServiceHandle<phot::PhotonVisibilityService const>().get() : nullptr)
    , fUseNhitsModel(fPVS && fPVS->UseNhitsModel())
    // for now, limit to the active volume only if semi-analytic model is used
    , fOnlyActiveVolume(usesSemiAnalyticModel())
  {
    SetProcessSubType(25); // TODO: unhardcode
    fTrackSecondariesFirst = false;
    fFiniteRiseTime = false;
    YieldFactor = 1.0;
    ExcitationRatio = 1.0;

    const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();

    scintillationByParticleType = larp->ScintByParticleType();

    if (verboseLevel > 0) { G4cout << GetProcessName() << " is created " << G4endl; }

    BuildThePhysicsTable();
    emSaturation = NULL;

    if (bPropagate) {
      assert(fPVS);

      // Loading the position of each optical channel, neccessary for the parametrizatiuons of Nhits and prop-time
      geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

      {
        auto log = mf::LogTrace("OpFastScintillation")
                   << "OpFastScintillation: active volume boundaries from " << fActiveVolumes.size()
                   << " volumes:";
        for (auto const& [iCryo, box] : util::enumerate(fActiveVolumes)) {
          log << "\n - C:" << iCryo << ": " << box.Min() << " -- " << box.Max() << " cm";
        } // for
        log << "\n  (scintillation photons are propagated "
            << (fOnlyActiveVolume ? "only from active volumes" : "from anywhere") << ")";
      } // local scope

      if (usesSemiAnalyticModel() && (geom.Ncryostats() > 1U)) {
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
      } // if

      geo::Point_t const Cathode_centre{geom.TPC(0, 0).GetCathodeCenter().X(),
                                        fActiveVolumes[0].CenterY(),
                                        fActiveVolumes[0].CenterZ()};
      mf::LogTrace("OpFastScintillation") << "Cathode_centre: " << Cathode_centre << " cm";

      // std::cout << "\nInitialize acos_arr with " << acos_bins+1
      //           << " hence with a resolution of " << 1./acos_bins << std::endl;
      // for(size_t i=0; i<=acos_bins; ++i){
      //   acos_arr[i] = std::acos(i/double(acos_bins));
      // }

      for (size_t const i : util::counter(fPVS->NOpChannels())) {
        geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpDet(i);
        fOpDetCenter.push_back(opDet.GetCenter());

        if (dynamic_cast<TGeoSphere const*>(opDet.Shape()) != nullptr) { // sphere/dome
          fOpDetType.push_back(1);                                       // Dome PMTs
          fOpDetLength.push_back(-1);
          fOpDetHeight.push_back(-1);
        }
        else if (opDet.isBar()) {  // box
          fOpDetType.push_back(0); //Arapucas
          fOpDetLength.push_back(opDet.Length());
          fOpDetHeight.push_back(opDet.Height());
        }
        else {                     // disk
          fOpDetType.push_back(2); // Disk PMTs
          //    std::cout<<"Radio: "<<geom.OpDetGeoFromOpDet(i).RMax()<<std::endl;
          fOpDetLength.push_back(-1);
          fOpDetHeight.push_back(-1);
        }
      }

      if (fPVS->IncludePropTime()) {
        std::cout << "Using parameterisation of timings." << std::endl;
        //OLD VUV time parapetrization (to be removed soon)
        //fPVS->SetDirectLightPropFunctions(functions_vuv, fd_break, fd_max, ftf1_sampling_factor);
        //fPVS->SetReflectedCOLightPropFunctions(functions_vis, ft0_max, ft0_break_point);
        //New VUV time parapetrization
        fPVS->LoadTimingsForVUVPar(fparameters,
                                   fstep_size,
                                   fmax_d,
                                   fmin_d,
                                   fvuv_vgroup_mean,
                                   fvuv_vgroup_max,
                                   finflexion_point_distance,
                                   fangle_bin_timing_vuv);

        // create vector of empty TF1s that will be replaces with the parameterisations that are generated as they are required
        // default TF1() constructor gives function with 0 dimensions, can then check numDim to qucikly see if a parameterisation has been generated
        const size_t num_params =
          (fmax_d - fmin_d) /
          fstep_size; // for d < fmin_d, no parameterisaton, a delta function is used instead
        size_t num_angles = std::round(90 / fangle_bin_timing_vuv);
        VUV_timing = std::vector(num_angles, std::vector(num_params, TF1()));

        // initialise vectors to contain range parameterisations sampled to in each case
        // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling is regenerated, this is the slow part!
        VUV_max = std::vector(num_angles, std::vector(num_params, 0.0));
        VUV_min = std::vector(num_angles, std::vector(num_params, 0.0));

        // VIS time parameterisation
        if (fPVS->StoreReflected()) {
          // load parameters
          fPVS->LoadTimingsForVISPar(fdistances_refl,
                                     fradial_distances_refl,
                                     fcut_off_pars,
                                     ftau_pars,
                                     fvis_vmean,
                                     fangle_bin_timing_vis);
        }
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
        fPVS->LoadVUVSemiAnalyticProperties(
          fIsFlatPDCorr, fIsDomePDCorr, fdelta_angulo_vuv, fradius);
        if (!fIsFlatPDCorr && !fIsDomePDCorr) {
          throw cet::exception("OpFastScintillation")
            << "Both isFlatPDCorr and isDomePDCorr parameters are false, at least one type of "
               "parameterisation is required for the semi-analytic light simulation."
            << "\n";
        }
        if (fIsFlatPDCorr) {
          fPVS->LoadGHFlat(fGHvuvpars_flat, fborder_corr_angulo_flat, fborder_corr_flat);
        }
        if (fIsDomePDCorr) {
          fPVS->LoadGHDome(fGHvuvpars_dome, fborder_corr_angulo_dome, fborder_corr_dome);
        }
        // cathode center coordinates required for corrections
        fcathode_centre = geom.TPC(0, 0).GetCathodeCenter();
        fcathode_centre[1] = fActiveVolumes[0].CenterY();
        fcathode_centre[2] =
          fActiveVolumes[0].CenterZ(); // to get full cathode dimension rather than just single tpc

        if (fPVS->StoreReflected()) {
          fStoreReflected = true;
          // Load corrections for VIS semi-anlytic hits
          std::cout << "Loading visible light corrections" << std::endl;
          fPVS->LoadVisSemiAnalyticProperties(fdelta_angulo_vis, fradius);
          if (fIsFlatPDCorr) {
            fPVS->LoadVisParsFlat(fvis_distances_x_flat, fvis_distances_r_flat, fvispars_flat);
          }
          if (fIsDomePDCorr) {
            fPVS->LoadVisParsDome(fvis_distances_x_dome, fvis_distances_r_dome, fvispars_dome);
          }

          // cathode dimensions
          fcathode_ydimension = fActiveVolumes[0].SizeY();
          fcathode_zdimension = fActiveVolumes[0].SizeZ();

          // set cathode plane struct for solid angle function
          fcathode_plane.h = fcathode_ydimension;
          fcathode_plane.w = fcathode_zdimension;
          fplane_depth = std::abs(fcathode_centre[0]);
        }
        else
          fStoreReflected = false;
      }
    }
    tpbemission = lar::providerFrom<detinfo::LArPropertiesService>()->TpbEm();
    std::vector<double> parent;
    parent.reserve(tpbemission.size());
    for (auto iter = tpbemission.begin(); iter != tpbemission.end(); ++iter) {
      parent.push_back(iter->second);
    }
    fTPBEm = std::make_unique<CLHEP::RandGeneral>(parent.data(), parent.size());
  }

  ////////////////
  // Destructors
  ////////////////
  OpFastScintillation::~OpFastScintillation()
  {
    if (theFastIntegralTable) theFastIntegralTable->clearAndDestroy();
    if (theSlowIntegralTable) theSlowIntegralTable->clearAndDestroy();
  }

  ////////////
  // Methods
  ////////////

  // AtRestDoIt
  // ----------
  //
  G4VParticleChange* OpFastScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)
  // This routine simply calls the equivalent PostStepDoIt since all the
  // necessary information resides in aStep.GetTotalEnergyDeposit()
  {
    return OpFastScintillation::PostStepDoIt(aTrack, aStep);
  }

  // PostStepDoIt
  // -------------
  //
  G4VParticleChange* OpFastScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
  // This routine is called for each tracking step of a charged particle
  // in a scintillator. A Poisson/Gauss-distributed number of photons is
  // generated according to the scintillation yield formula, distributed
  // evenly along the track segment and uniformly into 4pi.
  {
    aParticleChange.Initialize(aTrack);
    // Check that we are in a material with a properties table, if not
    // just return
    const G4Material* aMaterial = aTrack.GetMaterial();
    G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
    if (!aMaterialPropertiesTable) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

    G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
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
    double MeanNumberOfPhotons =
      larg4::IonizationAndScintillation::Instance()->NumberScintillationPhotons();
    // double stepEnergy          = larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit()/CLHEP::MeV;
    RecordPhotonsProduced(aStep, MeanNumberOfPhotons); //, stepEnergy);
    if (verboseLevel > 0) {
      G4cout << "\n Exiting from OpFastScintillation::DoIt -- NumberOfSecondaries = "
             << aParticleChange.GetNumberOfSecondaries() << G4endl;
    }
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  void OpFastScintillation::ProcessStep(const G4Step& step)
  {
    if (step.GetTotalEnergyDeposit() <= 0) return;

    OpDetPhotonTable::Instance()->AddEnergyDeposit(
      -1,
      -1,
      1.0,                                                 //scintillation yield
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
      ParticleListAction::GetCurrentOrigTrackID(),
      step.GetPreStepPoint()->GetPhysicalVolume()->GetName());
  }

  bool OpFastScintillation::RecordPhotonsProduced(const G4Step& aStep,
                                                  double MeanNumberOfPhotons) //, double stepEnergy)
  {
    // make sure that whatever happens afterwards, the energy deposition is stored
    art::ServiceHandle<sim::LArG4Parameters const> lgp;
    if (lgp->FillSimEnergyDeposits()) ProcessStep(aStep);

    // Get the pointer to the fast scintillation table
    OpDetPhotonTable* fst = OpDetPhotonTable::Instance();

    const G4Track* aTrack = aStep.GetTrack();

    G4StepPoint const* pPreStepPoint = aStep.GetPreStepPoint();
    // unused G4StepPoint const* pPostStepPoint = aStep.GetPostStepPoint();

    const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
    const G4Material* aMaterial = aTrack->GetMaterial();

    G4int materialIndex = aMaterial->GetIndex();
    G4int tracknumber = aStep.GetTrack()->GetTrackID();

    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
    //G4double      t0 = pPreStepPoint->GetGlobalTime() - fGlobalTimeOffset;
    G4double t0 = pPreStepPoint->GetGlobalTime();

    G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();

    G4MaterialPropertyVector* Fast_Intensity =
      aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");
    G4MaterialPropertyVector* Slow_Intensity =
      aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

    if (!Fast_Intensity && !Slow_Intensity) return 1;

    if (!bPropagate) return 0;

    // Get the visibility vector for this point
    assert(fPVS);

    G4int nscnt = 1;
    if (Fast_Intensity && Slow_Intensity) nscnt = 2;

    double Num = 0;
    double YieldRatio = 0;

    if (scintillationByParticleType) {
      // The scintillation response is a function of the energy
      // deposited by particle types.

      // Get the definition of the current particle
      G4ParticleDefinition* pDef = aParticle->GetDefinition();

      // Obtain the G4MaterialPropertyVectory containing the
      // scintillation light yield as a function of the deposited
      // energy for the current particle type

      // Protons
      if (pDef == G4Proton::ProtonDefinition()) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("PROTONYIELDRATIO");
      }
      // Muons
      else if (pDef == G4MuonPlus::MuonPlusDefinition() ||
               pDef == G4MuonMinus::MuonMinusDefinition()) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("MUONYIELDRATIO");
      }
      // Pions
      else if (pDef == G4PionPlus::PionPlusDefinition() ||
               pDef == G4PionMinus::PionMinusDefinition()) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("PIONYIELDRATIO");
      }
      // Kaons
      else if (pDef == G4KaonPlus::KaonPlusDefinition() ||
               pDef == G4KaonMinus::KaonMinusDefinition()) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("KAONYIELDRATIO");
      }
      // Alphas
      else if (pDef == G4Alpha::AlphaDefinition()) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("ALPHAYIELDRATIO");
      }
      // Electrons (must also account for shell-binding energy
      // attributed to gamma from standard PhotoElectricEffect)
      else if (pDef == G4Electron::ElectronDefinition() || pDef == G4Gamma::GammaDefinition()) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("ELECTRONYIELDRATIO");
      }
      // Default for particles not enumerated/listed above
      else {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("ELECTRONYIELDRATIO");
      }
      // If the user has not specified yields for (p,d,t,a,carbon)
      // then these unspecified particles will default to the
      // electron's scintillation yield
      if (YieldRatio == 0) {
        YieldRatio = aMaterialPropertiesTable->GetConstProperty("ELECTRONYIELDRATIO");
      }
    }

    geo::Point_t const ScintPoint = {x0[0] / CLHEP::cm, x0[1] / CLHEP::cm, x0[2] / CLHEP::cm};
    if (fOnlyActiveVolume && !isScintInActiveVolume(ScintPoint)) return 0;
    const phot::MappedCounts_t& Visibilities = fPVS->GetAllVisibilities(ScintPoint);

    phot::MappedCounts_t ReflVisibilities;

    // Store timing information in the object for use in propagationTime method
    if (fPVS->StoreReflected()) {
      ReflVisibilities = fPVS->GetAllVisibilities(ScintPoint, true);
      if (fPVS->StoreReflT0()) ReflT0s = fPVS->GetReflT0s(ScintPoint);
    }
    if (fPVS->IncludeParPropTime()) { ParPropTimeTF1 = fPVS->GetTimingTF1(ScintPoint); }

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
    <<ScintPoint<<","<<t0<<","
    <<aStep.GetDeltaPosition().mag()<<","
    <<MeanNumberOfPhotons<<","<<std::flush;

    double gen_photon_ctr=0;
    double det_photon_ctr=0;
    */
    for (G4int scnt = 1; scnt <= nscnt; scnt++) {
      G4double ScintillationTime = 0. * CLHEP::ns;
      G4double ScintillationRiseTime = 0. * CLHEP::ns;
      G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;
      if (scnt == 1) {
        if (nscnt == 1) {
          if (Fast_Intensity) {
            ScintillationTime = aMaterialPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
            if (fFiniteRiseTime) {
              ScintillationRiseTime =
                aMaterialPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");
            }
            ScintillationIntegral =
              (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
          }
          if (Slow_Intensity) {
            ScintillationTime = aMaterialPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");
            if (fFiniteRiseTime) {
              ScintillationRiseTime =
                aMaterialPropertiesTable->GetConstProperty("SLOWSCINTILLATIONRISETIME");
            }
            ScintillationIntegral =
              (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
          }
        } //endif nscnt=1
        else {
          if (YieldRatio == 0)
            YieldRatio = aMaterialPropertiesTable->GetConstProperty("YIELDRATIO");

          if (ExcitationRatio == 1.0) { Num = std::min(YieldRatio, 1.0) * MeanNumberOfPhotons; }
          else {
            Num = std::min(ExcitationRatio, 1.0) * MeanNumberOfPhotons;
          }
          ScintillationTime = aMaterialPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
          if (fFiniteRiseTime) {
            ScintillationRiseTime =
              aMaterialPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");
          }
          ScintillationIntegral =
            (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
        } //endif nscnt!=1
      }   //end scnt=1

      else {
        Num = MeanNumberOfPhotons - Num;
        ScintillationTime = aMaterialPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");
        if (fFiniteRiseTime) {
          ScintillationRiseTime =
            aMaterialPropertiesTable->GetConstProperty("SLOWSCINTILLATIONRISETIME");
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
      //if (!Visibilities && (fPVS->NOpChannels() > 0)) {
      //  throw cet::exception("OpFastScintillator")
      //    << "Photon library does not cover point " << ScintPoint << " cm.\n";
      //}

      if (!Visibilities && !usesSemiAnalyticModel()) continue;

      // detected photons from direct light
      std::map<size_t, int> DetectedNum;
      if (Visibilities && !usesSemiAnalyticModel()) {
        for (size_t const OpDet : util::counter(fPVS->NOpChannels())) {
          if (fOpaqueCathode && !isOpDetInSameTPC(ScintPoint, fOpDetCenter.at(OpDet))) continue;
          int const DetThis = std::round(G4Poisson(Visibilities[OpDet] * Num));
          if (DetThis > 0) DetectedNum[OpDet] = DetThis;
        }
      }
      else {
        detectedDirectHits(DetectedNum, Num, ScintPoint);
      }

      // detected photons from reflected light
      std::map<size_t, int> ReflDetectedNum;
      if (fPVS->StoreReflected()) {
        if (!usesSemiAnalyticModel()) {
          for (size_t const OpDet : util::counter(fPVS->NOpChannels())) {
            if (fOpaqueCathode && !isOpDetInSameTPC(ScintPoint, fOpDetCenter.at(OpDet))) continue;
            int const ReflDetThis = std::round(G4Poisson(ReflVisibilities[OpDet] * Num));
            if (ReflDetThis > 0) ReflDetectedNum[OpDet] = ReflDetThis;
          }
        }
        else {
          detectedReflecHits(ReflDetectedNum, Num, ScintPoint);
        }
      }

      std::vector<double> arrival_time_dist;
      // Now we run through each PMT figuring out num of detected photons
      for (size_t Reflected = 0; Reflected <= 1; ++Reflected) {
        // Only do the reflected loop if we have reflected visibilities
        if (Reflected && !fPVS->StoreReflected()) continue;

        std::map<size_t, int>::const_iterator itstart;
        std::map<size_t, int>::const_iterator itend;
        if (Reflected) {
          itstart = ReflDetectedNum.begin();
          itend = ReflDetectedNum.end();
        }
        else {
          itstart = DetectedNum.begin();
          itend = DetectedNum.end();
        }
        for (auto itdetphot = itstart; itdetphot != itend; ++itdetphot) {
          const size_t OpChannel = itdetphot->first;
          const int NPhotons = itdetphot->second;

          // Set up the OpDetBTR information
          sim::OpDetBacktrackerRecord tmpOpDetBTRecord(OpChannel);
          int thisG4TrackID = ParticleListAction::GetCurrentTrackID();
          double xyzPos[3];
          average_position(aStep, xyzPos);
          double Edeposited = 0;
          if (scintillationByParticleType) {
            //We use this when it is the only sensical information. It may be of limited use to end users.
            Edeposited = aStep.GetTotalEnergyDeposit();
          }
          else if (emSaturation) {
            //If Birk Coefficient used, log VisibleEnergies.
            Edeposited =
              larg4::IonizationAndScintillation::Instance()->VisibleEnergyDeposit() / CLHEP::MeV;
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
            G4double Time = t0 + scint_time(aStep, ScintillationTime, ScintillationRiseTime) +
                            arrival_time_dist[i] * CLHEP::ns;

            // Always store the BTR
            tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, Time, 1, xyzPos, Edeposited);

            // Store as lite photon or as OnePhoton
            if (lgp->UseLitePhotons()) {
              fst->AddLitePhoton(OpChannel, static_cast<int>(Time), 1, Reflected);
            }
            else {
              // The sim photon in this case stores its production point and time
              TVector3 PhotonPosition(x0[0], x0[1], x0[2]);

              float PhotonEnergy = 0;
              if (Reflected)
                PhotonEnergy = reemission_energy() * CLHEP::eV;
              else
                PhotonEnergy = 9.7 * CLHEP::eV; // 9.7 eV peak of VUV emission spectrum

              // Make a photon object for the collection
              sim::OnePhoton PhotToAdd;
              PhotToAdd.InitialPosition = PhotonPosition;
              PhotToAdd.Energy = PhotonEnergy;
              PhotToAdd.Time = Time;
              PhotToAdd.SetInSD = false;
              PhotToAdd.MotherTrackID = tracknumber;

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

    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    G4int numOfMaterials = G4Material::GetNumberOfMaterials();

    // create new physics table
    if (!theFastIntegralTable)
      theFastIntegralTable = std::make_unique<G4PhysicsTable>(numOfMaterials);
    if (!theSlowIntegralTable)
      theSlowIntegralTable = std::make_unique<G4PhysicsTable>(numOfMaterials);

    // loop for materials
    for (G4int i = 0; i < numOfMaterials; i++) {
      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector = new G4PhysicsOrderedFreeVector();
      G4PhysicsOrderedFreeVector* bPhysicsOrderedFreeVector = new G4PhysicsOrderedFreeVector();

      // Retrieve vector of scintillation wavelength intensity for
      // the material from the material's optical properties table.
      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();

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

            aPhysicsOrderedFreeVector->InsertValues(currentPM, currentCII);

            // Set previous values to current ones prior to loop
            G4double prevPM = currentPM;
            G4double prevCII = currentCII;
            G4double prevIN = currentIN;

            // loop over all (photon energy, intensity)
            // pairs stored for this material
            for (size_t i = 1; i < theFastLightVector->GetVectorLength(); i++) {
              currentPM = theFastLightVector->Energy(i);
              currentIN = (*theFastLightVector)[i];

              currentCII = 0.5 * (prevIN + currentIN);

              currentCII = prevCII + (currentPM - prevPM) * currentCII;

              aPhysicsOrderedFreeVector->InsertValues(currentPM, currentCII);

              prevPM = currentPM;
              prevCII = currentCII;
              prevIN = currentIN;
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

            bPhysicsOrderedFreeVector->InsertValues(currentPM, currentCII);

            // Set previous values to current ones prior to loop
            G4double prevPM = currentPM;
            G4double prevCII = currentCII;
            G4double prevIN = currentIN;

            // loop over all (photon energy, intensity)
            // pairs stored for this material
            for (size_t i = 1; i < theSlowLightVector->GetVectorLength(); i++) {
              currentPM = theSlowLightVector->Energy(i);
              currentIN = (*theSlowLightVector)[i];

              currentCII = 0.5 * (prevIN + currentIN);

              currentCII = prevCII + (currentPM - prevPM) * currentCII;

              bPhysicsOrderedFreeVector->InsertValues(currentPM, currentCII);

              prevPM = currentPM;
              prevCII = currentCII;
              prevIN = currentIN;
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
      G4Exception("OpFastScintillation::SetScintillationByParticleType",
                  "Scint02",
                  JustWarning,
                  "Redefinition: Birks Saturation is replaced by ScintillationByParticleType!");
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

  G4double OpFastScintillation::GetMeanLifeTime(const G4Track&, G4ForceCondition* condition)
  {
    *condition = Forced;
    return DBL_MAX;
  }

  G4double OpFastScintillation::scint_time(const G4Step& aStep,
                                           G4double ScintillationTime,
                                           G4double ScintillationRiseTime) const
  {
    G4StepPoint const* pPreStepPoint = aStep.GetPreStepPoint();
    G4StepPoint const* pPostStepPoint = aStep.GetPostStepPoint();
    G4double avgVelocity = (pPreStepPoint->GetVelocity() + pPostStepPoint->GetVelocity()) / 2.;
    G4double deltaTime = aStep.GetStepLength() / avgVelocity;
    if (ScintillationRiseTime == 0.0) {
      deltaTime = deltaTime - ScintillationTime * std::log(G4UniformRand());
    }
    else {
      deltaTime = deltaTime + sample_time(ScintillationRiseTime, ScintillationTime);
    }
    return deltaTime;
  }

  void OpFastScintillation::propagationTime(std::vector<double>& arrival_time_dist,
                                            G4ThreeVector x0,
                                            const size_t OpChannel,
                                            bool Reflected) //const
  {
    assert(fPVS);
    if (fPVS->IncludeParPropTime() && fPVS->IncludePropTime()) {
      throw cet::exception("OpFastScintillation")
        << "Cannot have both propagation time models simultaneously.";
    }
    else if (fPVS->IncludeParPropTime() &&
             !(ParPropTimeTF1 && (ParPropTimeTF1[OpChannel].GetNdim() == 1))) {
      //Warning: TF1::GetNdim()==1 will tell us if the TF1 is really defined or it is the default one.
      //This will fix a segfault when using timing and interpolation.
      G4cout << "WARNING: Requested parameterized timing, but no function found. Not applying "
                "propagation time."
             << G4endl;
    }
    else if (fPVS->IncludeParPropTime()) {
      if (Reflected)
        throw cet::exception("OpFastScintillation")
          << "No parameterized propagation time for reflected light";
      for (size_t i = 0; i < arrival_time_dist.size(); ++i) {
        arrival_time_dist[i] = ParPropTimeTF1[OpChannel].GetRandom();
      }
    }
    else if (fPVS->IncludePropTime()) {
      // Get VUV photons arrival time distribution from the parametrization
      geo::Point_t const& opDetCenter = fOpDetCenter.at(OpChannel);
      if (!Reflected) {
        const G4ThreeVector OpDetPoint(
          opDetCenter.X() * CLHEP::cm, opDetCenter.Y() * CLHEP::cm, opDetCenter.Z() * CLHEP::cm);
        double distance_in_cm = (x0 - OpDetPoint).mag() / CLHEP::cm; // this must be in CENTIMETERS!
        double cosine = std::abs(x0[0] * CLHEP::cm - OpDetPoint[0] * CLHEP::cm) / distance_in_cm;
        double theta = fast_acos(cosine) * 180. / CLHEP::pi;
        int angle_bin = theta / fangle_bin_timing_vuv;
        getVUVTimes(arrival_time_dist, distance_in_cm, angle_bin); // in ns
      }
      else {
        TVector3 const ScintPoint(x0[0] / CLHEP::cm, x0[1] / CLHEP::cm, x0[2] / CLHEP::cm); // in cm
        getVISTimes(arrival_time_dist, ScintPoint, geo::vect::toTVector3(opDetCenter));     // in ns
      }
    }
  }

  G4double OpFastScintillation::sample_time(const G4double tau1, const G4double tau2) const
  {
    // tau1: rise time and tau2: decay time
    while (1) {
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
    return fTPBEm->fire() * ((*(--tpbemission.end())).first - (*tpbemission.begin()).first) +
           (*tpbemission.begin()).first;
  }

  void OpFastScintillation::average_position(G4Step const& aStep, double* xyzPos) const
  {
    CLHEP::Hep3Vector prePoint = (aStep.GetPreStepPoint())->GetPosition();
    CLHEP::Hep3Vector postPoint = (aStep.GetPostStepPoint())->GetPosition();
    xyzPos[0] = (((prePoint.getX() + postPoint.getX()) / 2.0) / CLHEP::cm);
    xyzPos[1] = (((prePoint.getY() + postPoint.getY()) / 2.0) / CLHEP::cm);
    xyzPos[2] = (((prePoint.getZ() + postPoint.getZ()) / 2.0) / CLHEP::cm);
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
  void OpFastScintillation::generateParam(const size_t index, const size_t angle_bin)
  {
    // get distance
    double distance_in_cm = (index * fstep_size) + fmin_d;

    // time range
    const double signal_t_range = 5000.; // TODO: unhardcode

    // parameterisation TF1
    TF1 fVUVTiming;

    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
    double t_direct_min = distance_in_cm / fvuv_vgroup_max;

    // Defining the model function(s) describing the photon transportation timing vs distance
    // Getting the landau parameters from the time parametrization
    std::array<double, 3> pars_landau;
    interpolate3(pars_landau,
                 fparameters[0][0],
                 fparameters[2][angle_bin],
                 fparameters[3][angle_bin],
                 fparameters[1][angle_bin],
                 distance_in_cm,
                 true);
    // Deciding which time model to use (depends on the distance)
    // defining useful times for the VUV arrival time shapes
    if (distance_in_cm >= finflexion_point_distance) {
      double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
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
      pars_expo[1] =
        interpolate(fparameters[4][0], fparameters[5][angle_bin], distance_in_cm, true);
      pars_expo[0] =
        interpolate(fparameters[4][0], fparameters[6][angle_bin], distance_in_cm, true);
      pars_expo[0] *= pars_landau[2];
      pars_expo[0] = std::log(pars_expo[0]);
      // this is to find the intersection point between the two functions:
      TF1 fint = TF1("fint", finter_d, pars_landau[0], 4 * t_direct_mean, 5);
      double parsInt[5] = {
        pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
      fint.SetParameters(parsInt);
      double t_int = fint.GetMinimumX();
      double minVal = fint.Eval(t_int);
      // the functions must intersect - output warning if they don't
      if (minVal > 0.015) {
        std::cout << "WARNING: Parametrization of VUV light discontinuous for distance = "
                  << distance_in_cm << std::endl;
        std::cout << "WARNING: This shouldn't be happening " << std::endl;
      }
      double parsfinal[7] = {t_int,
                             pars_landau[0],
                             pars_landau[1],
                             pars_landau[2],
                             pars_expo[0],
                             pars_expo[1],
                             t_direct_min};
      fVUVTiming.SetParameters(parsfinal);
    }

    // set the number of points used to sample parameterisation
    // for shorter distances, peak is sharper so more sensitive sampling required
    int fsampling; // TODO: unhardcode
    if (distance_in_cm < 50) { fsampling = 10000; }
    else if (distance_in_cm < 100) {
      fsampling = 5000;
    }
    else {
      fsampling = 1000;
    }
    fVUVTiming.SetNpx(fsampling);

    // calculate max and min distance relevant to sample parameterisation
    // max
    const size_t nq_max = 1;
    double xq_max[nq_max];
    double yq_max[nq_max];
    xq_max[0] = 0.995; // include 99.5%
    fVUVTiming.GetQuantiles(nq_max, yq_max, xq_max);
    double max = yq_max[0];
    // min
    double min = t_direct_min;

    // store TF1 and min/max, this allows identical TF1 to be used every time sampling
    // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
    // all subsequent calls check if it has been generated previously and are ~100+ times quicker
    VUV_timing[angle_bin][index] = fVUVTiming;
    VUV_max[angle_bin][index] = max;
    VUV_min[angle_bin][index] = min;
  }

  // VUV arrival times calculation function
  void OpFastScintillation::getVUVTimes(std::vector<double>& arrivalTimes,
                                        const double distance,
                                        const size_t angle_bin)
  {
    if (distance < fmin_d) {
      // times are fixed shift i.e. direct path only
      double t_prop_correction = distance / fvuv_vgroup_mean;
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = t_prop_correction;
      }
    }
    else { // distance >= fmin_d
      // determine nearest parameterisation in discretisation
      int index = std::round((distance - fmin_d) / fstep_size);
      // check whether required parameterisation has been generated, generating if not
      if (VUV_timing[angle_bin][index].GetNdim() == 0) { generateParam(index, angle_bin); }
      // randomly sample parameterisation for each photon
      for (size_t i = 0; i < arrivalTimes.size(); ++i) {
        arrivalTimes[i] = VUV_timing[angle_bin][index].GetRandom(VUV_min[angle_bin][index],
                                                                 VUV_max[angle_bin][index]);
      }
    }
  }

  // VIS arrival times calculation functions
  void OpFastScintillation::getVISTimes(std::vector<double>& arrivalTimes,
                                        const TVector3& ScintPoint,
                                        const TVector3& OpDetPoint)
  {
    // *************************************************************************************************
    //     Calculation of earliest arrival times and corresponding unsmeared distribution
    // *************************************************************************************************

    // set plane_depth for correct TPC:
    double plane_depth;
    if (ScintPoint[0] < 0) { plane_depth = -fplane_depth; }
    else {
      plane_depth = fplane_depth;
    }

    // calculate point of reflection for shortest path
    TVector3 bounce_point(plane_depth, ScintPoint[1], ScintPoint[2]);

    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point - ScintPoint).Mag();
    double Visdist = (OpDetPoint - bounce_point).Mag();

    // calculate times taken by VUV part of path
    int angle_bin_vuv = 0; // on-axis by definition
    getVUVTimes(arrivalTimes, VUVdist, angle_bin_vuv);

    // add visible direct path transport time
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] += Visdist / fvis_vmean;
    }

    // *************************************************************************************************
    //      Smearing of arrival time distribution
    // *************************************************************************************************
    // calculate fastest time possible
    // vis part
    double vis_time = Visdist / fvis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < fmin_d) { vuv_time = VUVdist / fvuv_vgroup_max; }
    else {
      // find index of required parameterisation
      const size_t index = std::round((VUVdist - fmin_d) / fstep_size);
      // find shortest time
      vuv_time = VUV_min[angle_bin_vuv][index];
    }
    // sum
    double fastest_time = vis_time + vuv_time;

    // calculate angle theta between bound_point and optical detector
    double cosine_theta = std::abs(OpDetPoint[0] - bounce_point[0]) / Visdist;
    double theta = fast_acos(cosine_theta) * 180. / CLHEP::pi;

    // determine smearing parameters using interpolation of generated points:
    // 1). tau = exponential smearing factor, varies with distance and angle
    // 2). cutoff = largest smeared time allowed, preventing excessively large
    //     times caused by exponential distance to cathode
    double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
    // angular bin
    size_t theta_bin = theta / fangle_bin_timing_vis;
    // radial distance from centre of TPC (y,z plane)
    double r = std::sqrt(std::pow(ScintPoint[1] - fcathode_centre[1], 2) +
                         std::pow(ScintPoint[2] - fcathode_centre[2], 2));

    // cut-off and tau
    // cut-off
    // interpolate in d_c for each r bin
    std::vector<double> interp_vals(fcut_off_pars[theta_bin].size(), 0.0);
    for (size_t i = 0; i < fcut_off_pars[theta_bin].size(); i++) {
      interp_vals[i] =
        interpolate(fdistances_refl, fcut_off_pars[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double cutoff = interpolate(fradial_distances_refl, interp_vals, r, true);

    // tau
    // interpolate in x for each r bin
    std::vector<double> interp_vals_tau(ftau_pars[theta_bin].size(), 0.0);
    for (size_t i = 0; i < ftau_pars[theta_bin].size(); i++) {
      interp_vals_tau[i] =
        interpolate(fdistances_refl, ftau_pars[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double tau = interpolate(fradial_distances_refl, interp_vals_tau, r, true);

    if (tau < 0) { tau = 0; } // failsafe if tau extrapolate goes wrong

    // apply smearing:
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      double arrival_time_smeared;
      // if time is already greater than cutoff, do not apply smearing
      if (arrivalTimes[i] >= cutoff) { continue; }
      // otherwise smear
      else {
        unsigned int counter = 0;
        // loop until time generated is within cutoff limit
        // most are within single attempt, very few take more than two
        do {
          // don't attempt smearings too many times
          if (counter >= 10) {                      // TODO: unhardcode
            arrival_time_smeared = arrivalTimes[i]; // don't smear
            break;
          }
          else {
            // generate random number in appropriate range
            double x = gRandom->Uniform(0.5, 1.0); // TODO: unhardcode
            // apply the exponential smearing
            arrival_time_smeared =
              arrivalTimes[i] + (arrivalTimes[i] - fastest_time) * (std::pow(x, -tau) - 1);
          }
          counter++;
        } while (arrival_time_smeared > cutoff);
      }
      arrivalTimes[i] = arrival_time_smeared;
    }
  }

  // ---------------------------------------------------------------------------
  bool OpFastScintillation::usesSemiAnalyticModel() const
  {
    return fUseNhitsModel;
  } // OpFastScintillation::usesSemiAnalyticModel()

  // ---------------------------------------------------------------------------
  void OpFastScintillation::detectedDirectHits(std::map<size_t, int>& DetectedNum,
                                               const double Num,
                                               geo::Point_t const& ScintPoint) const
  {
    for (size_t const OpDet : util::counter(fPVS->NOpChannels())) {
      if (!isOpDetInSameTPC(ScintPoint, fOpDetCenter.at(OpDet))) continue;

      // set detector struct for solid angle function
      const OpFastScintillation::OpticalDetector op{fOpDetHeight.at(OpDet),
                                                    fOpDetLength.at(OpDet),
                                                    fOpDetCenter.at(OpDet),
                                                    fOpDetType.at(OpDet)};
      const int DetThis = VUVHits(Num, ScintPoint, op);
      if (DetThis > 0) {
        DetectedNum[OpDet] = DetThis;
        //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
        //   //   it->second<<" " << Num << " " << DetThisPMT;
        //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
      }
    }
  }

  void OpFastScintillation::detectedReflecHits(std::map<size_t, int>& ReflDetectedNum,
                                               const double Num,
                                               geo::Point_t const& ScintPoint) const
  {
    // 1). calculate total number of hits of VUV photons on
    // reflective foils via solid angle + Gaisser-Hillas
    // corrections:

    // set plane_depth for correct TPC:
    double plane_depth;
    if (ScintPoint.X() < 0.) { plane_depth = -fplane_depth; }
    else {
      plane_depth = fplane_depth;
    }

    // get scintpoint coords relative to centre of cathode plane
    std::array<double, 3> ScintPoint_relative = {std::abs(ScintPoint.X() - plane_depth),
                                                 std::abs(ScintPoint.Y() - fcathode_centre[1]),
                                                 std::abs(ScintPoint.Z() - fcathode_centre[2])};
    // calculate solid angle of cathode from the scintillation point
    double solid_angle_cathode = Rectangle_SolidAngle(fcathode_plane, ScintPoint_relative);

    // calculate distance and angle between ScintPoint and hotspot
    // vast majority of hits in hotspot region directly infront of scintpoint,
    // therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
    double distance_cathode = std::abs(plane_depth - ScintPoint.X());
    // calculate hits on cathode plane via geometric acceptance
    double cathode_hits_geo = std::exp(-1. * distance_cathode / fL_abs_vuv) *
                              (solid_angle_cathode / (4. * CLHEP::pi)) * Num;
    // determine Gaisser-Hillas correction including border effects
    // use flat correction
    double r = std::sqrt(std::pow(ScintPoint.Y() - fcathode_centre[1], 2) +
                         std::pow(ScintPoint.Z() - fcathode_centre[2], 2));
    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    if (fIsFlatPDCorr) {
      pars_ini[0] = fGHvuvpars_flat[0][0];
      pars_ini[1] = fGHvuvpars_flat[1][0];
      pars_ini[2] = fGHvuvpars_flat[2][0];
      pars_ini[3] = fGHvuvpars_flat[3][0];
      s1 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[0], 0, true);
      s2 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[1], 0, true);
      s3 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[2], 0, true);
    }
    else
      std::cout
        << "Error: flat optical detector VUV correction required for reflected semi-analytic hits."
        << std::endl;

    // add border correction
    pars_ini[0] = pars_ini[0] + s1 * r;
    pars_ini[1] = pars_ini[1] + s2 * r;
    pars_ini[2] = pars_ini[2] + s3 * r;
    pars_ini[3] = pars_ini[3];

    // calculate corrected number of hits
    double GH_correction = Gaisser_Hillas(distance_cathode, pars_ini);
    const double cathode_hits_rec = GH_correction * cathode_hits_geo;

    // detemine hits on each PD
    const std::array<double, 3> hotspot = {plane_depth, ScintPoint.Y(), ScintPoint.Z()};
    for (size_t const OpDet : util::counter(fPVS->NOpChannels())) {
      if (!isOpDetInSameTPC(ScintPoint, fOpDetCenter.at(OpDet))) continue;

      // set detector struct for solid angle function
      const OpFastScintillation::OpticalDetector op{fOpDetHeight.at(OpDet),
                                                    fOpDetLength.at(OpDet),
                                                    fOpDetCenter.at(OpDet),
                                                    fOpDetType.at(OpDet)};

      int const ReflDetThis = VISHits(ScintPoint, op, cathode_hits_rec, hotspot);
      if (ReflDetThis > 0) { ReflDetectedNum[OpDet] = ReflDetThis; }
    }
  }

  // VUV semi-analytic hits calculation
  int OpFastScintillation::VUVHits(const double Nphotons_created,
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
    if (opDet.type == 0) {
      // get scintillation point coordinates relative to arapuca window centre
      std::array<double, 3> ScintPoint_rel = {std::abs(ScintPoint[0] - OpDetPoint[0]),
                                              std::abs(ScintPoint[1] - OpDetPoint[1]),
                                              std::abs(ScintPoint[2] - OpDetPoint[2])};
      // calculate solid angle
      solid_angle = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, ScintPoint_rel);
    }
    // PMTs (dome)
    else if (opDet.type == 1) {
      solid_angle = Omega_Dome_Model(distance, theta);
    }
    // PMTs (disk approximation)
    else if (opDet.type == 2) {
      // offset in z-y plane
      double d = dist(&ScintPoint[1], &OpDetPoint[1], 2);
      // drift distance (in x)
      double h = std::abs(ScintPoint[0] - OpDetPoint[0]);
      // Solid angle of a disk
      solid_angle = Disk_SolidAngle(d, h, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk"
                << std::endl;
    }

    // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
    double hits_geo =
      std::exp(-1. * distance / fL_abs_vuv) * (solid_angle / (4 * CLHEP::pi)) * Nphotons_created;

    // apply Gaisser-Hillas correction for Rayleigh scattering distance
    // and angular dependence offset angle bin
    const size_t j = (theta / fdelta_angulo_vuv);

    // determine GH parameters, accounting for border effects
    // radial distance from centre of detector (Y-Z)
    double r = dist(ScintPoint, fcathode_centre, 2, 1);
    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    // flat PDs
    if ((opDet.type == 0 || opDet.type == 2) && fIsFlatPDCorr) {
      pars_ini[0] = fGHvuvpars_flat[0][j];
      pars_ini[1] = fGHvuvpars_flat[1][j];
      pars_ini[2] = fGHvuvpars_flat[2][j];
      pars_ini[3] = fGHvuvpars_flat[3][j];
      s1 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[0], theta, true);
      s2 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[1], theta, true);
      s3 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[2], theta, true);
    }
    // dome PDs
    else if (opDet.type == 1 && fIsDomePDCorr) {
      pars_ini[0] = fGHvuvpars_dome[0][j];
      pars_ini[1] = fGHvuvpars_dome[1][j];
      pars_ini[2] = fGHvuvpars_dome[2][j];
      pars_ini[3] = fGHvuvpars_dome[3][j];
      s1 = interpolate(fborder_corr_angulo_dome, fborder_corr_dome[0], theta, true);
      s2 = interpolate(fborder_corr_angulo_dome, fborder_corr_dome[1], theta, true);
      s3 = interpolate(fborder_corr_angulo_dome, fborder_corr_dome[2], theta, true);
    }
    else
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk. Or "
                   "corrections for chosen optical detector type missing."
                << std::endl;

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

  // VIS hits semi-analytic model calculation
  int OpFastScintillation::VISHits(geo::Point_t const& ScintPoint_v,
                                   OpticalDetector const& opDet,
                                   const double cathode_hits_rec,
                                   const std::array<double, 3> hotspot) const
  {
    // 1). calculate total number of hits of VUV photons on reflective
    // foils via solid angle + Gaisser-Hillas corrections.
    // Done outside as it doesn't depend on OpDetPoint

    // set plane_depth for correct TPC:
    double plane_depth;
    if (ScintPoint_v.X() < 0.) { plane_depth = -fplane_depth; }
    else {
      plane_depth = fplane_depth;
    }

    // 2). calculate number of these hits which reach the optical
    // detector from the hotspot via solid angle

    // the interface has been converted into geo::Point_t, the implementation not yet
    std::array<double, 3U> ScintPoint;
    std::array<double, 3U> OpDetPoint;
    geo::vect::fillCoords(ScintPoint, ScintPoint_v);
    geo::vect::fillCoords(OpDetPoint, opDet.OpDetPoint);

    // calculate distances and angles for application of corrections
    // distance from hotspot to optical detector
    double distance_vis = dist(&hotspot[0], &OpDetPoint[0], 3);
    //  angle between hotspot and optical detector
    double cosine_vis = std::abs(hotspot[0] - OpDetPoint[0]) / distance_vis;
    // double theta_vis = std::acos(cosine_vis) * 180. / CLHEP::pi;
    double theta_vis = fast_acos(cosine_vis) * 180. / CLHEP::pi;

    // calculate solid angle of optical detector
    double solid_angle_detector = 0;
    // rectangular aperture
    if (opDet.type == 0) {
      // get hotspot coordinates relative to opDet
      std::array<double, 3> emission_relative = {std::abs(hotspot[0] - OpDetPoint[0]),
                                                 std::abs(hotspot[1] - OpDetPoint[1]),
                                                 std::abs(hotspot[2] - OpDetPoint[2])};
      // calculate solid angle
      solid_angle_detector = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, emission_relative);
    }
    // dome aperture
    else if (opDet.type == 1) {
      solid_angle_detector = Omega_Dome_Model(distance_vis, theta_vis);
    }
    // disk aperture
    else if (opDet.type == 2) {
      // offset in z-y plane
      double d = dist(&hotspot[1], &OpDetPoint[1], 2);
      // drift distance (in x)
      double h = std::abs(hotspot[0] - OpDetPoint[0]);
      // calculate solid angle
      solid_angle_detector = Disk_SolidAngle(d, h, fradius);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk"
                << std::endl;
    }

    // calculate number of hits via geometeric acceptance
    double hits_geo = (solid_angle_detector / (2. * CLHEP::pi)) *
                      cathode_hits_rec; // 2*pi due to presence of reflective foils

    // determine correction factor, depending on PD type
    const size_t k = (theta_vis / fdelta_angulo_vis); // off-set angle bin
    double r =
      dist(ScintPoint, fcathode_centre, 2, 1); // radial distance from centre of detector (Y-Z)
    double d_c = std::abs(ScintPoint[0] - plane_depth); // distance to cathode
    double border_correction = 0;
    // flat PDs
    if ((opDet.type == 0 || opDet.type == 2) && fIsFlatPDCorr) {
      // interpolate in d_c for each r bin
      const size_t nbins_r = fvispars_flat[k].size();
      std::vector<double> interp_vals(nbins_r, 0.0);
      {
        size_t idx = 0;
        size_t size = fvis_distances_x_flat.size();
        if (d_c >= fvis_distances_x_flat[size - 2])
          idx = size - 2;
        else {
          while (d_c > fvis_distances_x_flat[idx + 1])
            idx++;
        }
        for (size_t i = 0; i < nbins_r; ++i) {
          interp_vals[i] = interpolate(fvis_distances_x_flat, fvispars_flat[k][i], d_c, false, idx);
        }
      }
      // interpolate in r
      border_correction = interpolate(fvis_distances_r_flat, interp_vals, r, false);
    }
    // dome PDs
    else if (opDet.type == 1 && fIsDomePDCorr) {
      // interpolate in d_c for each r bin
      const size_t nbins_r = fvispars_dome[k].size();
      std::vector<double> interp_vals(nbins_r, 0.0);
      {
        size_t idx = 0;
        size_t size = fvis_distances_x_dome.size();
        if (d_c >= fvis_distances_x_dome[size - 2])
          idx = size - 2;
        else {
          while (d_c > fvis_distances_x_dome[idx + 1])
            idx++;
        }
        for (size_t i = 0; i < nbins_r; ++i) {
          interp_vals[i] = interpolate(fvis_distances_x_dome, fvispars_dome[k][i], d_c, false, idx);
        }
      }
      // interpolate in r
      border_correction = interpolate(fvis_distances_r_dome, interp_vals, r, false);
    }
    else {
      std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = dome, 2 = disk. Or "
                   "corrections for chosen optical detector type missing."
                << std::endl;
    }

    // apply correction factor
    double hits_rec = border_correction * hits_geo / cosine_vis;
    double hits_vis = std::round(G4Poisson(hits_rec));

    return hits_vis;
  }

  bool OpFastScintillation::isOpDetInSameTPC(geo::Point_t const& ScintPoint,
                                             geo::Point_t const& OpDetPoint) const
  {
    // check optical channel is in same TPC as scintillation light, if not return 0 hits
    // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
    // check x coordinate has same sign or is close to zero, otherwise return 0 hits
    if (((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) &&
        std::abs(OpDetPoint.X()) > 10.) { // TODO: unhardcode
      return false;
    }
    return true;
  }

  bool OpFastScintillation::isScintInActiveVolume(geo::Point_t const& ScintPoint)
  {
    //semi-analytic approach only works in the active volume
    return fActiveVolumes[0].ContainsPosition(ScintPoint);
  }

  G4double OpFastScintillation::single_exp(const G4double t, const G4double tau2) const
  {
    return std::exp(-1.0 * t / tau2) / tau2;
  }

  G4double OpFastScintillation::bi_exp(const G4double t,
                                       const G4double tau1,
                                       const G4double tau2) const
  { // TODO: what's up with this? ... / tau2 / tau2 ...
    return std::exp(-1.0 * t / tau2) * (1 - std::exp(-1.0 * t / tau1)) / tau2 / tau2 *
           (tau1 + tau2);
  }

  G4double OpFastScintillation::Gaisser_Hillas(const double x, const double* par) const
  {
    double X_mu_0 = par[3];
    double Normalization = par[0];
    double Diff = par[1] - X_mu_0;
    double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
    double Exponential = std::exp((par[1] - x) / par[2]);
    return (Normalization * Term * Exponential);
  }

  double finter_d(double* x, double* par)
  {
    double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
    double y2 = TMath::Exp(par[3] + x[0] * par[4]);
    return TMath::Abs(y1 - y2);
  }

  double LandauPlusExpoFinal(double* x, double* par)
  {
    // par0 = joining point
    // par1 = Landau MPV
    // par2 = Landau widt
    // par3 = normalization
    // par4 = Expo cte
    // par5 = Expo tau
    double y1 = par[3] * TMath::Landau(x[0], par[1], par[2]);
    double y2 = TMath::Exp(par[4] + x[0] * par[5]);
    if (x[0] > par[0]) y1 = 0.;
    if (x[0] < par[0]) y2 = 0.;
    return (y1 + y2);
  }

  double finter_r(double* x, double* par)
  {
    double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
    double y2 = par[5] * TMath::Landau(x[0], par[3], par[4]);
    return TMath::Abs(y1 - y2);
  }

  double model_close(double* x, double* par)
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
    if (x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
    if (x[0] < par[0]) y2 = 0.;
    return (y1 + y2);
  }

  double model_far(double* x, double* par)
  {
    // par1 = Landau MPV
    // par2 = Landau width
    // par3 = normalization
    // par0 = t_min
    double y = par[3] * TMath::Landau(x[0], par[1], par[2]);
    if (x[0] <= par[0]) y = 0.;
    return y;
  }

  //======================================================================
  //   Returns interpolated value at x from parallel arrays ( xData, yData )
  //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
  //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
  double OpFastScintillation::interpolate(const std::vector<double>& xData,
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

  void OpFastScintillation::interpolate3(std::array<double, 3>& inter,
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

  // solid angle of circular aperture
  // TODO: allow greater tolerance in comparisons, by default its using:
  // std::numeric_limits<double>::epsilon(): 2.22045e-16
  // that's an unrealistic small number, better setting
  // constexpr double tolerance = 0.0000001; // 1 nm
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
          mf::LogWarning("OpFastScintillation")
            << "Elliptic Integral in Disk_SolidAngle() given parameters outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("OpFastScintillation")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters outside domain.\n"
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
          mf::LogWarning("OpFastScintillation")
            << "Elliptic Integral in Disk_SolidAngle() given parameters outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("OpFastScintillation")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters outside domain.\n"
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

  // solid angle of rectangular aperture
  double OpFastScintillation::Rectangle_SolidAngle(const double a,
                                                   const double b,
                                                   const double d) const
  {
    double aa = a / (2. * d);
    double bb = b / (2. * d);
    double aux = (1. + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
    // return 4 * std::acos(std::sqrt(aux));
    return 4. * fast_acos(std::sqrt(aux));
  }

  // TODO: allow greater tolerance in comparisons, see note above on Disk_SolidAngle()
  double OpFastScintillation::Rectangle_SolidAngle(Dims const& o,
                                                   const std::array<double, 3> v) const
  {
    // v is the position of the track segment with respect to
    // the center position of the arapuca window

    // arapuca plane fixed in x direction
    if (isApproximatelyZero(v[1]) && isApproximatelyZero(v[2])) {
      return Rectangle_SolidAngle(o.h, o.w, v[0]);
    }
    if (isDefinitelyGreaterThan(v[1], o.h * .5) && isDefinitelyGreaterThan(v[2], o.w * .5)) {
      double A = v[1] - o.h * .5;
      double B = v[2] - o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (A + o.h), 2. * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2. * A, 2. * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2. * (A + o.h), 2. * B, v[0]) +
                          Rectangle_SolidAngle(2. * A, 2. * B, v[0])) *
                         .25;
      return to_return;
    }
    if ((v[1] <= o.h * .5) && (v[2] <= o.w * .5)) {
      double A = -v[1] + o.h * .5;
      double B = -v[2] + o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (o.h - A), 2. * (o.w - B), v[0]) +
                          Rectangle_SolidAngle(2. * A, 2. * (o.w - B), v[0]) +
                          Rectangle_SolidAngle(2. * (o.h - A), 2. * B, v[0]) +
                          Rectangle_SolidAngle(2. * A, 2. * B, v[0])) *
                         .25;
      return to_return;
    }
    if (isDefinitelyGreaterThan(v[1], o.h * .5) && (v[2] <= o.w * .5)) {
      double A = v[1] - o.h * .5;
      double B = -v[2] + o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (A + o.h), 2. * (o.w - B), v[0]) -
                          Rectangle_SolidAngle(2. * A, 2. * (o.w - B), v[0]) +
                          Rectangle_SolidAngle(2. * (A + o.h), 2. * B, v[0]) -
                          Rectangle_SolidAngle(2. * A, 2. * B, v[0])) *
                         .25;
      return to_return;
    }
    if ((v[1] <= o.h * .5) && isDefinitelyGreaterThan(v[2], o.w * .5)) {
      double A = -v[1] + o.h * .5;
      double B = v[2] - o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (o.h - A), 2. * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2. * (o.h - A), 2. * B, v[0]) +
                          Rectangle_SolidAngle(2. * A, 2. * (B + o.w), v[0]) -
                          Rectangle_SolidAngle(2. * A, 2. * B, v[0])) *
                         .25;
      return to_return;
    }
    // error message if none of these cases, i.e. something has gone wrong!
    // std::cout << "Warning: invalid solid angle call." << std::endl;
    return 0.;
  }

  double OpFastScintillation::Omega_Dome_Model(const double distance, const double theta) const
  {
    // this function calculates the solid angle of a semi-sphere of radius b,
    // as a correction to the analytic formula of the on-axix solid angle,
    // as we move off-axis an angle theta. We have used 9-angular bins
    // with delta_theta width.

    // par0 = Radius correction close
    // par1 = Radius correction far
    // par2 = breaking distance betwween "close" and "far"

    double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221};
    double par1[9] = {0, 0, 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
    const double delta_theta = 10.;
    int j = int(theta / delta_theta);
    // PMT radius
    const double b = fradius; // cm
    // distance form which the model parameters break (empirical value)
    const double d_break = 5 * b; //par2

    if (distance >= d_break) {
      double R_apparent_far = b - par1[j];
      return (2 * CLHEP::pi * (1 - std::sqrt(1 - std::pow(R_apparent_far / distance, 2))));
    }
    else {
      double R_apparent_close = b - par0[j];
      return (2 * CLHEP::pi * (1 - std::sqrt(1 - std::pow(R_apparent_close / distance, 2))));
    }
  }

  // ---------------------------------------------------------------------------
  std::vector<geo::BoxBoundedGeo> OpFastScintillation::extractActiveVolumes(
    geo::GeometryCore const& geom)
  {
    std::vector<geo::BoxBoundedGeo> activeVolumes;
    activeVolumes.reserve(geom.Ncryostats());

    for (geo::CryostatGeo const& cryo : geom.IterateCryostats()) {

      // can't use it default-constructed since it would always include origin
      geo::BoxBoundedGeo box{cryo.TPC(0).ActiveBoundingBox()};

      for (geo::TPCGeo const& TPC : cryo.IterateTPCs())
        box.ExtendToInclude(TPC.ActiveBoundingBox());

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
