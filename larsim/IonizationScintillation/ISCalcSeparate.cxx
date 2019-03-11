////////////////////////////////////////////////////////////////////////
/// \file  ISCalcuationSeparate.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using separate algorithms for each
///
/// Wes, 18Feb2018: this is a copy of the original, for standalone purposes
///
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include "CLHEP/Vector/ThreeVector.h"

#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceCharge/SpaceCharge.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"


// Wes did this for TRACE debugging
//#include "trace.h"
//#define TNAME "ISCalcSeparate"


namespace larg4{

  //----------------------------------------------------------------------------
  ISCalcSeparate::ISCalcSeparate()
  {
  }

  //----------------------------------------------------------------------------
  ISCalcSeparate::~ISCalcSeparate()
  {
  }

  //----------------------------------------------------------------------------
  void ISCalcSeparate::Initialize(const detinfo::LArProperties* larp,   
					 const detinfo::DetectorProperties* detp,
					 const sim::LArG4Parameters* lgp,
					 const spacecharge::SpaceCharge* sce)
  {
    //TRACEN(TNAME,3,"Initializing called.");

    fLArProp = larp;
    fSCE = sce;
    fDetProp = detp;
    fLArG4Prop = lgp;
    
    // \todo get scintillation yield from LArG4Parameters or LArProperties
    fScintYieldFactor  = 1.;

    // the recombination coefficient is in g/(MeVcm^2), but
    // we report energy depositions in MeV/cm, need to divide
    // Recombk from the LArG4Parameters service by the density
    // of the argon we got above.
    fRecombA             = lgp->RecombA();
    fRecombk             = lgp->Recombk()/detp->Density(detp->Temperature());
    fModBoxA             = lgp->ModBoxA();
    fModBoxB             = lgp->ModBoxB()/detp->Density(detp->Temperature());
    fUseModBoxRecomb     = (bool)lgp->UseModBoxRecomb();  
    fGeVToElectrons      = lgp->GeVToElectrons();

    this->Reset();
    
    //TRACEN(TNAME,3,"Initialize: RecombA=%f, Recombk=%f, ModBoxA=%f, ModBoxB=%f, UseModBoxRecomb=%d",
    //	   fRecombA, fRecombk, fModBoxA, fModBoxB, fUseModBoxRecomb);
    
    return;
  }

  //----------------------------------------------------------------------------
  void ISCalcSeparate::Reset()
  {
    fEnergyDeposit   = 0.;
    fNumScintPhotons = 0.;
    fNumIonElectrons = 0.;

    return;
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  void ISCalcSeparate::CalculateIonization(float e, float ds,
					   float x, float y, float z){

    //TRACEN(TNAME,3,"Calculate ionization called with e=%f, ds=%f, (x,y,z)=(%f,%f,%f)",
    //       e,ds,x,y,z);
    //TRACEN(TNAME,3,"Initialized values: RecombA=%lf, Recombk=%lf, ModBoxA=%lf, ModBoxB=%lf, UseModBoxRecomb=%d",
    //	   fRecombA, fRecombk, fModBoxA, fModBoxB, fUseModBoxRecomb);

    double recomb = 0.;
    double dEdx   = (ds<=0.0)? 0.0: e/ds;
    double EFieldStep = EFieldAtStep(fDetProp->Efield(),x,y,z);
    
    // Guard against spurious values of dE/dx. Note: assumes density of LAr
    if(dEdx < 1.) dEdx = 1.;
    
    //TRACEN(TNAME,3,"dEdx=%f, EFieldStep=%f",dEdx,EFieldStep);
    
    if(fUseModBoxRecomb) {
      if(ds>0){
	double Xi = fModBoxB * dEdx / EFieldStep;
	recomb = log(fModBoxA + Xi) / Xi;
      }
      else 
	recomb = 0;
    }
    else{
      recomb = fRecombA / (1. + dEdx * fRecombk / EFieldStep);
    }
    
    //TRACEN(TNAME,3,"recomb=%f",recomb);
    
    // 1.e-3 converts fEnergyDeposit to GeV
    fNumIonElectrons = fGeVToElectrons * 1.e-3 * e * recomb;
    
    //TRACEN(TNAME,3,"n_electrons=%f",fNumIonElectrons);
    
    MF_LOG_DEBUG("ISCalcSeparate") 
      << " Electrons produced for " << fEnergyDeposit 
      << " MeV deposited with "     << recomb 
      << " recombination: "         << fNumIonElectrons << std::endl; 
  }


  //----------------------------------------------------------------------------
  void ISCalcSeparate::CalculateIonization(sim::SimEnergyDeposit const& edep){
    CalculateIonization(edep.Energy(),edep.StepLength(),
			edep.MidPointX(),edep.MidPointY(),edep.MidPointZ());
  }
  
  //----------------------------------------------------------------------------
  void ISCalcSeparate::CalculateScintillation(float e, int pdg)
  {
    double scintYield = fLArProp->ScintYield(true);
    if(fLArProp->ScintByParticleType()){

      MF_LOG_DEBUG("ISCalcSeparate") << "scintillating by particle type";

      switch(pdg) {

      case 2212:
	scintYield = fLArProp->ProtonScintYield(true);
	break;
      case 13:
      case -13:
	scintYield = fLArProp->MuonScintYield(true);
	break;
      case 211:
      case -211:
	scintYield = fLArProp->PionScintYield(true);
	break;
      case 321:
      case -321:
	scintYield = fLArProp->KaonScintYield(true);
	break;
      case 1000020040:
	scintYield = fLArProp->AlphaScintYield(true);
	break;
      case 11:
      case -11:
      case 22:
	scintYield = fLArProp->ElectronScintYield(true);
	break;
      default:
	scintYield = fLArProp->ElectronScintYield(true);

      }

      fNumScintPhotons =  scintYield * e;
    }
    else
      fNumScintPhotons = fScintYieldFactor * scintYield * e;
    
  }

  //----------------------------------------------------------------------------
  void ISCalcSeparate::CalculateScintillation(sim::SimEnergyDeposit const& edep)
  {
    CalculateScintillation(edep.Energy(),edep.PdgCode());
  }

  //----------------------------------------------------------------------------
  void ISCalcSeparate::CalculateIonizationAndScintillation(sim::SimEnergyDeposit const& edep)
  {
    fEnergyDeposit = edep.Energy();
    CalculateIonization(edep);
    CalculateScintillation(edep);
  }

  double ISCalcSeparate::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
  {
    return EFieldAtStep(efield,
			edep.MidPointX(),edep.MidPointY(),edep.MidPointZ());
  }
  
  double ISCalcSeparate::EFieldAtStep(double efield, float x, float y, float z)
  {
    double EField = efield;
    if (fSCE->EnableSimEfieldSCE())
      {
        fEfieldOffsets = fSCE->GetEfieldOffsets(geo::Point_t{x,y,z});
        EField = std::sqrt( (efield + efield*fEfieldOffsets.X())*(efield + efield*fEfieldOffsets.X()) +
			    (efield*fEfieldOffsets.Y()*efield*fEfieldOffsets.Y()) +
			    (efield*fEfieldOffsets.Z()*efield*fEfieldOffsets.Z()) );
      }
    return EField;
  }

}// namespace
