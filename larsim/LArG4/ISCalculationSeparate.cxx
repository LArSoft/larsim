////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using separate algorithms for each
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include "CLHEP/Vector/ThreeVector.h"

#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4EmSaturation.hh"

#include "larsim/LArG4/ISCalculationSeparate.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Simulation/LArVoxelCalculator.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"


namespace larg4{

  //----------------------------------------------------------------------------
  ISCalculationSeparate::ISCalculationSeparate(CLHEP::HepRandomEngine&)
  {
  }

  //----------------------------------------------------------------------------
  ISCalculationSeparate::~ISCalculationSeparate()
  {
  }

  //----------------------------------------------------------------------------
  void ISCalculationSeparate::Initialize()
  {
    art::ServiceHandle<sim::LArG4Parameters> lgpHandle;
    const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    double density       = detprop->Density(detprop->Temperature());
    fEfield       = detprop->Efield();
    fScintByParticleType = larp->ScintByParticleType();
    fGeVToElectrons      = lgpHandle->GeVToElectrons();
    
    // \todo get scintillation yield from LArG4Parameters or LArProperties
    fScintYieldFactor  = 1.;

    // the recombination coefficient is in g/(MeVcm^2), but
    // we report energy depositions in MeV/cm, need to divide
    // Recombk from the LArG4Parameters service by the density
    // of the argon we got above.
    fRecombA             = lgpHandle->RecombA();
    fRecombk             = lgpHandle->Recombk()/density;
    fModBoxA             = lgpHandle->ModBoxA();
    fModBoxB             = lgpHandle->ModBoxB()/density;
    fUseModBoxRecomb     = lgpHandle->UseModBoxRecomb();  

    // Use Birks Correction in the Scintillation process    
    fEMSaturation = G4LossTableManager::Instance()->EmSaturation();

    // determine the step size using the voxel sizes
    art::ServiceHandle<sim::LArVoxelCalculator> lvc;
    double maxsize = std::max(lvc->VoxelSizeX(), std::max(lvc->VoxelSizeY(), lvc->VoxelSizeZ())) * CLHEP::cm;

    fStepSize = 0.1 * maxsize;

    return;
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  void ISCalculationSeparate::Reset()
  {
    fEnergyDeposit   = 0.;
    fNumScintPhotons = 0.;
    fNumIonElectrons = 0.;

    return;
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  void ISCalculationSeparate::CalculateIonizationAndScintillation(const G4Step* step)
  {

    fEnergyDeposit = step->GetTotalEnergyDeposit()/CLHEP::MeV;

    // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
    // R = A/(1 + (dE/dx)*k)
    // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
    // from GeV/voxel width
    // A = 0.800 +/- 0.003
    // k = (0.097+/-0.001) g/(MeVcm^2)
    // the dx depends on the trajectory of the step
    // k should be divided by the density as our dE/dx is in MeV/cm,
    // the division is handled in the constructor when we set fRecombk
    // B.Baller: Add Modified Box recombination - ArgoNeuT result submitted to JINST

    G4ThreeVector totstep = step->GetPostStepPoint()->GetPosition();
    totstep -= step->GetPreStepPoint()->GetPosition();

    double dx     = totstep.mag()/CLHEP::cm;
    double recomb = 0.;
    double dEdx   = (dx == 0.0)? 0.0: fEnergyDeposit/dx;
    double EFieldStep = EFieldAtStep(fEfield,step);

    // Guard against spurious values of dE/dx. Note: assumes density of LAr
    if(dEdx < 1.) dEdx = 1.;

    if(fUseModBoxRecomb) {
      if (dx){
	double Xi = fModBoxB * dEdx / EFieldStep;
	recomb = log(fModBoxA + Xi) / Xi;
      }
      else 
	recomb = 0;
    } 
    else{
      recomb = fRecombA/(1. + dEdx * fRecombk / EFieldStep);
    }


    // 1.e-3 converts fEnergyDeposit to GeV
    fNumIonElectrons = fGeVToElectrons * 1.e-3 * fEnergyDeposit * recomb;

    MF_LOG_DEBUG("ISCalculationSeparate") << " Electrons produced for " << fEnergyDeposit 
				       << " MeV deposited with "     << recomb 
				       << " recombination: "         << fNumIonElectrons;

    // Now do the scintillation
    G4MaterialPropertiesTable* mpt = step->GetTrack()->GetMaterial()->GetMaterialPropertiesTable();
    if( !mpt) 
      throw cet::exception("ISCalculationSeparate") << "Cannot find materials property table"
						    << " for this step! "
						    << step->GetTrack()->GetMaterial() << "\n";

    // if not doing the scintillation by particle type, use the saturation
    double scintYield = mpt->GetConstProperty("SCINTILLATIONYIELD");

    if(fScintByParticleType){

      MF_LOG_DEBUG("ISCalculationSeparate") << "scintillating by particle type";

      // Get the definition of the current particle
      G4ParticleDefinition *pDef = step->GetTrack()->GetDynamicParticle()->GetDefinition();
      //G4MaterialPropertyVector *Scint_Yield_Vector = NULL;

      // Obtain the G4MaterialPropertyVectory containing the
      // scintillation light yield as a function of the deposited
      // energy for the current particle type
      
      // Protons
      if(pDef == G4Proton::ProtonDefinition()){
	scintYield = mpt->GetConstProperty("PROTONSCINTILLATIONYIELD");
      }
      // Muons
      else if(pDef == G4MuonPlus::MuonPlusDefinition() ||
	      pDef == G4MuonMinus::MuonMinusDefinition()){
	scintYield = mpt->GetConstProperty("MUONSCINTILLATIONYIELD");
      }
      // Pions
      else if(pDef == G4PionPlus::PionPlusDefinition() ||
	      pDef == G4PionMinus::PionMinusDefinition()){
	scintYield = mpt->GetConstProperty("PIONSCINTILLATIONYIELD");
      }
      // Kaons
      else if(pDef == G4KaonPlus::KaonPlusDefinition() ||
	      pDef == G4KaonMinus::KaonMinusDefinition()){
	scintYield = mpt->GetConstProperty("KAONSCINTILLATIONYIELD");
      }
      // Alphas
      else if(pDef == G4Alpha::AlphaDefinition()){
	scintYield = mpt->GetConstProperty("ALPHASCINTILLATIONYIELD");
      }
      // Electrons (must also account for shell-binding energy
      // attributed to gamma from standard PhotoElectricEffect)
      else if(pDef == G4Electron::ElectronDefinition() ||
	      pDef == G4Gamma::GammaDefinition()){
	scintYield = mpt->GetConstProperty("ELECTRONSCINTILLATIONYIELD");
      }       
      // Default for particles not enumerated/listed above
      else{
	scintYield = mpt->GetConstProperty("ELECTRONSCINTILLATIONYIELD");
      }

      // If the user has not specified yields for (p,d,t,a,carbon)
      // then these unspecified particles will default to the 
      // electron's scintillation yield
	   
      // Throw an exception if no scintillation yield is found
      if (!scintYield) 
	throw cet::exception("ISCalculationSeparate") << "Request for scintillation yield for energy "
						      << "deposit and particle type without correct "
						      << "entry in MaterialPropertiesTable\n"
						      << "ScintillationByParticleType requires at "
						      << "minimum that ELECTRONSCINTILLATIONYIELD is "
						      << "set by the user\n";

      fNumScintPhotons =  scintYield * fEnergyDeposit;
    }
    else if(fEMSaturation){
      // The default linear scintillation process
      //fEMSaturation->SetVerbose(1);
      fVisibleEnergyDeposition = fEMSaturation->VisibleEnergyDepositionAtAStep(step);
      fNumScintPhotons = fScintYieldFactor * scintYield * fVisibleEnergyDeposition;
      //fNumScintPhotons = fScintYieldFactor * scintYield * fEMSaturation->VisibleEnergyDepositionAtAStep(step);
      //I need a dump here
      //mf::LogInfo("EMSaturation") <<"\n\nfEMSaturation VisibleEnergyDepositionAtAStep(step): "<<fEMSaturation->VisibleEnergyDepositionAtAStep(step)<<"\n" <<"BirksCoefs: \n";
     // fEMSaturation->DumpBirksCoefficients();
      //mf::LogInfo("EMSaturation")<<"\n" <<"G4Birks: \n";
      //fEMSaturation->DumpG4BirksCoefficients();
      //mf::LogInfo("EMSaturation")<<"fScintYieldFactor: "<<fScintYieldFactor <<"\nscintYield: "<<scintYield<<"\nfEMVisAtStep: "<<fEMSaturation->VisibleEnergyDepositionAtAStep(step)<<"\nfNumScintPhotons: "<<fNumScintPhotons<<"\n";
      //mf::LogInfo("EMSaturation")<<"fTotalEnergyDeposit: "<< step->GetTotalEnergyDeposit()/CLHEP::MeV<<"\nfNumIonElectrons: "<<fNumIonElectrons<<"\n";
    }
    else{
      fNumScintPhotons = fScintYieldFactor * scintYield * fEnergyDeposit;
      fVisibleEnergyDeposition = 0.0; //This is set to zero because I have not made a correct implimentation of this value for anything but EMSaturation.
    }

    MF_LOG_DEBUG("ISCalculationSeparate") << "number photons: " << fNumScintPhotons 
				       << " energy: "        << fEnergyDeposit/CLHEP::MeV
				       << " saturation: " 
				       << fEMSaturation->VisibleEnergyDepositionAtAStep(step)
				       << " step length: "   << step->GetStepLength()/CLHEP::cm;


    return;
  }


}// namespace
