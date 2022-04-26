////////////////////////////////////////////////////////////////////////
// Class:       ISCalcCorrelated
// Plugin Type: algorithm
// File:        ISCalcCorrelated.h and ISCalcCorrelated.cxx
// Description: Interface to algorithm class for a specific calculation of
//              ionization electrons and scintillation photons, based on
//              simple microphysics arguments to establish an anticorrelation
//              between these two quantities.
// Input: 'sim::SimEnergyDeposit'
// Output: num of Photons and Electrons
// May 2020 by W Foreman
// Modified: Adding corrections for low electric field (LArQL model)
// Jun 2020 by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalcCorrelated.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "CLHEP/Random/RandGauss.h" 
#include "CLHEP/Random/RandFlat.h"

namespace larg4 {
  //----------------------------------------------------------------------------
  ISCalcCorrelated::ISCalcCorrelated(detinfo::DetectorPropertiesData const& detProp, CLHEP::HepRandomEngine& Engine)
    : fISTPC{*(lar::providerFrom<geo::Geometry>())}
    , fEngine(Engine)
  {
    std::cout << "IonizationAndScintillation/ISCalcCorrelated Initialize." << std::endl;
    art::ServiceHandle<sim::LArG4Parameters const> LArG4PropHandle;

    fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    fLArProp = lar::providerFrom<detinfo::LArPropertiesService>();
    fScintPreScale = fLArProp->ScintPreScale();

    //the recombination coefficient is in g/(MeVcm^2), but we report energy depositions in MeV/cm,
    //need to divide Recombk from the LArG4Parameters service by the density of the argon we got above.
    fRecombA = LArG4PropHandle->RecombA();
    fRecombk = LArG4PropHandle->Recombk() / detProp.Density(detProp.Temperature());
    fModBoxA = LArG4PropHandle->ModBoxA();
    fModBoxB = LArG4PropHandle->ModBoxB() / detProp.Density(detProp.Temperature());
    fUseModBoxRecomb = (bool)LArG4PropHandle->UseModBoxRecomb();
    fUseModLarqlRecomb  = (bool)LArG4PropHandle->UseModLarqlRecomb();
    fLarqlChi0A = LArG4PropHandle->LarqlChi0A();
    fLarqlChi0B = LArG4PropHandle->LarqlChi0B();
    fLarqlChi0C = LArG4PropHandle->LarqlChi0C();
    fLarqlChi0D = LArG4PropHandle->LarqlChi0D();
    fLarqlAlpha = LArG4PropHandle->LarqlAlpha();
    fLarqlBeta  = LArG4PropHandle->LarqlBeta();
    fGeVToElectrons = LArG4PropHandle->GeVToElectrons();

    // ionization work function
    fWion = 1. / fGeVToElectrons * 1e3; // MeV

    // ion+excitation work function
    fWph = LArG4PropHandle->Wph() * 1e-6; // MeV
  }

  //----------------------------------------------------------------------------
  ISCalcData
  ISCalcCorrelated::CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                    sim::SimEnergyDeposit const& edep)
  {

    CLHEP::RandGauss GaussGen(fEngine); 
    CLHEP::RandFlat UniformGen(fEngine);

    double const energy_deposit = edep.Energy();

    // calculate total quanta (ions + excitons)
    double Nq = energy_deposit / fWph;

    double num_ions = energy_deposit / fWion;
    double num_quanta = energy_deposit / fWph;

    float ds = edep.StepLength();
    double dEdx = (ds <= 0.0) ? 0.0 : energy_deposit / ds;
    double EFieldStep = EFieldAtStep(detProp.Efield(), edep);
    double recomb = 0.;

    //calculate recombination survival fraction value inside, otherwise zero

    if(EFieldStep > 0) {
      // Guard against spurious values of dE/dx. Note: assumes density of LAr
      //if (dEdx < 1.) dEdx = 1.;

      // calculate recombination survival fraction
      if (fUseModBoxRecomb) {
	if (ds > 0) {
	  double Xi = fModBoxB * dEdx / EFieldStep;
	  recomb = log(fModBoxA + Xi) / Xi;
	}
	else {
	  recomb = 0;
	}
      }
      else {
	recomb = fRecombA / (1. + dEdx * fRecombk / EFieldStep);
      }

    }//Efield

    if(fUseModLarqlRecomb){ //Use corrections from LArQL model
      recomb += EscapingEFraction(dEdx)*FieldCorrection(EFieldStep, dEdx); //Correction for low EF
    }

    // using this recombination, calculate number of ionization electrons
    int const num_electrons = BinomFluct(num_ions, recomb);
    // calculate scintillation photons
    int const num_photons = floor(num_quanta - num_electrons);

    MF_LOG_DEBUG("ISCalcCorrelated")
      << " Electrons produced for " << energy_deposit << " MeV deposited with " << recomb
      << " recombination: " << num_electrons << std::endl;
    MF_LOG_DEBUG("ISCalcCorrelated") << "number photons: " << num_photons;

    return {energy_deposit, num_electrons, num_photons, GetScintYieldRatio(edep)};
  }

  //----------------------------------------------------------------------------
   double
   ISCalcCorrelated::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
   {
     //electric field outside active volume set to zero                                                                                                                                                      
     if(!fISTPC.isScintInActiveVolume(edep.MidPoint())) return 0;

     //electric field inside active volume                                                                                                                                                                   
     if (!fSCE->EnableSimEfieldSCE()) return efield;

     auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());
     return efield * std::hypot(1 + eFieldOffsets.X(), eFieldOffsets.Y(), eFieldOffsets.Z());
   }



  double ISCalcCorrelated::EscapingEFraction(double const dEdx){ //LArQL chi0 function = fraction of escaping electrons
    return fLarqlChi0A/(fLarqlChi0B+exp(fLarqlChi0C+fLarqlChi0D*dEdx));
  }

  double ISCalcCorrelated::FieldCorrection(double const EF, double const dEdx){ //LArQL f_corr function = correction factor for electric field dependence
    return exp(-EF/(fLarqlAlpha*log(dEdx)+fLarqlBeta));
  }

  int
  ISCalcCorrelated::BinomFluct(int N0, double prob)
  {
    CLHEP::RandGauss GaussGen(fEngine);
    CLHEP::RandFlat UniformGen(fEngine);

    double mean = N0 * prob;
    double sigma = sqrt(N0 * prob * (1 - prob));
    int N1 = 0;
    if (prob == 0.00) { return N1; }
    if (prob == 1.00) { return N0; }

    if (N0 < 10) {
      for (int i = 0; i < N0; i++) {
        if (UniformGen.fire() < prob) { N1++; }
      }
    }
    else {
      N1 = int(floor(GaussGen.fire(mean, sigma) + 0.5));
    }
    if (N1 > N0) { N1 = N0; }
    if (N1 < 0) { N1 = 0; }
    return N1;
  }
  
} // namespace
