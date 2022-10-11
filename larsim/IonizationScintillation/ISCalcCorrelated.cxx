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
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "CLHEP/Random/RandBinomial.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace larg4 {
  //----------------------------------------------------------------------------
  ISCalcCorrelated::ISCalcCorrelated(detinfo::DetectorPropertiesData const& detProp,
                                     CLHEP::HepRandomEngine& Engine)
    : fISTPC{*(lar::providerFrom<geo::Geometry>())}
    , fSCE(lar::providerFrom<spacecharge::SpaceChargeService>())
    , fBinomialGen{CLHEP::RandBinomial(Engine)}
  {
    MF_LOG_INFO("ISCalcCorrelated") << "IonizationAndScintillation/ISCalcCorrelated Initialize.";

    fScintPreScale = lar::providerFrom<detinfo::LArPropertiesService>()->ScintPreScale();

    art::ServiceHandle<sim::LArG4Parameters const> LArG4PropHandle;

    // The recombination coefficient is in g/(MeVcm^2), but we report
    // energy depositions in MeV/cm, need to divide Recombk from the
    // LArG4Parameters service by the density of the argon we got
    // above.
    fRecombA = LArG4PropHandle->RecombA();
    fRecombk = LArG4PropHandle->Recombk() / detProp.Density(detProp.Temperature());
    fModBoxA = LArG4PropHandle->ModBoxA();
    fModBoxB = LArG4PropHandle->ModBoxB() / detProp.Density(detProp.Temperature());
    fUseModBoxRecomb = (bool)LArG4PropHandle->UseModBoxRecomb();
    fUseModLarqlRecomb = (bool)LArG4PropHandle->UseModLarqlRecomb();
    fUseBinomialFlucts = (bool)LArG4PropHandle->UseBinomialFlucts();
    fLarqlChi0A = LArG4PropHandle->LarqlChi0A();
    fLarqlChi0B = LArG4PropHandle->LarqlChi0B();
    fLarqlChi0C = LArG4PropHandle->LarqlChi0C();
    fLarqlChi0D = LArG4PropHandle->LarqlChi0D();
    fLarqlAlpha = LArG4PropHandle->LarqlAlpha();
    fLarqlBeta = LArG4PropHandle->LarqlBeta();
    fGeVToElectrons = LArG4PropHandle->GeVToElectrons();

    // ionization work function
    fWion = 1. / fGeVToElectrons * 1e3; // MeV

    // ion+excitation work function
    fWph = LArG4PropHandle->Wph() * 1e-6; // MeV
  }

  //----------------------------------------------------------------------------
  ISCalcData ISCalcCorrelated::CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                               sim::SimEnergyDeposit const& edep)
  {

    double const energy_deposit = edep.Energy();

    // calculate total quanta (ions + excitons)
    double num_ions = energy_deposit / fWion;
    double num_quanta = energy_deposit / fWph;

    double ds = edep.StepLength();
    double dEdx = (ds <= 0.0) ? 0.0 : energy_deposit / ds;
    dEdx = (dEdx < 1.) ? 1. : dEdx;
    double EFieldStep = EFieldAtStep(detProp.Efield(), edep);
    double recomb = 0.;

    //calculate recombination survival fraction value inside, otherwise zero
    if (EFieldStep > 0.) {
      // calculate recombination survival fraction
      // ...using Modified Box model
      if (fUseModBoxRecomb) {
        double Xi = fModBoxB * dEdx / EFieldStep;
        recomb = std::log(fModBoxA + Xi) / Xi;
      }
      // ... or using Birks/Doke
      else {
        recomb = fRecombA / (1. + dEdx * fRecombk / EFieldStep);
      }
    }

    if (fUseModLarqlRecomb) { //Use corrections from LArQL model
      recomb += EscapingEFraction(dEdx) * FieldCorrection(EFieldStep, dEdx); //Correction for low EF
    }

    // Guard against unphysical recombination values
    if (recomb < 0.) {
      mf::LogWarning("ISCalcCorrelated")
        << "Recombination survival fraction is lower than 0.: " << recomb << ", fixing it to 0.";
      recomb = 0.;
    }
    else if (recomb > 1.) {
      mf::LogWarning("ISCalcCorrelated")
        << "Recombination survival fraction is higher than 1.: " << recomb << ", fixing it to 1.";
      recomb = 1.;
    }

    // using this recombination, calculate number of ionization electrons
    double num_electrons =
      (fUseBinomialFlucts) ? fBinomialGen.fire(num_ions, recomb) : (num_ions * recomb);

    // calculate scintillation photons
    double num_photons = (num_quanta - num_electrons) * fScintPreScale;

    MF_LOG_DEBUG("ISCalcCorrelated")
      << "With " << energy_deposit << " MeV of deposited energy, "
      << "and a recombination of " << recomb << ", \nthere are " << num_electrons
      << " electrons, and " << num_photons << " photons.";

    return {energy_deposit, num_electrons, num_photons, GetScintYieldRatio(edep)};
  }

  //----------------------------------------------------------------------------
  double ISCalcCorrelated::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
  {
    // electric field outside active volume set to zero
    if (!fISTPC.isScintInActiveVolume(edep.MidPoint())) return 0.;

    // electric field inside active volume
    if (!fSCE->EnableSimEfieldSCE()) return efield;

    auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());
    return efield * std::hypot(1 + eFieldOffsets.X(), eFieldOffsets.Y(), eFieldOffsets.Z());
  }

  //----------------------------------------------------------------------------
  // LArQL chi0 function = fraction of escaping electrons
  double ISCalcCorrelated::EscapingEFraction(double const dEdx) const
  {
    return fLarqlChi0A / (fLarqlChi0B + std::exp(fLarqlChi0C + fLarqlChi0D * dEdx));
  }

  //----------------------------------------------------------------------------
  // LArQL f_corr function = correction factor for electric field dependence
  double ISCalcCorrelated::FieldCorrection(double const EF, double const dEdx) const
  {
    return std::exp(-EF / (fLarqlAlpha * std::log(dEdx) + fLarqlBeta));
  }

} // namespace larg4
