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
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalcCorrelated.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace larg4 {
  //----------------------------------------------------------------------------
  ISCalcCorrelated::ISCalcCorrelated(detinfo::DetectorPropertiesData const& detProp)
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
    fGeVToElectrons = LArG4PropHandle->GeVToElectrons();

    // ionization work function
    fWion = 1. / fGeVToElectrons * 1e3; // MeV

    // ion+excitation work function (\todo: get from LArG4Parameters or LArProperties?)
    fWph = 19.5 * 1e-6; // MeV
  }

  //----------------------------------------------------------------------------
  ISCalcData
  ISCalcCorrelated::CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                    sim::SimEnergyDeposit const& edep)
  {
    double const energy_deposit = edep.Energy();

    // calculate total quanta (ions + excitons)
    double Nq = energy_deposit / fWph;

    float ds = edep.StepLength();
    double dEdx = (ds <= 0.0) ? 0.0 : energy_deposit / ds;
    double EFieldStep = EFieldAtStep(detProp.Efield(), edep);
    double recomb = 0.;

    // Guard against spurious values of dE/dx. Note: assumes density of LAr
    if (dEdx < 1.) dEdx = 1.;

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

    // using this recombination, calculate number of ionization electrons
    double const num_electrons = (energy_deposit / fWion) * recomb;

    // calculate scintillation photons
    double const num_photons = (Nq - num_electrons) * fScintPreScale;

    MF_LOG_DEBUG("ISCalcCorrelated")
      << " Electrons produced for " << energy_deposit << " MeV deposited with " << recomb
      << " recombination: " << num_electrons << std::endl;
    MF_LOG_DEBUG("ISCalcCorrelated") << "number photons: " << num_photons;

    return {energy_deposit, num_electrons, num_photons, GetScintYieldRatio(edep)};
  }

  //----------------------------------------------------------------------------
  double
  ISCalcCorrelated::GetScintYieldRatio(sim::SimEnergyDeposit const& edep)
  {
    // For ISCalcCorrelated, the ScintByParticleType option only controls
    // the scintillation yield ratio, which is the ratio of fast light (singlet
    // component) to the total light (singlet+triplet components).
    //
    // TODO: move this to ISCalc, since it is the same function used in the
    //       other ionization/scintillation calculation algs

    if (!fLArProp->ScintByParticleType()) return fLArProp->ScintYieldRatio();

    switch (edep.PdgCode()) {
    case 2212: return fLArProp->ProtonScintYieldRatio();
    case 13:
    case -13: return fLArProp->MuonScintYieldRatio();
    case 211:
    case -211: return fLArProp->PionScintYieldRatio();
    case 321:
    case -321: return fLArProp->KaonScintYieldRatio();
    case 1000020040: return fLArProp->AlphaScintYieldRatio();
    case 11:
    case -11:
    case 22: return fLArProp->ElectronScintYieldRatio();
    default: return fLArProp->ElectronScintYieldRatio();
    }
  }

  //----------------------------------------------------------------------------
  double
  ISCalcCorrelated::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
  {
    if (!fSCE->EnableSimEfieldSCE()) return efield;
    auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());
    return efield * std::hypot(1 + eFieldOffsets.X(), eFieldOffsets.Y(), eFieldOffsets.Z());
  }

} // namespace
