////////////////////////////////////////////////////////////////////////
// Class:       ISCalcSeparate
// Plugin Type: algorithm
// File:        ISCalcSeparate.h and ISCalcSeparate.cxx
// Description:
// Interface to algorithm class for a specific calculation of ionization electrons and scintillation photons
// assuming there is no correlation between the two
// Input: 'sim::SimEnergyDeposit'
// Output: num of Photons and Electrons
// Sept.16 by Mu Wei
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalcSeparate.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <numeric>

namespace larg4 {
  //----------------------------------------------------------------------------
  ISCalcSeparate::ISCalcSeparate()
  {
    fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    // the recombination coefficient is in g/(MeVcm^2), but we report
    // energy depositions in MeV/cm, need to divide Recombk from the
    // LArG4Parameters service by the density of the argon we got
    // above; this is done in 'CalcIon' function below.
    art::ServiceHandle<sim::LArG4Parameters const> LArG4PropHandle;
    fRecombA = LArG4PropHandle->RecombA();
    fRecombk = LArG4PropHandle->Recombk();
    fModBoxA = LArG4PropHandle->ModBoxA();
    fModBoxB = LArG4PropHandle->ModBoxB();
    fUseModBoxRecomb = (bool)LArG4PropHandle->UseModBoxRecomb();
    fGeVToElectrons = LArG4PropHandle->GeVToElectrons();
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  double ISCalcSeparate::CalcIon(detinfo::DetectorPropertiesData const& detProp,
                                 sim::SimEnergyDeposit const& edep)
  {
    float e = edep.Energy();
    float ds = edep.StepLength();

    double recomb = 0.;
    double dEdx = (ds <= 0.0) ? 0.0 : e / ds;
    double EFieldStep = EFieldAtStep(detProp.Efield(), edep);

    // Guard against spurious values of dE/dx. Note: assumes density of LAr
    if (dEdx < 1.) { dEdx = 1.; }

    if (fUseModBoxRecomb) {
      if (ds > 0) {
        double const scaled_modboxb = fModBoxB / detProp.Density(detProp.Temperature());
        double const Xi = scaled_modboxb * dEdx / EFieldStep;
        recomb = log(fModBoxA + Xi) / Xi;
      }
      else {
        recomb = 0;
      }
    }
    else {
      double const scaled_recombk = fRecombk / detProp.Density(detProp.Temperature());
      recomb = fRecombA / (1. + dEdx * scaled_recombk / EFieldStep);
    }

    // 1.e-3 converts fEnergyDeposit to GeV
    auto const numIonElectrons = fGeVToElectrons * 1.e-3 * e * recomb;

    MF_LOG_DEBUG("ISCalcSeparate")
      << " Electrons produced for " << edep.Energy() << " MeV deposited with " << recomb
      << " recombination: " << numIonElectrons << std::endl;
    return numIonElectrons;
  }

  //----------------------------------------------------------------------------
  std::pair<double, double> ISCalcSeparate::CalcScint(sim::SimEnergyDeposit const& edep)
  {
    double numScintPhotons = GetScintYield(edep, true) * edep.Energy();
    return {numScintPhotons, GetScintYieldRatio(edep)};
  }

  //----------------------------------------------------------------------------
  ISCalcData ISCalcSeparate::CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                             sim::SimEnergyDeposit const& edep)
  {
    auto const numElectrons = CalcIon(detProp, edep);
    auto const [numPhotons, scintYieldRatio] = CalcScint(edep);
    return {edep.Energy(), numElectrons, numPhotons, scintYieldRatio};
  }
  //----------------------------------------------------------------------------
  double ISCalcSeparate::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
  {
    if (not fSCE->EnableSimEfieldSCE()) { return efield; }

    auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());
    return std::hypot(
      efield + efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
  }

} // namespace
