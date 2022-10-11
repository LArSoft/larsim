////////////////////////////////////////////////////////////////////////
/// \file  ISCalcuation.cc
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using separate algorithms for each
///
/// Wes, 18Feb2018: this is a copy of the original, for standalone purposes
///
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/ElectronDrift/ISCalculationSeparate.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include <cmath>
#include <ostream>

namespace {
  constexpr double scint_yield_factor{1.};
}

namespace detsim {

  //----------------------------------------------------------------------------
  ISCalculationSeparate::ISCalculationSeparate(fhicl::ParameterSet const& pset)
    : fRecombA{pset.get<double>("RecombA")}
    , fRecombk{pset.get<double>("Recombk")}
    , fModBoxA{pset.get<double>("ModBoxA")}
    , fModBoxB{pset.get<double>("ModBoxB")}
    , fUseModBoxRecomb{pset.get<bool>("UseModBoxRecomb")}
    , fGeVToElectrons{pset.get<double>("GeVToElectrons")}
  {
    fLArProp = lar::providerFrom<detinfo::LArPropertiesService>();
    fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    // the recombination coefficient is in g/(MeVcm^2), but
    // we report energy depositions in MeV/cm, need to divide
    // Recombk from the LArG4Parameters service by the density
    // of the argon we got above.  FIXME: is this being done?

    std::cout << "fLArG4Prop->GeVToElectrons(): " << fGeVToElectrons << std::endl;
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  double ISCalculationSeparate::CalculateIonization(detinfo::DetectorPropertiesData const& detProp,
                                                    sim::SimEnergyDeposit const& edep) const
  {
    float e = edep.Energy();
    float ds = edep.StepLength();
    float x = edep.MidPointX();
    float y = edep.MidPointY();
    float z = edep.MidPointZ();
    double recomb = 0.;
    double dEdx = (ds <= 0.0) ? 0.0 : e / ds;
    double EFieldStep = EFieldAtStep(detProp.Efield(), x, y, z);

    // Guard against spurious values of dE/dx. Note: assumes density of LAr
    if (dEdx < 1.) dEdx = 1.;

    if (fUseModBoxRecomb) {
      if (ds > 0) {
        double Xi = fModBoxB * dEdx / EFieldStep;
        recomb = log(fModBoxA + Xi) / Xi;
      }
      else
        recomb = 0;
    }
    else {
      recomb = fRecombA / (1. + dEdx * fRecombk / EFieldStep);
    }

    // 1.e-3 converts fEnergyDeposit to GeV
    double const numIonElectrons = fGeVToElectrons * 1.e-3 * e * recomb;

    MF_LOG_DEBUG("ISCalculationSeparate")
      << " Electrons produced for " << edep.Energy() << " MeV deposited with " << recomb
      << " recombination: " << numIonElectrons << std::endl;

    return numIonElectrons;
  }

  //----------------------------------------------------------------------------
  double ISCalculationSeparate::CalculateScintillation(sim::SimEnergyDeposit const& edep) const
  {
    float const e = edep.Energy();
    int const pdg = edep.PdgCode();
    double scintYield = fLArProp->ScintYield(true);
    if (fLArProp->ScintByParticleType()) {

      MF_LOG_DEBUG("ISCalculationSeparate") << "scintillating by particle type";

      switch (pdg) {

      case 2212: scintYield = fLArProp->ProtonScintYield(true); break;
      case 13:
      case -13: scintYield = fLArProp->MuonScintYield(true); break;
      case 211:
      case -211: scintYield = fLArProp->PionScintYield(true); break;
      case 321:
      case -321: scintYield = fLArProp->KaonScintYield(true); break;
      case 1000020040: scintYield = fLArProp->AlphaScintYield(true); break;
      case 11:
      case -11:
      case 22: scintYield = fLArProp->ElectronScintYield(true); break;
      default: scintYield = fLArProp->ElectronScintYield(true);
      }

      return scintYield * e;
    }

    return scint_yield_factor * scintYield * e;
  }
  //----------------------------------------------------------------------------
  ISCalculationSeparate::Data ISCalculationSeparate::CalculateIonizationAndScintillation(
    detinfo::DetectorPropertiesData const& detProp,
    sim::SimEnergyDeposit const& edep) const
  {
    return {edep.Energy(), CalculateIonization(detProp, edep), CalculateScintillation(edep)};
  }

  double ISCalculationSeparate::EFieldAtStep(double const efield,
                                             sim::SimEnergyDeposit const& edep) const
  {
    return EFieldAtStep(efield, edep.MidPointX(), edep.MidPointY(), edep.MidPointZ());
  }

  double ISCalculationSeparate::EFieldAtStep(double const efield,
                                             float const x,
                                             float const y,
                                             float const z) const
  {
    if (not fSCE->EnableSimEfieldSCE()) { return efield; }

    auto const offsets = fSCE->GetEfieldOffsets(geo::Point_t{x, y, z});
    return std::hypot(efield + efield * offsets.X(), efield * offsets.Y(), efield * offsets.Z());
  }

} // namespace
