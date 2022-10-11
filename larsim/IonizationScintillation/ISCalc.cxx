////////////////////////////////////////////////////////////////////////
//Class: ISCalc
//File: ISCalc.h and ISCalc.cxx
//Description: Interface to the ion and scintillation calculation algorithm
//2022
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalc.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace larg4 {

  //----------------------------------------------------------------------------
  ISCalc::ISCalc() : fLArProp{lar::providerFrom<detinfo::LArPropertiesService>()} {}

  //----------------------------------------------------------------------------
  double ISCalc::GetScintYield(sim::SimEnergyDeposit const& edep, bool prescale)
  {
    if (!fLArProp->ScintByParticleType()) return fLArProp->ScintYield(true);
    switch (edep.PdgCode()) {
    case 2212: return fLArProp->ProtonScintYield(true);
    case 13:
    case -13: return fLArProp->MuonScintYield(true);
    case 211:
    case -211: return fLArProp->PionScintYield(true);
    case 321:
    case -321: return fLArProp->KaonScintYield(true);
    case 1000020040: return fLArProp->AlphaScintYield(true);
    case 11:
    case -11:
    case 22: return fLArProp->ElectronScintYield(true);
    default: return fLArProp->ElectronScintYield(true);
    }
  }

  //----------------------------------------------------------------------------
  double ISCalc::GetScintYieldRatio(sim::SimEnergyDeposit const& edep)
  {
    // ScintByParticleType option only controls the scintillation
    // yield ratio, which is the ratio of fast light (singlet
    // component) to the total light (singlet+triplet components).

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

}
