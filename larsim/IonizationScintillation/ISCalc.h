////////////////////////////////////////////////////////////////////////
// Class:       ISCalc
// Plugin Type: Algorithm
// File:        ISCalc.h and ISCalc.cxx
// Description: Interface to the ion and scintillation calculation algorithm
// Input: 'sim::SimEnergyDeposit'
// Output: number of photons and electrons
// Spet.16 by Mu Wei
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_ISCALC_H
#define LARG4_ISCALC_H

#include "lardataobj/Simulation/SimEnergyDeposit.h"
namespace detinfo {
  class DetectorPropertiesData;
}

namespace larg4 {
  struct ISCalcData {
    double energyDeposit;           // total energy deposited in the step
    double numElectrons;            // number of ionization electrons for this step
    double numPhotons;              // number of scintillation photons for this step
    double scintillationYieldRatio; // liquid argon scintillation yield ratio
  };

  class ISCalc {
  public:
    virtual ~ISCalc() = default;
    virtual ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                       sim::SimEnergyDeposit const& edep) = 0;
    virtual double EFieldAtStep(
      double efield,
      sim::SimEnergyDeposit const& edep) = 0; //value of field with any corrections for this step
  };
}
#endif // LARG4_ISCALC_H
