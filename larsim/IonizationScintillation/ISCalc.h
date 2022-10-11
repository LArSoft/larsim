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

namespace detinfo {
  class DetectorPropertiesData;
  class LArProperties;
}

namespace sim {
  class SimEnergyDeposit;
}

namespace larg4 {
  struct ISCalcData {
    double energyDeposit;           // total energy deposited in the step
    double numElectrons;            // number of ionization electrons for this step
    double numPhotons;              // number of scintillation photons for this step
    double scintillationYieldRatio; // liquid argon scintillation yield ratio
  };

  class ISCalc {
  private:
    const detinfo::LArProperties* fLArProp;

  public:
    ISCalc();
    virtual ~ISCalc() = default;
    virtual ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                       sim::SimEnergyDeposit const& edep) = 0;
    virtual double EFieldAtStep(
      double efield,
      sim::SimEnergyDeposit const& edep) = 0; //value of field with any corrections for this step
    double GetScintYield(sim::SimEnergyDeposit const& edep, bool prescale);
    double GetScintYieldRatio(sim::SimEnergyDeposit const& edep);
  };
}
#endif // LARG4_ISCALC_H
