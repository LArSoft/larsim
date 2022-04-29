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

#ifndef IS_ISCALCSEPARATE_H
#define IS_ISCALCSEPARATE_H

#include "larsim/IonizationScintillation/ISCalc.h"

namespace spacecharge {
  class SpaceCharge;
}

namespace detinfo {
  class LArProperties;
}

#include <utility>

namespace larg4 {
  class ISCalcSeparate : public ISCalc {
  public:
    ISCalcSeparate();
    void Reset();

    double EFieldAtStep(double efield,
                        sim::SimEnergyDeposit const& edep)
      override; //value of field with any corrections for this step
    ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                               sim::SimEnergyDeposit const& edep) override;

  private:
    double fGeVToElectrons; ///< from LArG4Parameters service
    double fRecombA;        ///< from LArG4Parameters service
    double fRecombk;        ///< from LArG4Parameters service
    double fModBoxA;        ///< from LArG4Parameters service
    double fModBoxB;        ///< from LArG4Parameters service
    bool fUseModBoxRecomb;  ///< from LArG4Parameters service

    const spacecharge::SpaceCharge* fSCE;
    const detinfo::LArProperties* fLArProp;

    double CalcIon(detinfo::DetectorPropertiesData const& detProp,
                   sim::SimEnergyDeposit const& edep);
    std::pair<double, double> CalcScint(sim::SimEnergyDeposit const& edep);
    double GetScintYieldRatio(sim::SimEnergyDeposit const& edep);
  };
}
#endif // LARG4_ISCALCULATIONSEPARATE_H
