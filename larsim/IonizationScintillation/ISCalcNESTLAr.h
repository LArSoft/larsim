////////////////////////////////////////////////////////////////////////
// Class:       ISCalcNESTLAr
// Plugin Type: Algorithm
// File:        ISCalcNESTLAr.cxx
// Description:
// Aug. 30 by Mu Wei
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_ISCALCNESTLAr_H
#define LARG4_ISCALCNESTLAr_H

#include "larsim/IonizationScintillation/ISCalc.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
namespace detinfo {
  class DetectorPropertiesData;
}

#include "CLHEP/Random/RandEngine.h"

namespace larg4 {
  class ISCalcNESTLAr : public ISCalc {
  public:
    explicit ISCalcNESTLAr(CLHEP::HepRandomEngine& fEngine);
    void Reset();

    double EFieldAtStep(double efield,
                        sim::SimEnergyDeposit const& edep)
      override; //value of field with any corrections for this step
    ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                               sim::SimEnergyDeposit const& edep) override;

  private:
    CLHEP::HepRandomEngine& fEngine; // random engine
    const spacecharge::SpaceCharge* fSCE;
    const detinfo::LArProperties* fLArProp;

    int BinomFluct(int N0, double prob);
    double CalcElectronLET(double E);
    double GetScintYieldRatio(sim::SimEnergyDeposit const& edep);
  };
}
#endif // LARG4_ISCALCNESTLAr_H
