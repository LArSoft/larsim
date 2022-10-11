////////////////////////////////////////////////////////////////////////
// Class:       ISCalcCorrelated
// Plugin Type: algorithm
// File:        ISCalcCorrelated.h and ISCalcCorrelated.cxx
// Description: Interface to algorithm class for a specific calculation of
//              ionization electrons and scintillation photons, based on
//              simple microphysics arguments to establish an anticorrelation
//              between these two quantities.
// Input: 'sim::SimEnergyDeposit'
// Output: Number of Photons and Electrons
// May 2020 by W Foreman
// Modified: Adding corrections for low electric field (LArQL model)
// Jun 2020 by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#ifndef IS_ISCALCCORRELATED_H
#define IS_ISCALCCORRELATED_H

#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/IonizationScintillation/ISTPC.h"

#include "CLHEP/Random/RandBinomial.h"

namespace spacecharge {
  class SpaceCharge;
}

namespace CLHEP {
  class HepRandomEngine;
}

namespace larg4 {
  class ISCalcCorrelated : public ISCalc {
  public:
    ISCalcCorrelated(detinfo::DetectorPropertiesData const& detProp,
                     CLHEP::HepRandomEngine& Engine);

    double EFieldAtStep(double efield,
                        sim::SimEnergyDeposit const& edep)
      override; //value of field with any corrections for this step
    ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                               sim::SimEnergyDeposit const& edep) override;

  private:
    ISTPC fISTPC;
    const spacecharge::SpaceCharge* fSCE;
    CLHEP::RandBinomial fBinomialGen;

    double fGeVToElectrons;  ///< from LArG4Parameters service
    double fWion;            ///< W_ion (23.6 eV) == 1/fGeVToElectrons
    double fWph;             ///< from LArG4Parameters service
    double fScintPreScale;   ///< scintillation pre-scaling factor from LArProperties service
    double fRecombA;         ///< from LArG4Parameters service
    double fRecombk;         ///< from LArG4Parameters service
    double fModBoxA;         ///< from LArG4Parameters service
    double fModBoxB;         ///< from LArG4Parameters service
    double fLarqlChi0A;      ///< from LArG4Parameters service
    double fLarqlChi0B;      ///< from LArG4Parameters service
    double fLarqlChi0C;      ///< from LArG4Parameters service
    double fLarqlChi0D;      ///< from LArG4Parameters service
    double fLarqlAlpha;      ///< from LArG4Parameters service
    double fLarqlBeta;       ///< from LArG4Parameters service
    bool fUseModBoxRecomb;   ///< from LArG4Parameters service
    bool fUseModLarqlRecomb; ///< from LArG4Parameters service
    bool fUseBinomialFlucts; ///< from LArG4Parameters service

    double EscapingEFraction(
      double const dEdx) const; //LArQL chi0 function = fraction of escaping electrons
    double FieldCorrection(double const EF, double const dEdx)
      const; //LArQL f_corr function = correction factor for electric field dependence
  };
}
#endif
