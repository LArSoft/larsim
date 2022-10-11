////////////////////////////////////////////////////////////////////////
// \file  ISCalculationCorrelated.h
// \brief Interface to algorithm class for a specific calculation of
//        ionization electrons and scintillation photons, based on
//        simple microphysics arguments to establish an anticorrelation
//        between these two quantities.
//
//        To enable this in simulation, change LArG4Parameters variable
//        in your fhicl file:
//
//        services.LArG4Parameters.IonAndScintCalculator: "Correlated"
//
// \author wforeman @ iit.edu
// Modified: Adding corrections for low electric field (LArQL model)
// Mar 2021 by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATIONCORRELATED_H
#define LARG4_ISCALCULATIONCORRELATED_H

#include "larsim/LegacyLArG4/ISCalculation.h"
namespace detinfo {
  class DetectorPropertiesData;
}

// forward declarations
class G4Step;

namespace larg4 {

  class ISCalculationCorrelated : public ISCalculation {
  public:
    explicit ISCalculationCorrelated(detinfo::DetectorPropertiesData const& detProp);

    void Reset();
    void CalculateIonizationAndScintillation(const G4Step* step);
    double StepSizeLimit() const { return fStepSize; }

  private:
    double fStepSize;        ///< maximum step to take
    double fEfield;          ///< value of electric field from LArProperties service
    double fWion;            ///< W_ion (23.6 eV) == 1/fGeVToElectrons
    double fWph;             ///< W_ph (19.5 eV)
    double fScintPreScale;   ///< scintillation pre-scale from LArProperties service
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

    double EscapingEFraction(
      double const dEdx); //LArQL chi0 function = fraction of escaping electrons
    double FieldCorrection(
      double const EF,
      double const dEdx); //LArQL f_corr function = correction factor for electric field dependence
  };
}
#endif // LARG4_ISCALCULATIONCORRELATED_H
