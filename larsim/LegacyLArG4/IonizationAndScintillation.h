////////////////////////////////////////////////////////////////////////
/// \file  IonizationAndScintillation.h
/// \brief Singleton to access a unified treatment of ionization and
///        scintillation in LAr
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_IONIZATIONANDSCINTILLATION_H
#define LARG4_IONIZATIONANDSCINTILLATION_H

#include <memory>
#include <string>

#include "larsim/LegacyLArG4/ISCalculation.h"

class G4Step;
class TH1F;
class TH2F;

namespace CLHEP {
  class HepRandomEngine;
}

namespace detinfo {
  class DetectorPropertiesData;
}

namespace larg4 {

  // The Ionization and Scintillation singleton
  class IonizationAndScintillation {
  public:
    static IonizationAndScintillation* CreateInstance(
      detinfo::DetectorPropertiesData const& detProp,
      CLHEP::HepRandomEngine& engine);
    static IonizationAndScintillation* Instance();

    // Method to reset the internal variables held in the ISCalculation
    // This method should be called at the start of any G4Step
    void Reset(const G4Step* step);

    double EnergyDeposit() const { return fISCalc->EnergyDeposit(); }
    double VisibleEnergyDeposit() const { return fISCalc->VisibleEnergyDeposit(); }
    double NumberIonizationElectrons() const { return fISCalc->NumberIonizationElectrons(); }
    double NumberScintillationPhotons() const { return fISCalc->NumberScintillationPhotons(); }
    double StepSizeLimit() const { return fISCalc->StepSizeLimit(); }

  private:
    IonizationAndScintillation(detinfo::DetectorPropertiesData const& detProp,
                               CLHEP::HepRandomEngine& engine);

    std::unique_ptr<larg4::ISCalculation>
      fISCalc;                    ///< object to calculate ionization and scintillation
                                  ///< produced by an energy deposition
    std::string fISCalculator;    ///< name of calculator to use, NEST or Separate
    G4Step const* fStep{nullptr}; ///< pointer to the current G4 step
    int fStepNumber{-1};          ///< last StepNumber checked
    int fTrkID{-1};               ///< last TrkID checked

    TH1F* fElectronsPerStep{nullptr};   ///< histogram of electrons per step
    TH1F* fStepSize{nullptr};           ///< histogram of the step sizes
    TH1F* fPhotonsPerStep{nullptr};     ///< histogram of the photons per step
    TH1F* fEnergyPerStep{nullptr};      ///< histogram of the energy deposited per step
    TH1F* fElectronsPerLength;          ///< histogram of electrons per cm
    TH1F* fPhotonsPerLength;            ///< histogram of photons per cm
    TH1F* fElectronsPerEDep;            ///< histogram of electrons per MeV deposited
    TH1F* fPhotonsPerEDep;              ///< histogram of photons per MeV deposited
    TH2F* fElectronsVsPhotons{nullptr}; ///< histogram of electrons vs photons per step
    CLHEP::HepRandomEngine& fEngine;    ///< random engine (needed for NEST)
  };

} // namespace larg4

#endif // LARG4_IONIZATIONANDSCINTILLATION
