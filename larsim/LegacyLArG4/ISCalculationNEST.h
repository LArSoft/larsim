////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationNEST.h
/// \brief Interface to algorithm class for a specific calculation of
///        ionization electrons and scintillation photons using NEST
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATIONNEST_H
#define LARG4_ISCALCULATIONNEST_H

#include "larsim/LegacyLArG4/ISCalculation.h"
#include "larsim/LegacyLArG4/NestAlg.h"

// forward declarations
class G4Step;
namespace CLHEP {
  class HepRandomEngine;
}

namespace larg4 {

  class ISCalculationNEST : public ISCalculation {
  public:
    ISCalculationNEST(CLHEP::HepRandomEngine& engine);

    void Reset();
    void CalculateIonizationAndScintillation(const G4Step* step);
    double StepSizeLimit() const { return fStepSize; }

  private:
    NestAlg fNest;    ///< the fast optical simulation process
    double fStepSize; ///< maximum step to take
  };
}
#endif // LARG4_ISCALCULATIONNEST_H
