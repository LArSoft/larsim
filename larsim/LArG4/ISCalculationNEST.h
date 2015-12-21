////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationNEST.h
/// \brief Interface to algorithm class for a specific calculation of 
///        ionization electrons and scintillation photons using NEST
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATIONNEST_H
#define LARG4_ISCALCULATIONNEST_H

#include "larsim/LArG4/ISCalculation.h"
#include "larsim/LArG4/NestAlg.h"

namespace larg4 {

  class ISCalculationNEST : public ISCalculation {

 public:

   ISCalculationNEST();
   virtual ~ISCalculationNEST();

   void   Initialize();
   void   Reset();
   void   CalculateIonizationAndScintillation(const G4Step* step);
   double EnergyDeposit()              const { return fEnergyDeposit;   }
   int    NumberIonizationElectrons()  const { return fNumIonElectrons; }
   int    NumberScintillationPhotons() const { return fNumScintPhotons; }
   double StepSizeLimit()              const { return fStepSize;        }
   
 private:

   NestAlg* fNest;     ///< the fast optical simulation process
   double   fStepSize; ///< maximum step to take
 };
}
#endif // LARG4_ISCALCULATIONNEST_H

