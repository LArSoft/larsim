////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationCorrelated.h
/// \brief Interface to algorithm class for a specific calculation of
///        ionization electrons and scintillation photons, based on 
///        simple microphysics arguments to establish an anticorrelation
///        between these two quantities.
///
///        To enable this in simulation, change LArG4Parameters variable
///        in your fhicl file:
///       
///        services.LArG4Parameters.IonAndScintCalculator: "Correlated" 
///
/// \author wforeman @ iit.edu
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATIONCORRELATED_H
#define LARG4_ISCALCULATIONCORRELATED_H

#include "larsim/LArG4/ISCalculation.h"

// forward declarations
class G4Step;
namespace CLHEP { class HepRandomEngine; }

namespace larg4 {

 class ISCalculationCorrelated : public ISCalculation {

 public:

   ISCalculationCorrelated(CLHEP::HepRandomEngine&);
   virtual ~ISCalculationCorrelated();

   void   Initialize();
   void   Reset();
   void   CalculateIonizationAndScintillation(const G4Step* step);
   double StepSizeLimit()              const { return fStepSize;            }

 private:

   double                fStepSize;            ///< maximum step to take
   double                fEfield;              ///< value of electric field from LArProperties service
   double                fWion;                ///< W_ion (23.6 eV) == 1/fGeVToElectrons
   double                fWph;                 ///< W_ph (19.5 eV) 
   double                fRecombA;             ///< from LArG4Parameters service
   double                fRecombk;             ///< from LArG4Parameters service
   double                fModBoxA;             ///< from LArG4Parameters service
   double                fModBoxB;             ///< from LArG4Parameters service
   bool                  fUseModBoxRecomb;     ///< from LArG4Parameters service
 };
}
#endif // LARG4_ISCALCULATIONCORRELATED_H
