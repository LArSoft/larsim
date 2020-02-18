////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.h
/// \brief Interface to algorithm class for a specific calculation of
///        ionization electrons and scintillation photons assuming there
///        is no correlation between the two
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATIONSEPARATE_H
#define LARG4_ISCALCULATIONSEPARATE_H

#include "larsim/LegacyLArG4/ISCalculation.h"

// forward declarations
class G4EmSaturation;
class G4Step;
namespace CLHEP { class HepRandomEngine; }

namespace larg4 {

 class ISCalculationSeparate : public ISCalculation {

 public:

   ISCalculationSeparate(CLHEP::HepRandomEngine&);
   virtual ~ISCalculationSeparate();

   void   Initialize();
   void   Reset();
   void   CalculateIonizationAndScintillation(const G4Step* step);
   double StepSizeLimit()              const { return fStepSize;            }

 private:

   double                fStepSize;            ///< maximum step to take
   double                fEfield;              ///< value of electric field from LArProperties service
   double 	   	 fGeVToElectrons;      ///< conversion factor from LArProperties service
   double 	   	 fRecombA;             ///< from LArG4Parameters service
   double 	   	 fRecombk;             ///< from LArG4Parameters service
   double 	   	 fModBoxA;             ///< from LArG4Parameters service
   double 	   	 fModBoxB;             ///< from LArG4Parameters service
   bool   	   	 fUseModBoxRecomb;     ///< from LArG4Parameters service
   bool   	   	 fScintByParticleType; ///< from LArProperties service
   double 	   	 fScintYieldFactor;    ///< scintillation yield factor
   G4EmSaturation* 	 fEMSaturation;        ///< pointer to EM saturation
 };
}
#endif // LARG4_ISCALCULATIONSEPARATE_H
