////////////////////////////////////////////////////////////////////////
/// \file  ISCalculation.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATION_H
#define LARG4_ISCALCULATION_H

#include "Geant4/G4Step.hh"
#include "Geant4/G4UserLimits.hh"

#include "larsim/LArG4/OpFastScintillation.hh"

namespace larg4{

 class ISCalculation{

 public:

   ISCalculation();
   virtual ~ISCalculation();

   virtual void                 Initialize()                                            = 0;
   virtual void                 Reset()                            		      = 0;
   virtual void                 CalculateIonizationAndScintillation(const G4Step* step) = 0;
   virtual double               EnergyDeposit()              const                      = 0;
   virtual int    	      NumberIonizationElectrons()  const 		      = 0;
   virtual int    	      NumberScintillationPhotons() const 		      = 0;
   virtual double               StepSizeLimit()              const 		      = 0;

   //Method to get electric field
   double EFieldAtStep(double fEfield, const G4Step* step) const; //value of field with any corrections for this step  

 protected:

   double fEnergyDeposit;   ///< total energy deposited in the step
   int    fNumIonElectrons; ///< number of ionization electrons for this step
   int    fNumScintPhotons; ///< number of scintillation photons for this step   
   
 };
}
#endif // LARG4_ISCALCULATION_H

