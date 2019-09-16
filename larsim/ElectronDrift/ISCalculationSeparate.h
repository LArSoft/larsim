////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.h
/// \brief Calculation of ionization electrons and scintillation photons 
///        assuming there is no correlation between the two
///
///
/// \version $Id:  $
/// \author  wenzel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef IS_ISCALCULATION_H
#define IS_ISCALCULATION_H

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// forward declaration
namespace detinfo { class LArProperties; class DetectorProperties; }
namespace fhicl { class ParameterSet; }
namespace sim { class SimEnergyDeposit; class LArG4Parameters; }
namespace spacecharge { class SpaceCharge; }

namespace detsim {

 class ISCalculationSeparate{

 public:

   explicit ISCalculationSeparate(fhicl::ParameterSet const& pset);
   virtual ~ISCalculationSeparate();

   //   void   Initialize(const detinfo::LArProperties* larp,
   // 		     const detinfo::DetectorProperties* detp,
   //		     const sim::LArG4Parameters* lgp,
   //		     const spacecharge::SpaceCharge* sce);
   //void   Reset();
   void   CalculateIonizationAndScintillation(sim::SimEnergyDeposit const& edep);
   double EnergyDeposit()              const { return fEnergyDeposit;       }
   double NumberIonizationElectrons()  const { return fNumIonElectrons;     }
   double NumberScintillationPhotons() const { return fNumScintPhotons;     }

   double EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep); //value of field with any corrections for this step

 private:

   double 	   	 fRecombA;         
   double 	   	 fRecombk;         
   double 	   	 fModBoxA;         
   double 	   	 fModBoxB;         
   bool   	   	 fUseModBoxRecomb; 
   double                fGeVToElectrons;      ///< from LArG4Parameters service
   double 	   	 fScintYieldFactor;    ///< scintillation yield factor

   double fEnergyDeposit;   ///< total energy deposited in the step
   double fNumIonElectrons; ///< number of ionization electrons for this step
   double fNumScintPhotons; ///< number of scintillation photons for this step

   geo::Vector_t fEfieldOffsets;
   double EFieldAtStep(double efield, float x, float y, float z);

   const detinfo::LArProperties*      fLArProp;
   const spacecharge::SpaceCharge*    fSCE;
   const detinfo::DetectorProperties* fDetProp;
   const sim::LArG4Parameters*        fLArG4Prop;

   //   void CalculateIonization(float e, float ds, float x, float y, float z);
   void CalculateIonization(sim::SimEnergyDeposit const& edep);

   void CalculateScintillation(float e, int pdg);
   void CalculateScintillation(sim::SimEnergyDeposit const& edep);

 };
}
#endif // LARG4_ISCALCULATION_H
