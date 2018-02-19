////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.h
/// \brief Interface to algorithm class for a specific calculation of 
///        ionization electrons and scintillation photons assuming there
///        is no correlation between the two
///
/// Wes, 18Feb2018: this is a copy of the original, for standalone purposes
///
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef IS_ISCALCULATIONSEPARATE_H
#define IS_ISCALCULATIONSEPARATE_H

#include <vector>

// forward declaration
namespace detinfo { class LArProperties; class DetectorProperties; }
namespace sim { class SimEnergyDeposit; class LArG4Parameters; }
namespace spacecharge { class SpaceCharge; }

namespace larg4 {

 class ISCalculationSeparate{

 public:

   ISCalculationSeparate();
   virtual ~ISCalculationSeparate();

   void   Initialize(const detinfo::LArProperties* larp,   
		     const detinfo::DetectorProperties* detp,
		     const sim::LArG4Parameters* lgp,
		     const spacecharge::SpaceCharge* sce);
   void   Reset();
   void   CalculateIonizationAndScintillation(sim::SimEnergyDeposit const& edep);
   void   CalculateIonizationAndScintillation(sim::SimEnergyDeposit const& edep, std::vector<double> const&);
   double EnergyDeposit()              const { return fEnergyDeposit;       }
   int    NumberIonizationElectrons()  const { return fNumIonElectrons;     }
   int    NumberScintillationPhotons() const { return fNumScintPhotons;     }

   double EFieldAtStep(double fEfield, sim::SimEnergyDeposit const& edep); //value of field with any corrections for this step  
   double EFieldAtStep(double fEfield, sim::SimEnergyDeposit const& edep, std::vector<double> const&); //value of field with any corrections for this step  

 private:

   double                fEfield;              ///< value of electric field from LArProperties service
   double 	   	 fGeVToElectrons;      ///< conversion factor from LArProperties service	  
   double 	   	 fRecombA;             ///< from LArG4Parameters service			  
   double 	   	 fRecombk;             ///< from LArG4Parameters service			  
   double 	   	 fModBoxA;             ///< from LArG4Parameters service			  
   double 	   	 fModBoxB;             ///< from LArG4Parameters service			  
   bool   	   	 fUseModBoxRecomb;     ///< from LArG4Parameters service			  
   bool   	   	 fScintByParticleType; ///< from LArProperties service			  
   double 	   	 fScintYieldFactor;    ///< scintillation yield factor                             

   double fEnergyDeposit;   ///< total energy deposited in the step
   int    fNumIonElectrons; ///< number of ionization electrons for this step
   int    fNumScintPhotons; ///< number of scintillation photons for this step   
   
   std::vector<double> fEfieldOffsets;
   double EFieldAtStep(double fEfield, float x, float y, float z);
   double EFieldAtStep(double fEfield, float x, float y, float z, std::vector<double> const&);

   const detinfo::LArProperties*      fLArProp;   
   const spacecharge::SpaceCharge*    fSCE;
   
   void CalculateIonization(float e, float ds, float x, float y, float z);
   void CalculateIonization(sim::SimEnergyDeposit const& edep);

   void CalculateIonization(float e, float ds, float x, float y, float z, std::vector<double> const&);
   void CalculateIonization(sim::SimEnergyDeposit const& edep, std::vector<double> const&);

   void CalculateScintillation(float e, int pdg);
   void CalculateScintillation(sim::SimEnergyDeposit const& edep);
   
 };
}
#endif // LARG4_ISCALCULATIONSEPARATE_H

