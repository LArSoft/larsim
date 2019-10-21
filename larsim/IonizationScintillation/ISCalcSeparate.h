////////////////////////////////////////////////////////////////////////
// Class:       ISCalcSeparate
// Plugin Type: algorithm
// File:        ISCalcSeparate.h and ISCalcSeparate.cxx
// Description:
// Interface to algorithm class for a specific calculation of ionization electrons and scintillation photons
// assuming there is no correlation between the two
// Input: 'sim::SimEnergyDeposit'
// Output: num of Photons and Electrons
// Sept.16 by Mu Wei
////////////////////////////////////////////////////////////////////////

#ifndef IS_ISCALCSEPARATE_H
#define IS_ISCALCSEPARATE_H

#include "larsim/IonizationScintillation/ISCalc.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include <vector>
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 
{
    class ISCalcSeparate : public ISCalc
    {
    public:
        ISCalcSeparate();
        virtual ~ISCalcSeparate();
        void   Initialize();
        void   Reset();
        
        double EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep); //value of field with any corrections for this step
        void   CalcIonAndScint(sim::SimEnergyDeposit const& edep);
        
    private:
        double  fGeVToElectrons;      ///< from LArG4Parameters service
        double  fScintYieldFactor;    ///< scintillation yield factor
        double  fRecombA;             ///< from LArG4Parameters service
        double  fRecombk;             ///< from LArG4Parameters service
        double  fModBoxA;             ///< from LArG4Parameters service
        double  fModBoxB;             ///< from LArG4Parameters service
        bool    fUseModBoxRecomb;     ///< from LArG4Parameters service
        
        const spacecharge::SpaceCharge*    fSCE;
        const detinfo::DetectorProperties* fDetProp;
        const detinfo::LArProperties*      fLArProp;
        
        void CalcIon(sim::SimEnergyDeposit const& edep);
        void CalcScint(sim::SimEnergyDeposit const& edep);
        double GetScintYieldRatio(sim::SimEnergyDeposit const& edep);
    };
}
#endif // LARG4_ISCALCULATIONSEPARATE_H
