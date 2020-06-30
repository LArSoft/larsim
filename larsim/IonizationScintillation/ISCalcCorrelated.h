////////////////////////////////////////////////////////////////////////
// Class:       ISCalcCorrelated
// Plugin Type: algorithm
// File:        ISCalcCorrelated.h and ISCalcCorrelated.cxx
// Description: Interface to algorithm class for a specific calculation of
//              ionization electrons and scintillation photons, based on 
//              simple microphysics arguments to establish an anticorrelation
//              between these two quantities.
// Input: 'sim::SimEnergyDeposit'
// Output: num of Photons and Electrons
// May 2020 by W Foreman
////////////////////////////////////////////////////////////////////////

#ifndef IS_ISCALCCORRELATED_H
#define IS_ISCALCCORRELATED_H

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
    class ISCalcCorrelated : public ISCalc
    {
    public:
        ISCalcCorrelated();
        void   Reset();
        
        double EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep); //value of field with any corrections for this step
        void   CalcIonAndScint(sim::SimEnergyDeposit const& edep);
        
    private:
        double  fGeVToElectrons;      ///< from LArG4Parameters service
        double  fWion;                ///< W_ion (23.6 eV) == 1/fGeVToElectrons
        double  fWph;                 ///< W_ph (19.5 eV) 
        double  fScintPreScale;       ///< scintillation pre-scaling factor from LArProperties service
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
#endif
