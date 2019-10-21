////////////////////////////////////////////////////////////////////////
// Class:       ISCalcNESTLAr
// Plugin Type: Algorithm
// File:        ISCalcNESTLAr.cxx
// Description:
// Aug. 30 by Mu Wei
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_ISCALCNESTLAr_H
#define LARG4_ISCALCNESTLAr_H

#include "larsim/IonizationScintillation/ISCalc.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CLHEP/Units/SystemOfUnits.h"
// Random number engine
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#define LAr_Z       (18)
#define Density_LAr (1.393)
#define LAr_W_Value (19.5*CLHEP::eV)

namespace larg4
{
    class ISCalcNESTLAr : public ISCalc
    {
    public:
        explicit ISCalcNESTLAr(CLHEP::HepRandomEngine& fEngine);
        virtual ~ISCalcNESTLAr();
        void   Initialize();
        void   Reset();
        
        double EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep); //value of field with any corrections for this step
        void   CalcIonAndScint(sim::SimEnergyDeposit const& edep);
        
    private:
        double  fResolutionScale;   
        double  fScintYield;             // quanta (electrons or photons) yield per eV
        double  fYieldFactor;            // quenching factor
        double  fExcitationRatio;        // N_ex/N_i, the dimensionless ratio of initial excitons to ions
        
        const spacecharge::SpaceCharge*    fSCE;
        const detinfo::DetectorProperties* fDetProp;
        const detinfo::LArProperties*      fLArProp;
        CLHEP::HepRandomEngine&            fEngine; // random engine

        int    BinomFluct (int N0, double prob );
        double CalcElectronLET (double E );
        double GetScintYieldRatio(sim::SimEnergyDeposit const& edep);
    };
}
#endif // LARG4_ISCALCNESTLAr_H
