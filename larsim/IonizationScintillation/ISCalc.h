////////////////////////////////////////////////////////////////////////
// Class:       ISCalc
// Plugin Type: Algorithm
// File:        ISCalc.h and ISCalc.cxx
// Description: Interface to the ion and scintillation calculation algorithm
// Input: 'sim::SimEnergyDeposit'
// Output: number of photons and electrons
// Spet.16 by Mu Wei
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_ISCALC_H
#define LARG4_ISCALC_H

#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace larg4
{
    class ISCalc
    {
    public:
        ISCalc();
        virtual ~ISCalc();
        virtual void       Initialize()                                                         = 0;
        virtual void       Reset()                                                              = 0;
        virtual void       CalcIonAndScint(sim::SimEnergyDeposit const& edep)                   = 0;
        virtual double     EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)       = 0; //value of field with any corrections for this step
        
        double             EnergyDeposit() const                       {return fEnergyDeposit;}
        double             NumOfElectrons() const                      {return fNumIonElectrons;}
        double             NumOfPhotons() const                        {return fNumScintPhotons;}
//        double             NumOfFastPhotons() const                    {return fNumFastScintPhotons;}
//        double             NumOfSlowPhotons() const                    {return fNumSlowScintPhotons;}
        double             ScintillationYieldRatio() const             {return fScintillationYieldRatio;}
        
    protected:
        double  fEnergyDeposit;             // total energy deposited in the step
        double  fNumIonElectrons;           // number of ionization electrons for this step
        double  fNumScintPhotons;           // number of scintillation photons for this step
//        double  fNumFastScintPhotons;       // number of fast scintillation photons for this step
//        double  fNumSlowScintPhotons;       // number of slow scintillation photons for this step
        double  fScintillationYieldRatio;   // liquid argon scintillation yield ratio
    };
}
#endif // LARG4_ISCALC_H
