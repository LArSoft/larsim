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

#include "larsim/IonizationScintillation/ISCalcSeparate.h"

namespace larg4
{
    //----------------------------------------------------------------------------
    ISCalcSeparate::ISCalcSeparate()
    {
    }
    
    //----------------------------------------------------------------------------
    ISCalcSeparate::~ISCalcSeparate()
    {
    }
    
    //----------------------------------------------------------------------------
    void ISCalcSeparate::Initialize()
    {
        std::cout << "ISCalcSeparate Initialize." << std::endl;
        art::ServiceHandle<sim::LArG4Parameters const> LArG4PropHandle;
        
        fSCE       = lar::providerFrom<spacecharge::SpaceChargeService>();
        fLArProp   = lar::providerFrom<detinfo::LArPropertiesService>();
        fDetProp   = lar::providerFrom<detinfo::DetectorPropertiesService>();
        
        fScintYieldFactor  = 1.; // true scintillation yield will be got from LArProperties
        
        //the recombination coefficient is in g/(MeVcm^2), but we report energy depositions in MeV/cm,
        //need to divide Recombk from the LArG4Parameters service by the density of the argon we got above.
        fRecombA          = LArG4PropHandle->RecombA();
        fRecombk          = LArG4PropHandle->Recombk()/fDetProp->Density(fDetProp->Temperature());
        fModBoxA          = LArG4PropHandle->ModBoxA();
        fModBoxB          = LArG4PropHandle->ModBoxB()/fDetProp->Density(fDetProp->Temperature());
        fUseModBoxRecomb  = (bool)LArG4PropHandle->UseModBoxRecomb();
        fGeVToElectrons   = LArG4PropHandle->GeVToElectrons();
        
        return;
    }
    
    //----------------------------------------------------------------------------
    void ISCalcSeparate::Reset()
    {
        fEnergyDeposit           = 0.;
        fNumScintPhotons         = 0.;
        fNumIonElectrons         = 0.;
//        fNumFastScintPhotons     = 0.;
//        fNumSlowScintPhotons     = 0.;
        fScintillationYieldRatio = 0.;
        
        return;
    }
    
    //----------------------------------------------------------------------------
    // fNumIonElectrons returns a value that is not corrected for life time effects
    void ISCalcSeparate::CalcIon(sim::SimEnergyDeposit const& edep)
    {
        float e           = edep.Energy();
        float ds          = edep.StepLength();
        
        double recomb     = 0.;
        double dEdx       = (ds<=0.0)? 0.0: e/ds;
        double EFieldStep = EFieldAtStep(fDetProp->Efield(), edep);
        
        // Guard against spurious values of dE/dx. Note: assumes density of LAr
        if(dEdx < 1.)
        {
            dEdx = 1.;
        }
        
        if(fUseModBoxRecomb) 
        {
            if(ds>0)
            {
                double Xi = fModBoxB * dEdx / EFieldStep;
                recomb    = log(fModBoxA + Xi) / Xi;
            }
            else
            {
                recomb = 0;
            }
        }
        else
        {
            recomb = fRecombA / (1. + dEdx * fRecombk / EFieldStep);
        }
        
        // 1.e-3 converts fEnergyDeposit to GeV
        fNumIonElectrons = fGeVToElectrons * 1.e-3 * e * recomb;
        
        MF_LOG_DEBUG("ISCalcSeparate")  << " Electrons produced for " << fEnergyDeposit
                                        << " MeV deposited with "     << recomb
                                        << " recombination: "         << fNumIonElectrons 
                                        << std::endl;
        return;
    }
    
    //----------------------------------------------------------------------------
    void ISCalcSeparate::CalcScint(sim::SimEnergyDeposit const& edep)
    {
        float e    = edep.Energy();
        int   pdg  = edep.PdgCode();
        
        double scintYield = fLArProp->ScintYield(true);
        if(fLArProp->ScintByParticleType())
        {
            MF_LOG_DEBUG("ISCalcSeparate") << "scintillating by particle type";
            
            switch(pdg)
            {
                case 2212:
                    scintYield = fLArProp->ProtonScintYield(true);
                    break;
                case  13:
                case -13:
                    scintYield = fLArProp->MuonScintYield(true);
                    break;
                case  211:
                case -211:
                    scintYield = fLArProp->PionScintYield(true);
                    break;
                case  321:
                case -321:
                    scintYield = fLArProp->KaonScintYield(true);
                    break;
                case 1000020040:
                    scintYield = fLArProp->AlphaScintYield(true);
                    break;
                case  11:
                case -11:
                case  22:
                    scintYield = fLArProp->ElectronScintYield(true);
                    break;
                default:
                    scintYield = fLArProp->ElectronScintYield(true);
            }
            
            fNumScintPhotons = scintYield * e;
        }
        else
        {
            fNumScintPhotons = fScintYieldFactor * scintYield * e;
        }
        
        fScintillationYieldRatio = GetScintYieldRatio(edep);
//        fNumFastScintPhotons     = fNumScintPhotons * fScintillationYieldRatio;
//        fNumSlowScintPhotons     = fNumScintPhotons - fNumFastScintPhotons;
        return;
    }
    
    //----------------------------------------------------------------------------
    void ISCalcSeparate::CalcIonAndScint(sim::SimEnergyDeposit const& edep)
    {
        fEnergyDeposit = edep.Energy();
        CalcIon(edep);
        CalcScint(edep);
        
        return;
    }
    //----------------------------------------------------------------------------
    double ISCalcSeparate::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
    {
        geo::Point_t pos = edep.MidPoint();
        double EField    = efield;
        
        geo::Vector_t eFieldOffsets;
        
        if (fSCE->EnableSimEfieldSCE())
        {
            eFieldOffsets = fSCE->GetEfieldOffsets(pos);
            EField = std::sqrt((efield + efield*eFieldOffsets.X())*(efield + efield*eFieldOffsets.X())
                              +(efield*eFieldOffsets.Y()*efield*eFieldOffsets.Y()) 
                              +(efield*eFieldOffsets.Z()*efield*eFieldOffsets.Z()) );
        }
        
        return EField;
    }
    //----------------------------------------------------------------------------
    double ISCalcSeparate::GetScintYieldRatio(sim::SimEnergyDeposit const& edep)
    {
        if (!fLArProp->ScintByParticleType())
        {
            return fLArProp->ScintYieldRatio();
        }
        switch (edep.PdgCode())
        {
            case  2212:
               return fLArProp->ProtonScintYieldRatio();
            case  13:
            case -13:
                return fLArProp->MuonScintYieldRatio();
            case  211:
            case -211:
                return fLArProp->PionScintYieldRatio();
            case  321:
            case -321:
                return fLArProp->KaonScintYieldRatio();
            case 1000020040:
                return fLArProp->AlphaScintYieldRatio();
            case  11:
            case -11:
            case  22:
                return fLArProp->ElectronScintYieldRatio();
            default:
                return fLArProp->ElectronScintYieldRatio();
        }
    }
    
}// namespace
