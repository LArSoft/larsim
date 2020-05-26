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

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/IonizationScintillation/ISCalcCorrelated.h"

namespace larg4
{
    //----------------------------------------------------------------------------
    ISCalcCorrelated::ISCalcCorrelated()
    {
        std::cout << "IonizationAndScintillation/ISCalcCorrelated Initialize." << std::endl;
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
        
        // ionization work function
        fWion             = 1./fGeVToElectrons * 1e3; // MeV

        // ion+excitation work function (\todo: get from LArG4Parameters or LArProperties?)
        fWph              = 19.5 * 1e-6; // MeV

    }

    //----------------------------------------------------------------------------
    void ISCalcCorrelated::Reset()
    {
        fEnergyDeposit           = 0.;
        fNumScintPhotons         = 0.;
        fNumIonElectrons         = 0.;
        fScintillationYieldRatio = 0.;

        return;
    }

    //----------------------------------------------------------------------------
    void ISCalcCorrelated::CalcIonAndScint(sim::SimEnergyDeposit const& edep)
    {   
        fEnergyDeposit    = edep.Energy();

        // calculate total quanta (ions + excitons)
        double Nq = fEnergyDeposit / fWph;  

        float  ds         = edep.StepLength();
        double dEdx       = (ds<=0.0)? 0.0: fEnergyDeposit/ds;
        double EFieldStep = EFieldAtStep(fDetProp->Efield(), edep);
        double recomb     = 0.;

        // Guard against spurious values of dE/dx. Note: assumes density of LAr
        if(dEdx < 1.) dEdx = 1.;

        // calculate recombination survival fraction
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

        // using this recombination, calculate number of ionization electrons
        fNumIonElectrons = ( fEnergyDeposit / fWion ) * recomb;

        // calculate scintillation photons
        fNumScintPhotons = Nq - fNumIonElectrons;
    
        MF_LOG_DEBUG("ISCalcCorrelated")  << " Electrons produced for " << fEnergyDeposit
                                << " MeV deposited with "     << recomb
                                << " recombination: "         << fNumIonElectrons
                                << std::endl;
        MF_LOG_DEBUG("ISCalcCorrelated") << "number photons: " << fNumScintPhotons;

        return;

    }


    //----------------------------------------------------------------------------
    double ISCalcCorrelated::GetScintYieldRatio(sim::SimEnergyDeposit const& edep)
    {
        // For ISCalcCorrelated, the ScintByParticleType option only controls 
        // the scintillation yield ratio, which is the ratio of fast light (singlet
        // component) to the total light (singlet+triplet components).
        //
        if (!fLArProp->ScintByParticleType()) 
          return fLArProp->ScintYieldRatio();
        
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
    
    //----------------------------------------------------------------------------
    double ISCalcCorrelated::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
    {
      if (!fSCE->EnableSimEfieldSCE()) return efield;
      auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());
      return efield * std::hypot(1+eFieldOffsets.X(), eFieldOffsets.Y(), eFieldOffsets.Z());
    }

}// namespace
