////////////////////////////////////////////////////////////////////////
// Class:       ISCalcNESTLAr
// Plugin Type: Algorithm
// File:        ISCalcNESTLAr.cxx
// Description:
// Aug. 30 by Mu Wei
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"

namespace larg4
{
    //----------------------------------------------------------------------------
    ISCalcNESTLAr::ISCalcNESTLAr(CLHEP::HepRandomEngine& Engine)
    : fEngine (Engine)
    {
    }
    //----------------------------------------------------------------------------
    ISCalcNESTLAr::~ISCalcNESTLAr()
    {
    }
    
    //----------------------------------------------------------------------------
    void ISCalcNESTLAr::Initialize()
    {
        std::cout << "ISCalcNESTLAr Initialize." << std::endl;
        
        fSCE       = lar::providerFrom<spacecharge::SpaceChargeService>();
        fLArProp   = lar::providerFrom<detinfo::LArPropertiesService>();
        fDetProp   = lar::providerFrom<detinfo::DetectorPropertiesService>();
        
        fScintYield        = 1.0 / (19.5*CLHEP::eV);
        fResolutionScale   = 0.107;                   // Doke 1976        
        
        return;
    }
    
    //----------------------------------------------------------------------------
    void ISCalcNESTLAr::Reset()
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
    void ISCalcNESTLAr::CalcIonAndScint(sim::SimEnergyDeposit const& edep)
    {
        CLHEP::RandGauss GaussGen(fEngine);
        CLHEP::RandFlat  UniformGen(fEngine);
        
        fNumIonElectrons  = 0.0;
        fNumScintPhotons  = 0.0;
        fYieldFactor      = 1.0;              // default, for electronic recoils 
        fExcitationRatio  = 0.21;             // ratio for light particle in LAr, such as e-, mu-, Aprile et. al book
        
        fEnergyDeposit    = edep.Energy();
        if( fEnergyDeposit < 1*CLHEP::eV ) // too small energy deposition
        {
            return;
        }
        
        int pdgcode = edep.PdgCode();
        
        geo::Length_t startx = edep.StartX();
        geo::Length_t starty = edep.StartY();
        geo::Length_t startz = edep.StartZ();
        geo::Length_t endx   = edep.EndX();
        geo::Length_t endy   = edep.EndY();
        geo::Length_t endz   = edep.EndZ();
//        geo::Point_t  midpos = edep.MidPoint();
                
//        double delta    = 1*CLHEP::mm;
//        double R0           = 1.568*CLHEP::um;    //Mozumder 1995

        double DokeBirks[3];
//        double ThomasImel = 0.00;

        double eField = EFieldAtStep(fDetProp->Efield(), edep);
        if(eField) 
        {
//            ThomasImel   = 0.156977 * pow(eField, -0.1);
            DokeBirks[0] = 0.07 * pow((eField/1.0e3), -0.85);
            DokeBirks[2] = 0.00;
        }
        else 
        {
//            ThomasImel   = 0.099;
            DokeBirks[0] = 0.0003;
            DokeBirks[2] = 0.75;
        }
        
        double Density = fDetProp->Density()/(CLHEP::g/CLHEP::cm3); // argon density at the temperature from Temperature()
        
        //biExc = 0.6;
      
        // nuclear recoil quenching "L" factor: total yield is
        // reduced for nuclear recoil as per Lindhard theory
        double epsilon = 11.5 *(fEnergyDeposit/CLHEP::keV) * pow(LAr_Z,(-7./3.));
//      double gamma   = 3. * pow(epsilon, 0.15) + 0.7 * pow( epsilon, 0.6 ) + epsilon;
//      double kappa   = 0.133*pow(z1,(2./3.))*pow(a2,(-1./2.))*(2./3.);
        
        if ( pdgcode == 2112 || pdgcode == -2112 ) //nuclear recoil
        {
//            fYieldFactor = (kappa * gamma) / ( 1 + kappa * gamma ); //Lindhard factor
            fYieldFactor = 0.23 * ( 1 + exp( -5 * epsilon ) ); //liquid argon L_eff
//            if ( eField == 0 ) 
//            {
//                ThomasImel = 0.25; //special TIB parameters for nuclear recoil only, in LAr
//            } 
            fExcitationRatio = 0.69337 + 0.3065 * exp( -0.008806 *pow( eField, 0.76313 ) );
        }
        
        // determine ultimate number of quanta from current E-deposition (ph+e-) total mean number of exc/ions
        //the total number of either quanta produced is equal to product of the
        //work function, the energy deposited, and yield reduction, for NR
        double MeanNumQuanta = fScintYield * fEnergyDeposit;      
        double sigma         = sqrt(fResolutionScale * MeanNumQuanta); //Fano
        int NumQuanta        = int(floor(GaussGen.fire(MeanNumQuanta,sigma)+0.5));
        double LeffVar       = GaussGen.fire(fYieldFactor, 0.25 * fYieldFactor);
        if (LeffVar > 1)
        {
            LeffVar = 1.0; 
        }
        
        if (LeffVar < 0)
        {
            LeffVar = 0.0; 
        }
        
        if ( fYieldFactor < 1 ) //nuclear reocils
        {
            NumQuanta = BinomFluct(NumQuanta, LeffVar);
        }
        
        //if Edep below work function, can't make any quanta, and if NumQuanta
        //less than zero because Gaussian fluctuated low, update to zero
        if(fEnergyDeposit < 1/fScintYield || NumQuanta < 0)
        {
            NumQuanta = 0;
        }
      
        // next section binomially assigns quanta to excitons and ions
        int NumExcitons = BinomFluct(NumQuanta, fExcitationRatio/(1 + fExcitationRatio));
        int NumIons     = NumQuanta - NumExcitons;
        
        // this section calculates recombination following the modified Birks'Law of Doke, deposition by deposition,
        // may be overridden later in code if a low enough energy necessitates switching to the
        // Thomas-Imel box model for recombination instead (determined by site)
        double dE  = fEnergyDeposit/CLHEP::MeV;
        double dx  = 0.0;
        double LET = 0.0;
        double recombProb;
        
        if ( pdgcode != 11 && pdgcode != -11 &&  pdgcode != 13 && pdgcode != -13 ) //e-: 11, e+: -11, mu-: 13, mu+: -13
        {
            //in other words, if it's a gamma,ion,proton,alpha,pion,et al. do not
            //use the step length provided by Geant4 because it's not relevant,
            //instead calculate an estimated LET and range of the electrons that
            //would have been produced if Geant4 could track them
            LET = CalcElectronLET( 1000*dE );
            
            if(LET)
            {
                dx = dE/(Density * LET); //find the range based on the LET
            }
            
            if(abs(pdgcode) == 2112)  //nuclear recoils
            {
                dx = 0;
            }
        }
        else //normal case of an e-/+ energy deposition recorded by Geant
        { 
            dx = std::sqrt((startx-endx)*(startx-endx)+(starty-endy)*(starty-endy)+(startz-endz)*(startz-endz))/CLHEP::cm;
            if(dx)
            {
                LET = (dE/dx)*(1/Density); //lin. energy xfer (prop. to dE/dx)
            }
            if ( LET > 0 && dE > 0 && dx > 0 ) 
            {
                double ratio = CalcElectronLET( dE*1e3 )/LET;
                if ( ratio < 0.7 && pdgcode == 11 ) 
                {
                    dx  /= ratio;
                    LET *= ratio; 
                }
            }
        }
        
        DokeBirks[1] = DokeBirks[0]/(1-DokeBirks[2]);                        //B=A/(1-C) (see paper)    r
        recombProb   = (DokeBirks[0]*LET)/(1+DokeBirks[1]*LET)+DokeBirks[2]; //Doke/Birks' Law as spelled out in the NEST pape
        recombProb  *= (Density/Density_LAr);
        
        //check against unphysicality resulting from rounding errors
        if(recombProb < 0)
        {
            recombProb = 0;
        }
        if(recombProb > 1)
        {
            recombProb = 1;
        }
        
        //use binomial distribution to assign photons, electrons, where photons
        //are excitons plus recombined ionization electrons, while final
        //collected electrons are the "escape" (non-recombined) electrons
        int NumPhotons   = NumExcitons + BinomFluct(NumIons,recombProb);
        int NumElectrons = NumQuanta - NumPhotons;
        
        fNumIonElectrons          = NumElectrons;
        fNumScintPhotons          = NumPhotons;
        fScintillationYieldRatio  = GetScintYieldRatio(edep);
//        fNumFastScintPhotons      = fNumScintPhotons * fScintillationYieldRatio;
//        fNumSlowScintPhotons      = fNumScintPhotons - fNumFastScintPhotons;
        
        return;
    }    

    //----------------------------------------------------------------------------
    int ISCalcNESTLAr::BinomFluct ( int N0, double prob ) 
    {
        CLHEP::RandGauss GaussGen(fEngine);
        CLHEP::RandFlat  UniformGen(fEngine);
        
        double mean  = N0*prob;
        double sigma = sqrt(N0*prob*(1-prob));
        int N1       = 0;
        if ( prob == 0.00 )
        {
            return N1;
        }
        if ( prob == 1.00 )
        {
            return N0;
        }
        
        if ( N0 < 10 ) 
        {
            for(int i = 0; i < N0; i++) 
            {
                if(UniformGen.fire() < prob)
                {
                    N1++;
                }
            }
        }
        else 
        {
             N1 = int(floor(GaussGen.fire(mean,sigma)+0.5));
        }
        if ( N1 > N0 )
        {
            N1 = N0;
        }
        if ( N1 < 0 )
        {
            N1 = 0;
        }
        return N1;
    }
    
    //----------------------------------------------------------------------------
    double ISCalcNESTLAr::CalcElectronLET ( double E ) 
    {
        double LET;
        
        if ( E >= 1 )
        {
            LET = 116.70 - 162.97*log10(E) + 99.361*pow(log10(E),2)
                - 33.405*pow(log10(E),3) + 6.5069*pow(log10(E),4)
                - 0.69334*pow(log10(E),5) +.031563*pow(log10(E),6);
        }
        else if ( E>0 && E<1 )
        {
            LET = 100;
        }
        else
        {
            LET = 0;
        }
        
        return LET;
    }
    
    //----------------------------------------------------------------------------
    double ISCalcNESTLAr::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
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
    
    //......................................................................    
    double ISCalcNESTLAr::GetScintYieldRatio(sim::SimEnergyDeposit const& edep)
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
}
