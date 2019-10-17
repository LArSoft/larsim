////////////////////////////////////////////////////////////////////////
// Class:       ScintTimeLAr
// Plugin Type: tool
// File:        ScintTimeLAr_tool.cc ScintTimeLAr.h
// Description:
// Generate a random number for timing of LAr scintillation
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTimeLAr.h"

namespace phot
{
    //......................................................................    
    ScintTimeLAr::ScintTimeLAr(fhicl::ParameterSet const& pset)
    : LogLevel{pset.get<int>("LogLevel")}
    , SRTime{pset.get<double>("SlowRisingTime", 0.0)}
    , SDTime{pset.get<double>("SlowDecayTime",  0.0)}
    , FRTime{pset.get<double>("FastRisingTime", 0.0)}
    , FDTime{pset.get<double>("FastDecayTime",  0.0)}
    {
        if ( LogLevel >= 1 ) 
        {
            std::cout << "ScintTimeLAr Tool configure:" << std::endl;
            std::cout << "Fast rising time: "   << FRTime
                      << ", Fast decay time: "  << FDTime
                      << ", Slow rising time: " << SRTime
                      << ", Slow decay time: "  << SDTime
                      << std::endl;
        }
    }
        
    //......................................................................    
    double ScintTimeLAr::single_exp(double t, double tau2)
    {
        return std::exp((-1.0 * t) / tau2) / tau2;
    }
    
    //......................................................................    
    double ScintTimeLAr::bi_exp(double t, double tau1, double tau2)
    {
        return (((std::exp((-1.0 * t) / tau2) * (1.0 - std::exp((-1.0 * t) / tau1))) / tau2) / tau2) * (tau1 + tau2);
    }
    
    //......................................................................    
    // Returns the time within the time distribution of the scintillation process, when the photon was created.
    // Scintillation light has an exponential decay which is given by the decay time, tau2,
    // and an exponential increase, which here is given by the rise time, tau1.
    // randflatScintTimeLAr is passed to use the saved seed from the RandomNumberSaver in order to be able to reproduce the same results.
    void ScintTimeLAr::GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine)
    {        
        double tau1;
        double tau2;
        
        if(is_fast == 1)  // fast scinitllation
        {
            tau1 = FRTime;
            tau2 = FDTime;
        }
        else
        {
            tau1 = SRTime;
            tau2 = SDTime;
        }
        
        CLHEP::RandFlat randflatscinttime{engine};
        
        if ((tau1 == 0.0) || (tau1 == -1.0))
        {
            timing = -tau2 * std::log(randflatscinttime());
            return;
        }
        
        //ran1, ran2 = random numbers for the algorithm        
        while (1)
        {
            auto ran1 = randflatscinttime();
            auto ran2 = randflatscinttime();
            auto d = (tau1 + tau2) / tau2;
            auto t = -tau2 * std::log(1 - ran1);
            auto g = d * single_exp(t, tau2);
            if (ran2 <= bi_exp(t, tau1, tau2) / g)
            {
                timing = t;
                return;
            }
        }
    }
}
    
DEFINE_ART_CLASS_TOOL(phot::ScintTimeLAr)
