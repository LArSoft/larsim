////////////////////////////////////////////////////////////////////////
// Class:       ScintTime
// Plugin Type: tool
// File:        ScintTime_tool.cc ScintTime.h
// Description:
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#ifndef ScintTime_H
#define ScintTime_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//#include <string>

namespace phot
{
    class ScintTime
    {
    public:
        ScintTime();
        virtual ~ScintTime() = default;

        virtual void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine)      = 0;
        double GetScintTime() const        {return timing;}
        
    protected:
        double timing;
    };
}

#endif
