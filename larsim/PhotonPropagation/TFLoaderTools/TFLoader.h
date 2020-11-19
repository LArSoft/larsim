////////////////////////////////////////////////////////////////////////
// Class:       TFLoader
// Plugin Type: tool
// File:        TFLoader_tool.cc TFLoader.h
// Oct. 21, 2020 by Mu Wei (wmu@fnal.gov)
////////////////////////////////////////////////////////////////////////
#ifndef TFLoader_H
#define TFLoader_H

#include <memory>
#include <vector>
#include <string>

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

namespace phot
{
    class TFLoader
    {
    public:
        TFLoader();
        virtual ~TFLoader() = default;

        virtual void Initialization()       = 0;
        virtual void CloseSession()      = 0;
        virtual void Predict(std::vector<double> pars)    = 0;
        std::vector<double> GetPrediction() const   {return prediction;}
        
    protected:
        std::vector<double>   prediction;
    };
}
#endif
