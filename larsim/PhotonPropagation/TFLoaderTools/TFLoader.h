////////////////////////////////////////////////////////////////////////
// Class:       TFLoader
// Plugin Type: tool
// File:        TFLoader.cc TFLoader.h
// Aug. 20, 2022 by Mu Wei (wmu@fnal.gov)
////////////////////////////////////////////////////////////////////////
#ifndef TFLoader_H
#define TFLoader_H

#include <memory>
#include <vector>
#include <string>

//#include "tensorflow/core/public/session.h"
//#include "tensorflow/core/platform/env.h"

#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"

//namespace tensorflow
//{
//    class Session;
//    class Tensor;
//    struct SavedModelBundle;
//}

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

