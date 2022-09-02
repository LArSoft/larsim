////////////////////////////////////////////////////////////////////////
// Class:       TFLoaderMLP
// Plugin Type: tool
// File:        TFLoaderMLP_tool.cc TFLoaderMLP.h
// Aug. 20, 2022 by Mu Wei (wmu@fnal.gov)
////////////////////////////////////////////////////////////////////////
#ifndef TFLoaderMLP_H
#define TFLoaderMLP_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/PhotonPropagation/TFLoaderTools/TFLoader.h"

namespace phot
{
    class TFLoaderMLP : public TFLoader
    {
    public:
        explicit TFLoaderMLP(fhicl::ParameterSet const& pset);
        void Initialization();
        void CloseSession();
        void Predict(std::vector<double> pars)  ;
        
    private:
        std::string              ModelName;    //Full path to the model (.pb) file;
        std::vector<std::string> InputsName;   //Name of the input layer;
        std::string              OutputName;   //Name of the output layer;
        
        tensorflow::SavedModelBundleLite*  modelbundle;
//        tensorflow::SessionOptions*        sessionoptions;
//        tensorflow::Session*               session;
        tensorflow::Status                 status;
    };
}
#endif



