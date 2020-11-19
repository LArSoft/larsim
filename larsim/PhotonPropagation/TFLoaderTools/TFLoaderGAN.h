////////////////////////////////////////////////////////////////////////
// Class:       TFLoaderGAN
// Plugin Type: tool
// File:        TFLoaderGAN_tool.cc TFLoaderGAN.h
// Oct. 21, 2020 by Mu Wei (wmu@fnal.gov)
////////////////////////////////////////////////////////////////////////
#ifndef TFLoaderGAN_H
#define TFLoaderGAN_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/PhotonPropagation/TFLoaderTools/TFLoader.h"


namespace phot
{
    class TFLoaderGAN : public TFLoader
    {
    public:
        explicit TFLoaderGAN(fhicl::ParameterSet const& pset);
        void Initialization();
        void CloseSession();
        void Predict(std::vector<double> pars)  ;
        
    private:
        std::string           ModelName;    //Full path to the model (.pb) file;
        std::string           InputsName;   //Name of the input layer;
        std::string           OutputName;   //Name of the output layer;
        
        tensorflow::Session*  session;
        tensorflow::Status    status;
    };
}
#endif
