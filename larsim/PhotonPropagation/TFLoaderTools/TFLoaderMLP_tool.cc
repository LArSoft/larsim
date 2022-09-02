////////////////////////////////////////////////////////////////////////
// Class:       TFLoaderMLP
// Plugin Type: tool
// File:        TFLoaderMLP_tool.cc TFLoaderMLP.h
// Aug. 20, 2022 by Mu Wei (wmu@fnal.gov)
////////////////////////////////////////////////////////////////////////
#include "larsim/PhotonPropagation/TFLoaderTools/TFLoaderMLP.h"

namespace phot
{
    //......................................................................    
    TFLoaderMLP::TFLoaderMLP(fhicl::ParameterSet const& pset)
    : ModelName{pset.get<std::string>("ModelName")}
    , InputsName{pset.get<std::vector<std::string>>("InputsName")}
    , OutputName{pset.get<std::string>("OutputName")}
    {
    }
    
    //......................................................................    
    void TFLoaderMLP::Initialization()
    {
        int num_input = int(InputsName.size());
        if(num_input != 3)
        {
            std::cout << "Input name error! exit!" << std::endl;
            return;
        }
        std::cout << "Loading TF Model from: " << ModelName << ", Input Layer: ";
        for(int i = 0; i < num_input; ++ i)
        {
            std::cout << InputsName[i] << " ";
        }
        std::cout << ", Output Layer: " << OutputName << "\n";

        //Load SavedModel
        modelbundle    = new tensorflow::SavedModelBundleLite();

        status = tensorflow::LoadSavedModel(tensorflow::SessionOptions(),
                                            tensorflow::RunOptions(),
                                            ModelName,
                                            {tensorflow::kSavedModelTagServe},
                                            modelbundle);
        
        //Initialize a tensorflow session
//        status = tensorflow::NewSession(tensorflow::SessionOptions(), &session);
        if (!status.ok())
        {
            std::cout << "Failed to load SavedModel, status: " << status.ToString() << std::endl;
            return;
        }
        
        //Read in the protobuf graph
//        tensorflow::GraphDef graph_def;
//        status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), ModelName, &graph_def);
//        if (!status.ok())
//        {
//            std::cout << status.ToString() << std::endl;
//            return;
//        }
//        
//        //Add the graph to the session
//        status = session->Create(graph_def);
//        if (!status.ok())
//        {
//            std::cout << status.ToString() << std::endl;
//            return;
//        }
        
        std::cout << "TF SavedModel loaded successfully." << std::endl;
        return;
    }
    
    //......................................................................    
    void TFLoaderMLP::CloseSession()
    {
        if (status.ok())
        {
            std::cout << "Close TF session." << std::endl;
//            session->Close();
        }

        delete modelbundle;
        
//        delete session;
        return;
    }
    
    //......................................................................    
    void TFLoaderMLP::Predict(std::vector<double> pars)
    {
        //std::cout << "TFLoader MLP:: Predicting... " << std::endl;
        int num_input = int(pars.size());
        if(num_input != 3)
        {
            std::cout << "Input parameter error! exit!" << std::endl;
            return;
        }
        //Clean prediction
        std::vector<double>().swap(prediction);

        //Define inputs
        tensorflow::Tensor pos_x(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 1}));
        tensorflow::Tensor pos_y(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 1}));
        tensorflow::Tensor pos_z(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 1}));
        auto dst_x = pos_x.flat<float>().data();
        auto dst_y = pos_y.flat<float>().data();
        auto dst_z = pos_z.flat<float>().data();
        copy_n(pars.begin(),   1, dst_x);
        copy_n(pars.begin()+1, 1, dst_y);
        copy_n(pars.begin()+2, 1, dst_z);
        std::vector<std::pair<std::string, tensorflow::Tensor>> inputs = {{InputsName[0], pos_x}, {InputsName[1], pos_y}, {InputsName[2], pos_z}};
        //Define outps
        std::vector<tensorflow::Tensor> outputs;
        
        //Run the session
        status = modelbundle->GetSession()->Run(inputs, {OutputName}, {}, &outputs);
//        status = session->Run(inputs, {OutputName}, {}, &outputs);
        if (!status.ok())
        {
            std::cout << status.ToString() << std::endl;
            return;
        }
        
        //Grab the outputs        
        unsigned int pdr = outputs[0].shape().dim_size(1);        
        //std::cout << "TFLoader MLP::Num of optical channels: " << pdr << std::endl;
        
        for (unsigned int i = 0; i < pdr; i++)
        {
            double value = outputs[0].flat<float>()(i);
            //std::cout << value << ", ";
            prediction.push_back(value);
        }
        //std::cout << std::endl;
        return;
    }
}
DEFINE_ART_CLASS_TOOL(phot::TFLoaderMLP)


