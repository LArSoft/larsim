////////////////////////////////////////////////////////////////////////
// Class:       TFLoaderGAN
// Plugin Type: tool
// File:        TFLoaderGAN_tool.cc TFLoaderGAN.h
// Oct. 21, 2020 by Mu Wei (wmu@fnal.gov)
////////////////////////////////////////////////////////////////////////
#include "larsim/PhotonPropagation/TFLoaderTools/TFLoaderGAN.h"

namespace phot
{
    //......................................................................    
    TFLoaderGAN::TFLoaderGAN(fhicl::ParameterSet const& pset)
    : ModelName{pset.get<std::string>("ModelName")}
    , InputsName{pset.get<std::string>("InputsName")}
    , OutputName{pset.get<std::string>("OutputName")}
    {
    }
    
    //......................................................................    
    void TFLoaderGAN::Initialization()
    {
        std::cout << "Load TF Model: " << ModelName << ", Input Layer: " << InputsName << ", Output Layer: " << OutputName << "\n";
        
        //Initialize a tensorflow session
        status = tensorflow::NewSession(tensorflow::SessionOptions(), &session);
        if (!status.ok())
        {
            std::cout << status.ToString() << std::endl;
            return;
        }
        
        //Read in the protobuf graph
        tensorflow::GraphDef graph_def;
        status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), ModelName, &graph_def);
        if (!status.ok())
        {
            std::cout << status.ToString() << std::endl;
            return;
        }
        
        //Add the graph to the session
        status = session->Create(graph_def);
        if (!status.ok())
        {
            std::cout << status.ToString() << std::endl;
            return;
        }
        
        std::cout << "TFLoader tool loaded successfully" << std::endl;
        return;
    }
    
    //......................................................................    
    void TFLoaderGAN::CloseSession()
    {
        if (status.ok())
        {
            std::cout << "Close TF session" << std::endl;
            session->Close();
        }
        
        delete session;
        return;
    }
    
    //......................................................................    
    void TFLoaderGAN::Predict(std::vector<double> pars)
    {
        //Clean prediction
        std::vector<double>().swap(prediction);
        
        //Define inputs
        tensorflow::Tensor input(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 3}));
        auto dst = input.flat<float>().data();
        copy_n(pars.begin(), 3, dst);        
        std::vector<std::pair<std::string, tensorflow::Tensor>> inputs = {{ InputsName, input}};
        
        //Define outps
        std::vector<tensorflow::Tensor> outputs;
        
        //Run the session
        status = session->Run(inputs, {OutputName}, {}, &outputs);
        if (!status.ok())
        {
            std::cout << status.ToString() << std::endl;
            return;
        }
        
        //Grab the outputs        
        unsigned int img_num = outputs[0].shape().dim_size(0);
        unsigned int img_row = outputs[0].shape().dim_size(1);
        unsigned int img_col = outputs[0].shape().dim_size(2);
        unsigned int img_cha = outputs[0].shape().dim_size(3);
        
        if(img_num !=1 || img_cha != 1)
        {
            std::cout << "Num of image or channel error!" << std::endl;
        }
        
        for (unsigned int row_i = 0; row_i < img_row; row_i++)
        {
            for (unsigned int col_i = 0; col_i < img_col; col_i++)
            {
                // Get vaule through .flat()
                unsigned int offset = row_i*img_col + col_i;
                double value = outputs[0].flat<float>()(offset);
//                // Get value through .tensor()
//                double value = outputs[0].tensor<float, 4>()(0, row_i, col_i, 0);
//                std::cout << "output[0](0, " << row_i << ", " << col_i << ", 0) = " << value << std::endl;
                prediction.push_back(value);
            }
        }        
        return;
    }
}
DEFINE_ART_CLASS_TOOL(phot::TFLoaderGAN)
