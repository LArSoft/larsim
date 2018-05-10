#ifndef WEIGHTMANAGER_CXX
#define WEIGHTMANAGER_CXX

#include <sstream>
#include <map>
#include <set>
#include "WeightManager.h"


namespace evwgh {

  WeightManager::WeightManager(const std::string name)
    : _name(name)
  {
    _configured = false;
  }

  const std::string& WeightManager::Name() const
  { return _name; }

  size_t WeightManager::Configure(fhicl::ParameterSet const & p, art::EngineCreator &engine) 
  {
   
    ::art::ServiceHandle<rndm::NuRandomService> seedservice;

    // Get list of weight functions
    std::vector<std::string> rw_func = p.get<std::vector<std::string>>("weight_functions");

    // Loop over all the functions and register them
    for (auto ifunc=rw_func.begin(); ifunc!=rw_func.end(); ifunc++) 
    {
      fhicl::ParameterSet const &ps_func = p.get<fhicl::ParameterSet> (*ifunc);
      std::string func_type = ps_func.get<std::string>("type");
  
      WeightCalc* wcalc=WeightCalcFactory::Create(func_type+"WeightCalc");
      if ( wcalc == NULL )
        throw cet::exception(__FUNCTION__) << "Function " << *ifunc << " requested in fcl file has not been registered!" << std::endl;
      if ( fWeightCalcMap.find(*ifunc)!= fWeightCalcMap.end() )
        throw cet::exception(__FUNCTION__) << "Function " << *ifunc << " has been requested multiple times in fcl file!" << std::endl;
  
      //mf::LogInfo("") << "Configuring weight calculator " << *ifunc;
  
      // Create random engine for each rw function (name=*ifunc) (and seed it with random_seed set in the fcl)
      seedservice->createEngine(engine, "HepJamesRandom", *ifunc, ps_func, "random_seed");
  
      wcalc->SetName(*ifunc);
      wcalc->Configure(p);
      Weight_t* winfo=new Weight_t();
      winfo->fWeightCalcType=func_type;
      winfo->fWeightCalc=wcalc;
      winfo->fNmultisims=ps_func.get<int>("number_of_multisims");
  
      std::pair<std::string, Weight_t*> pwc(*ifunc,winfo);
      fWeightCalcMap.insert(pwc);
    }

    _configured = true;

    return fWeightCalcMap.size();
  }


 
  // 
  // CORE FUNCTION
  //
  MCEventWeight WeightManager::Run(art::Event & e, const int inu)
  {     

    if (!_configured)
      throw cet::exception(__PRETTY_FUNCTION__) << "Have not configured yet!" << std::endl;

    //
    // Loop over all functions ang calculate weights
    //
    MCEventWeight mcwgh;
    for (auto it = fWeightCalcMap.begin() ;it != fWeightCalcMap.end(); it++) {
      
      auto const & weights = it->second->GetWeight(e);
      
      if(weights.size() == 0){
        std::vector<double> empty;
        std::pair<std::string, std::vector <double> > p("empty",empty);
        mcwgh.fWeight.insert(p);
      } 
      else{
        std::pair<std::string, std::vector<double> > 
          p(it->first+"_"+it->second->fWeightCalcType,
            weights[inu]);
        mcwgh.fWeight.insert(p);
      }
    }

    return mcwgh;
  }



  void WeightManager::PrintConfig() {
    
    return; 
  }

}

#endif
