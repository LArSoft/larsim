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
