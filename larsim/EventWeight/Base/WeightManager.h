/**
 * \file WeightManager.h
 *
 * 
 * \brief Allows to interface to EventWeight calculators
 *
 * @author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 */

#ifndef WEIGHTMANAGER_H
#define WEIGHTMANAGER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Utilities/InputTag.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "lardataobj/Simulation/sim.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Weight_t.h"
#include "MCEventWeight.h"
#include "WeightCalc.h"
#include "WeightCalcFactory.h"

namespace evwgh {
  /**
     \class WeightManager
  */
  class WeightManager {

  public:
    
    /// Default constructor
    WeightManager(const std::string name="WeightManager");
    
    /// Default destructor
    ~WeightManager(){}

    /// Name getter
    const std::string& Name() const;

    /** 
      * @brief Configuration function
      * @param cfg the input parameters for settings
      * @param the enging creator for the random seed (usually passed with *this)
       CONFIGURE FUNCTION: created the weights algorithms in the following way: \n
       0) Looks at the weight_functions fcl parameter to get the name of the calculators   \n
       1) Creates the Calculators requested in step 0, and assigne a different random seed to each one \n
       3) The future call WeightManager::Run will run the calculators           \n
    */
    template <typename Module>
    size_t Configure(fhicl::ParameterSet const & cfg, Module& module);

    /**
      * @brief Core function (previous call to Configure is needed)
      * @param e the art event
      * @param inu the index of the simulated neutrino in the event
       CORE FUNCTION: executes algorithms to assign a weight to the event as requested users. \n
       WeightManager::Configure needs to be called first \n
       The execution takes following steps:             \n
       0) Loos over all the previously emplaced calculators \n
       1) For each of them calculates the weights (more weight can be requested per calculator) \n
       3) Returns a map from "calculator name" to vector of weights calculated which is available inside MCEventWeight
     */
    MCEventWeight Run(art::Event &e, const int inu);

    /**
      * @brief Returns the map between calculator name and Weight_t product
      */
    std::map<std::string, Weight_t*> GetWeightCalcMap() { return fWeightCalcMap; }

    /// Reset
    void Reset()
    { _configured = false; }

    void PrintConfig();


  private:

    std::map<std::string, Weight_t*> fWeightCalcMap; ///< A set of custom weight calculators

    bool _configured = false; ///< Readiness flag

    std::string _name; ///< Name

  };

  template <typename Module>
  size_t WeightManager::Configure(fhicl::ParameterSet const & p, Module& module)
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
      seedservice->createEngine(module, "HepJamesRandom", *ifunc, ps_func, "random_seed");

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

}

#endif
