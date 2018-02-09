/**
 * \file WeightManager.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class WeightManager
 *
 * @author Marco Del Tutto
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
#include "art/Framework/Core/EngineCreator.h"
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

    /// Configuration
    size_t Configure(fhicl::ParameterSet const & cfg, art::EngineCreator &);

    /**
       CORE FUNCTION: executes algorithms to find a match of TPC object and flash provided by users. \n
       The execution takes following steps:             \n
       0) TPC filter algorithm if provided (optional)   \n
       1) Flash filter algorithm if provided (optional) \n
       3) Flash matching algorithm (required)           \n
       4) Returns match information for created TPC object & flash pair which respects the outcome of 3)
     */
    MCEventWeight Run(art::Event &e, const int inu);

    /// Clears locally kept TPC object (QClusterArray_t) and flash (FlashArray_t), both provided by a user
    void Reset()
    { _configured = false; }

    void PrintConfig();


  private:

    std::map<std::string, Weight_t*> fWeightCalcMap; ///< A set of custom weight calculators

    /// Readiness flag
    bool _configured = false;

    /// Name
    std::string _name;


  };
}

#endif

