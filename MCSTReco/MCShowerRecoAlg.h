#ifndef MCSHOWERRECOALG_H
#define MCSHOWERRECOALG_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/Geometry.h"

#include "MCRecoPart.h"
#include "MCShowerRecoPart.h"
#include "MCRecoEdep.h"
#include "MCBase/MCShower.h"

// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCShowerRecoAlg {

  public:

    /// Default constructor with fhicl parameters
    MCShowerRecoAlg(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCShowerRecoAlg(){};

    void Reconstruct(const MCRecoPart& part_v,const MCRecoEdep& edep_v);

    const std::vector<sim::MCShower>& MCShower() const { return fMCShower; }

  protected:

    MCShowerRecoPart fPartAlg;

    bool             fDebugMode;
    
    std::vector<sim::MCShower> fMCShower;

    double fMinShowerEnergy;
    unsigned int fMinNumDaughters;

  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
