#ifndef MCSHOWERRECOALG_H
#define MCSHOWERRECOALG_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"

#include "MCRecoPart.h"
#include "MCShowerRecoPart.h"
#include "MCRecoEdep.h"
#include "lardataobj/MCBase/MCShower.h"

// STL
#include <set>
#include <vector>

#include <sstream>
#include <memory>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCShowerRecoAlg {

  public:

    /// Default constructor with fhicl parameters
    explicit MCShowerRecoAlg(fhicl::ParameterSet const& pset);

    std::unique_ptr<std::vector<sim::MCShower>> Reconstruct(MCRecoPart& part_v,MCRecoEdep& edep_v);

  protected:

    MCShowerRecoPart fPartAlg;
    bool             fDebugMode;
    double fMinShowerEnergy;
    unsigned int fMinNumDaughters;

  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
