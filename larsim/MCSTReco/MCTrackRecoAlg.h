#ifndef MCTRACKRECOALG_H
#define MCTRACKRECOALG_H

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
#include "MCRecoEdep.h"
#include "lardataobj/MCBase/MCTrack.h"

// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCTrackRecoAlg {

  public:

    /// Default constructor with fhicl parameters
    explicit MCTrackRecoAlg(fhicl::ParameterSet const& pset);
    std::unique_ptr<std::vector<sim::MCTrack>> Reconstruct(MCRecoPart& part_v, MCRecoEdep& edep_v);

  protected:
    bool             fDebugMode;

  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
