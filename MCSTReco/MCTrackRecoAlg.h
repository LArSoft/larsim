#ifndef MCTRACKRECOALG_H
#define MCTRACKRECOALG_H

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
#include "MCRecoEdep.h"
#include "MCBase/MCTrack.h"

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
