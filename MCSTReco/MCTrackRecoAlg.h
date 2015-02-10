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
    MCTrackRecoAlg(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCTrackRecoAlg(){};

    void Reconstruct(const MCRecoPart& part_v, const MCRecoEdep& edep_v);

    const std::vector<sim::MCTrack>& MCTrack() const { return fMCTrack; }

  protected:

    bool             fDebugMode;
    
    std::vector<sim::MCTrack> fMCTrack;

  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
