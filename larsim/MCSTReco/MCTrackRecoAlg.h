#ifndef MCTRACKRECOALG_H
#define MCTRACKRECOALG_H

// ART includes
namespace fhicl {
  class ParameterSet;
}

// LArSoft
#include "lardataobj/MCBase/MCTrack.h"

namespace sim {
  class MCRecoEdep;
  class MCRecoPart;
}

// STL
#include <memory>
#include <vector>

namespace sim {

  class MCTrackRecoAlg {

  public:
    /// Default constructor with fhicl parameters
    explicit MCTrackRecoAlg(fhicl::ParameterSet const& pset);
    std::unique_ptr<std::vector<sim::MCTrack>> Reconstruct(MCRecoPart& part_v, MCRecoEdep& edep_v);

  protected:
    bool fDebugMode;

  }; // class MCShowerHitRecoAlg

} //namespace cluster
#endif
