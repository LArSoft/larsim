#ifndef MCSHOWERRECOALG_H
#define MCSHOWERRECOALG_H

// ART includes
namespace fhicl { class ParameterSet; }

// LArSoft
#include "MCShowerRecoPart.h"
#include "lardataobj/MCBase/MCShower.h"
namespace sim {
  class MCRecoEdep;
  class MCRecoPart;
  class MCShower;
}

// STL
#include <vector>
#include <memory>

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
