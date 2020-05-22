////////////////////////////////////////////////////////////////////////
/// \file SimListUtils.h
///
/// Utility functions for getting the various list types in the Simulation
/// package from the Event Data Model
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SIMLISTUTILS_H
#define SIMLISTUTILS_H

#include <string>

#include "art/Framework/Principal/fwd.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/Simulation/LArVoxelList.h"

namespace sim {
  class SimListUtils {
  public:
    static sim::LArVoxelList GetLArVoxelList(const art::Event& evt, std::string moduleLabel);
    static sim::SimPhotonsCollection GetSimPhotonsCollection(const art::Event& evt,
                                                             std::string moduleLabel);

  }; // class SimListUtils
} // namespace sim
#endif // SIMLISTUTILS_H
