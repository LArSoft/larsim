////////////////////////////////////////////////////////////////////////
/// \file SimListUtils.h
///
/// Utility functions for getting the various list types in the Simulation
/// package from the Event Data Model
///
/// \version $Id: SingleGen.h,v 1.2 2010/02/15 19:10:40 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SIMLISTUTILS_H
#define SIMLISTUTILS_H

#include <vector>
#include <string>
#include "larsim/Simulation/LArVoxelList.h"
#include "larsim/Simulation/SimPhotons.h"

namespace art { class Event; }


namespace sim{
  class SimListUtils {
  public:
    SimListUtils();
    virtual ~SimListUtils();
    
    static sim::LArVoxelList         GetLArVoxelList        (const art::Event& evt, std::string moduleLabel); 
    static sim::SimPhotonsCollection GetSimPhotonsCollection(const art::Event& evt, std::string moduleLabel);
    
  }; // class SimListUtils
} //namespace sim
#endif // SIMLISTUTILS_H
