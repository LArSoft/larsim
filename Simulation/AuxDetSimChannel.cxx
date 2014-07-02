////////////////////////////////////////////////////////////////////////
///
/// \file  Simulation/AuxDetSimChannel.cxx
///
/// \author  miceli@fnal.gov
///
////////////////////////////////////////////////////////////////////////


#include "Simulation/AuxDetSimChannel.h"
#include "SimpleTypesAndConstants/PhysicalConstants.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace sim{

  // Default constructor
  //-------------------------------------------------
  AuxDetIDE::AuxDetIDE()
  : trackID        (util::kBogusI)
  , energyDeposited(util::kBogusF)
  , entryX         (util::kBogusF)
  , entryY         (util::kBogusF)
  , entryZ         (util::kBogusF)
  , entryT         (util::kBogusF)
  , exitX          (util::kBogusF)
  , exitY          (util::kBogusF)
  , exitZ          (util::kBogusF)
  , exitT          (util::kBogusF)
  , exitMomentumX  (util::kBogusF)
  , exitMomentumY  (util::kBogusF)
  , exitMomentumZ  (util::kBogusF)
  {}

	AuxDetSimChannel::AuxDetSimChannel()
  : fAuxDetID(UINT_MAX)
  {
  }
  
  AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID, std::set<sim::AuxDetIDE> inputAuxDetIDEs)
  {
    fAuxDetID = inputAuxDetID;
    fAuxDetIDEs = inputAuxDetIDEs;
  }

}//namespace sim
