////////////////////////////////////////////////////////////////////////
///
/// \file  Simulation/AuxDetSimChannel.cxx
///
/// \author  miceli@fnal.gov
///
////////////////////////////////////////////////////////////////////////


// our header
#include "Simulation/AuxDetSimChannel.h"

// C/C++ standard library
#include <limits> // std::numeric_limits<>

// LArSoft headers
#include "SimpleTypesAndConstants/PhysicalConstants.h" // util::kBogusX

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
    : fAuxDetID(std::numeric_limits<uint32_t>::max())
  {
  }
  
  AuxDetSimChannel::AuxDetSimChannel
    (uint32_t inputAuxDetID, const std::vector<sim::AuxDetIDE>& inputAuxDetIDEs):
    fAuxDetID(inputAuxDetID), fAuxDetIDEs(inputAuxDetIDEs)
  {}

  AuxDetSimChannel::AuxDetSimChannel
    (uint32_t inputAuxDetID, std::vector<sim::AuxDetIDE>&& inputAuxDetIDEs):
    fAuxDetID(inputAuxDetID), fAuxDetIDEs(inputAuxDetIDEs)
  {}

}//namespace sim
