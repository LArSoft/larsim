////////////////////////////////////////////////////////////////////////
///
/// \file  Simulation/AuxDetSimChannel.cxx
///
/// \author  miceli@fnal.gov
///
////////////////////////////////////////////////////////////////////////


// our header
#include "larsim/Simulation/AuxDetSimChannel.h"

// C/C++ standard library
#include <limits> // std::numeric_limits<>
#include <stdexcept>

// LArSoft headers
#include "larcore/SimpleTypesAndConstants/PhysicalConstants.h" // util::kBogusX

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

  // Copy with offset constructor
  //-------------------------------------------------
  AuxDetIDE::AuxDetIDE(AuxDetIDE const& ide, int offset)
    : trackID        (ide.trackID+offset)
    , energyDeposited(ide.energyDeposited)
    , entryX         (ide.entryX)
    , entryY         (ide.entryY)
    , entryZ         (ide.entryZ)
    , entryT         (ide.entryT)
    , exitX          (ide.exitX)
    , exitY          (ide.exitY)
    , exitZ          (ide.exitZ)
    , exitT          (ide.exitT)
    , exitMomentumX  (ide.exitMomentumX)
    , exitMomentumY  (ide.exitMomentumY)
    , exitMomentumZ  (ide.exitMomentumZ)
    {}

  //----------------------------------------------------------------------------
  AuxDetSimChannel::AuxDetSimChannel()
    : fAuxDetID(std::numeric_limits<uint32_t>::max())
    , fAuxDetSensitiveID(std::numeric_limits<uint32_t>::max())
  {
  }

  //----------------------------------------------------------------------------
  AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID,
				     uint32_t inputAuxDetSensitiveID)
      : fAuxDetID(inputAuxDetID)
      , fAuxDetSensitiveID(inputAuxDetSensitiveID)
    {}
  
  //----------------------------------------------------------------------------
  AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID, 
				     const std::vector<sim::AuxDetIDE>& inputAuxDetIDEs,
				     uint32_t inputAuxDetSensitiveID)
      : fAuxDetID(inputAuxDetID)
      , fAuxDetSensitiveID(inputAuxDetSensitiveID)
      , fAuxDetIDEs(inputAuxDetIDEs)
  {}

  //----------------------------------------------------------------------------
  AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID, 
				     std::vector<sim::AuxDetIDE>&& inputAuxDetIDEs,
				     uint32_t inputAuxDetSensitiveID)
    : fAuxDetID(inputAuxDetID)
    , fAuxDetSensitiveID(inputAuxDetSensitiveID)
    , fAuxDetIDEs(inputAuxDetIDEs)
  {}

  //----------------------------------------------------------------------------
  std::pair<int,int> AuxDetSimChannel::MergeAuxDetSimChannel(const AuxDetSimChannel& chan, 
							     int offset)
  {
    if(this->fAuxDetID != chan.AuxDetID() && this->fAuxDetSensitiveID != chan.AuxDetSensitiveID())
      throw std::runtime_error("ERROR AuxDetSimChannel Merge: Trying to merge different channels!");

    std::pair<int,int> range_trackID(std::numeric_limits<int>::max(),
				     std::numeric_limits<int>::min());

    for(auto const& ide : AuxDetIDEs()){
      this->fAuxDetIDEs.emplace_back(ide,offset);

      if( ide.trackID+offset < range_trackID.first  )
	range_trackID.first = ide.trackID+offset;
      if( ide.trackID+offset > range_trackID.second )
	range_trackID.second = ide.trackID+offset;
    }

    return range_trackID;
  }

  //----------------------------------------------------------------------------
  bool AuxDetSimChannel::operator<  (const sim::AuxDetSimChannel& other) const 
  { 
    if(fAuxDetID < other.AuxDetID() ) return true;

    return fAuxDetSensitiveID < other.AuxDetSensitiveID(); 
  }

  //----------------------------------------------------------------------------
  bool AuxDetSimChannel::operator== (const sim::AuxDetSimChannel& other) const 
  { 
    return (fAuxDetID == other.AuxDetID() && fAuxDetSensitiveID == other.AuxDetSensitiveID()); 
  }

  
}//namespace sim
