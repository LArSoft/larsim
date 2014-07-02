////////////////////////////////////////////////////////////////////////
/// \file  AuxDetSimChannel.h
///
/// \brief object containing MC truth information necessary for making RawDigits 
/// and doing back tracking
///
/// \author  miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef SIM_AUXDETSIMCHANNEL_H
#define SIM_AUXDETSIMCHANNEL_H

#include <vector>
#include <set>
#include <stdint.h>

#include "SimpleTypesAndConstants/geo_types.h"

namespace sim {

  class AuxDetIDE{
    
  public:
    AuxDetIDE();
    
    int   trackID;                       ///< Geant4 supplied track ID
    mutable float energyDeposited;       ///< total energy deposited for this track ID and time
    float entryX;                        ///< Entry position X of particle
    float entryY;                        ///< Entry position Y of particle
    float entryZ;                        ///< Entry position Z of particle
    float entryT;                        ///< Entry time of particle
    mutable float exitX;                 ///< Exit position X of particle
    mutable float exitY;                 ///< Exit position Y of particle
    mutable float exitZ;                 ///< Exit position Z of particle
    mutable float exitT;                 ///< Exit time of particle
    mutable float exitMomentumX;         ///< Exit X-Momentum of particle
    mutable float exitMomentumY;         ///< Exit Y-Momentum of particle
    mutable float exitMomentumZ;         ///< Exit Z-Momentum of particle

#ifndef __GCCXML__
    bool operator<  (const AuxDetIDE& other) const;
    bool operator== (const AuxDetIDE& other) const;

#endif
};

  class AuxDetSimChannel{
    
  public:
    AuxDetSimChannel();
    AuxDetSimChannel(uint32_t inputAuxDetID, std::set<sim::AuxDetIDE> inputAuxDetIDEs);
    
  private:
    uint32_t                 fAuxDetID;   ///< geo->AuxDet(auxDetID), integer used to retrieve AuxDetGeo object
    std::set<sim::AuxDetIDE> fAuxDetIDEs; ///< one fAuxDetIDE for each G4 track id

#ifndef __GCCXML__
  public:
    
    //getters
    uint32_t AuxDetID() const;

    std::set<sim::AuxDetIDE> const& AuxDetIDEs() const;

    //setters
//    void SetAuxDetIDEs(std::set<sim::AuxDetIDE> inputAuxDetIDEs) {fAuxDetIDEs = inputAuxDetIDEs;};

#endif
		
		
  };

} // namespace sim

#ifndef __GCCXML__

inline bool sim::AuxDetIDE::operator<  (const AuxDetIDE& other) const { return trackID < other.trackID;  }
inline bool sim::AuxDetIDE::operator== (const AuxDetIDE& other) const { return other.trackID == trackID; }
inline uint32_t  sim::AuxDetSimChannel::AuxDetID()              const { return fAuxDetID;                }
inline std::set<sim::AuxDetIDE> const& sim::AuxDetSimChannel::AuxDetIDEs() const { return fAuxDetIDEs; }
#endif

#endif // SIM_AUXDETSIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
