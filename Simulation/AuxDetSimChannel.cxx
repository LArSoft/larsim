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
  : fAuxDetID(0)
  {}
  
  AuxDetSimChannel::AuxDetSimChannel(uint32_t fAuxDetID)
  : fAuxDetID(fAuxDetID)
  {}
  
  void AuxDetSimChannel::AddParticleStep(
                       int   inputTrackID,
                       float inputEnergyDeposited,
                       float inputEntryX,
                       float inputEntryY,
                       float inputEntryZ,
                       float /*inputEntryT*/,
                       float inputExitX,
                       float inputExitY,
                       float inputExitZ,
                       float inputExitT,
                       float inputExitMomentumX,
                       float inputExitMomentumY,
                       float inputExitMomentumZ){
    
    AuxDetIDE auxDetIDE;
    auxDetIDE.trackID         = inputTrackID;
    auxDetIDE.energyDeposited = inputEnergyDeposited;
    auxDetIDE.entryX          = inputEntryX;
    auxDetIDE.entryY          = inputEntryY;
    auxDetIDE.entryZ          = inputEntryZ;
    auxDetIDE.entryT          = inputEntryY;
    auxDetIDE.exitX           = inputExitX;
    auxDetIDE.exitY           = inputExitY;
    auxDetIDE.exitZ           = inputExitZ;
    auxDetIDE.exitT           = inputExitT;
    auxDetIDE.exitMomentumX   = inputExitMomentumX;
    auxDetIDE.exitMomentumY   = inputExitMomentumY;
    auxDetIDE.exitMomentumZ   = inputExitMomentumZ;
    
    std::set<AuxDetIDE>::iterator setItr = fAuxDetIDEs.find(auxDetIDE);
    
    if(setItr != fAuxDetIDEs.end()){ //If trackID is already in the map, update it
      
      (*setItr).energyDeposited += inputEnergyDeposited;
      (*setItr).exitX            = inputExitX;
      (*setItr).exitY            = inputExitY;
      (*setItr).exitZ            = inputExitZ;
      (*setItr).exitT            = inputExitT;
      (*setItr).exitMomentumX    = inputExitMomentumX;
      (*setItr).exitMomentumY    = inputExitMomentumY;
      (*setItr).exitMomentumZ    = inputExitMomentumZ;
      
    }
    else{ //if trackID is not in the set yet, add it
      
      auto insertResult = fAuxDetIDEs.insert(auxDetIDE);
      if(!insertResult.second)
        throw cet::exception("BadAuxDetIDEInsert") << "Track ID: "
                                                   << inputTrackID
                                                   << " already in fAuxDetIDEs set";
    }//else
  
  
  }//AddParticleStep
  
}//namespace sim
