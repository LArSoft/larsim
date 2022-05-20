////////////////////////////////////////////////////////////////////////////
//
// ParticleInventoryService.cc
// Author: JStock
// EMail:  jason.stock@mines.sdsmt.edu
// 2017-09-12
//
// Maintinence Notes: When the ParticleInventory is initialized, none of the prep work (previously
// the BackTracker rebuild stage) will be done. Each function needs to check and make sure the
// needed data products have been loaded. To see what objects a function uses, you will have to
// check the appropriate part of ParticleInventory. After this, you will need to manually write the check
// into whatever function you are writing. You will also want to include a call to prepare the needed items
// if your check fails.
//
// Example:
// std::set<int> ParticleInventoryService::GetSetOfTrackIds(){
//   if(!this->priv_ParticleListReady()){this->priv_PrepParticleList();} //The GetTrackIds in ParticleInventory needs the ParticleList.
//                                                                         So, we check if it's ready, and if it isn't we ready it.
//   return ParticleInventory::GetSetOfTrackIds();
// }
//
// If you have any questions about how to incorperate something in here, let me know. I know this is a rather odd
// use model. The rationale is to allow the BackTracker service to be lazy, while at the same time allowing gallery
// to use backtracker functions (the gallery implimentation is not lazy).
////////////////////////////////////////////////////////////////////////////

#include "larsim/MCCheater/ParticleInventoryService.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

DEFINE_ART_SERVICE(cheat::ParticleInventoryService)
