////////////////////////////////////////////////////////////////////////////////////////
//
// \file BackTrackerService_service.cc
// \brief A service for backtracking reconstruction information to its truth
// information
//
// \author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
//
////////////////////////////////////////////////////////////////////////////////////////

#include "larsim/MCCheater/BackTrackerService.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

DEFINE_ART_SERVICE(cheat::BackTrackerService)
