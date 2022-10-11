////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parameters_service.cc
/// \brief Store parameters for running LArG4
///
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This service exists to pass parameters to various different
// classes in LArG4, which are not necessary directly called by
// the LArG4_module class.
//
// Ben Jones, MIT, March 2010

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "larsim/Simulation/LArG4Parameters.h"

DEFINE_ART_SERVICE(sim::LArG4Parameters)
