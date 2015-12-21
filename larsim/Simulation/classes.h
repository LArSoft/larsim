//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by klg
//
// Notes:
// 1) The system is not able to deal with
//    art::Wrapper<std::vector<std::string> >;
//    The problem is somewhere inside root's reflex mechanism
//    and Philippe Canal says that it is ( as of March 2010) a
//    known problem.  He also says that they do not have any
//    plans to fix it soon.  We can always work around it 
//    by putting the string inside another object.

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"

// nutools includes
#include "SimulationBase/MCTruth.h"

// Simulation includes
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/SimPhotons.h"
#include "larsim/Simulation/BeamGateInfo.h"
#include "larsim/Simulation/AuxDetSimChannel.h"

namespace {

  std::vector<sim::IDE> dummy1;
  std::pair<unsigned short, std::vector<sim::IDE> > dummy2;
}

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class art::Wrapper< std::vector<int>                  >;
template class art::Wrapper< std::vector<sim::SimChannel>      >;
template class art::Wrapper< std::vector<sim::SimPhotons>      >;
template class art::Wrapper< std::vector<sim::OnePhoton>       >;
template class art::Wrapper< std::vector<sim::SimPhotonsLite>  >;
template class art::Wrapper< std::vector<sim::BeamGateInfo>    >;
template class art::Wrapper< std::vector<sim::AuxDetSimChannel> >;

/// \todo: Remove this line after DUNE makes new files that do not contain
///        DUNE10ktPhotons.  That data product is now called SimPhotonsLite
///        11 September 2013 brebel
template class art::Wrapper< std::vector<sim::DUNE10ktPhotons> >;

