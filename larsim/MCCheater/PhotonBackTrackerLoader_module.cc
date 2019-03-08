////////////////////////////////////////////////////////////////////////
// Class:       PhotonBackTrackerLoader
// Module Type: producer
// File:        PhotonBackTrackerLoader.h
//
// Generated at Thu Jun 14 06:49:31 2012 by Brian Rebel using artmod
// from art v1_00_11.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cheat {
  class PhotonBackTrackerLoader;
}

class cheat::PhotonBackTrackerLoader : public art::EDProducer {
public:
  explicit PhotonBackTrackerLoader(fhicl::ParameterSet  const& p);
  virtual ~PhotonBackTrackerLoader();

  virtual void produce(art::Event & e);


private:

  // Declare member data here.

};

//------------------------------------------------------------------------------
cheat::PhotonBackTrackerLoader::PhotonBackTrackerLoader(fhicl::ParameterSet  const& p)
  : EDProducer{p}
{
  // Call appropriate Produces<>() functions here.
}

//------------------------------------------------------------------------------
cheat::PhotonBackTrackerLoader::~PhotonBackTrackerLoader()
{
  // Clean up dynamic memory and other resources here.
}

//------------------------------------------------------------------------------
// the sole purpose of this module is to issue the Rebuild command to the
// PhotonBackTracker service.  It should be put after all simulation data producing
// modules have run in the job, and only in jobs that create the simulation and
// then make use of the PhotonBackTracker in either cheating reconstruction modules
// or analyzers
void cheat::PhotonBackTrackerLoader::produce(art::Event & e)
{
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;
  pi_serv->Rebuild(e);
  pbt_serv->Rebuild(e);

  return;
}


DEFINE_ART_MODULE(cheat::PhotonBackTrackerLoader)
