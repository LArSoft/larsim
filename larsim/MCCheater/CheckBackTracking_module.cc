////////////////////////////////////////////////////////////////////////
//
// CheckBackTracking module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include <string>

// LArSoft includes
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nug4/ParticleNavigation/ParticleList.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace cheat {
  class CheckBackTracking : public art::EDAnalyzer {
  public:
    explicit CheckBackTracking(fhicl::ParameterSet const& pset);

    void analyze(art::Event const& evt);

 private:

    std::string fHitModuleLabel;    ///< label for module creating recob::Hit objects
    std::string fG4ModuleLabel;     ///< label for module running G4 and making particles, etc

  };
}

namespace cheat{

  //--------------------------------------------------------------------
  CheckBackTracking::CheckBackTracking(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fHitModuleLabel{pset.get< std::string >("HitModuleLabel",    "ffthit"  )}
    , fG4ModuleLabel {pset.get< std::string >("G4ModuleLabel",     "largeant")}
  {}

  //--------------------------------------------------------------------
  void CheckBackTracking::analyze(art::Event const& evt)
  {

    // grab the hits that have been reconstructed
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fHitModuleLabel, hitcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitcol);

    // loop over the hits and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();

    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;

    // make a collection of the distinct eve ID values
    std::set<int> eveIDs;

    while( itr != hits.end() ){

      // print the truth information for this hit
      mf::LogInfo("CheckBackTracking") << *((*itr).get()) << "\n channel is: " << (*itr)->Channel();

      std::vector<sim::TrackIDE> trackides = bt_serv->HitToTrackIDEs(*itr);
      std::vector<sim::TrackIDE> eveides   = bt_serv->HitToEveTrackIDEs(*itr);
      std::vector<double>          xyz       = bt_serv->HitToXYZ(*itr);

      mf::LogInfo("CheckBackTracking") << "hit weighted mean position is ("
				       << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")";

      for(size_t t = 0; t < trackides.size(); ++t){

	// find the Eve particle for the current trackID
	int eveID = pi_serv->ParticleList().EveId( trackides[t].trackID );

	mf::LogInfo("CheckBackTracking") << "track id: " << trackides[t].trackID
					 << " contributed " << trackides[t].energy << "/"
					 << trackides[t].energyFrac
					 << " to the current hit and has eveID: "
					 << eveID;
      }

      for(size_t e = 0; e < eveides.size(); ++e){
	mf::LogInfo("CheckBackTracking") << "eve id: " << eveides[e].trackID
					 << " contributed " << eveides[e].energy << "/"
					 << eveides[e].energyFrac
					 << " to the current hit";

	if(eveIDs.find(eveides[e].trackID) == eveIDs.end()) eveIDs.insert(eveides[e].trackID);
      }

      itr++;
    }// end loop over hits

    // loop over the eveID values and calculate the purity and efficiency for each
    std::set<int>::iterator setitr = eveIDs.begin();
    while( setitr != eveIDs.end() ){

      std::set<int> id;
      id.insert(*setitr);
      mf::LogInfo("CheckBackTracking") << "eve ID: " << *setitr
				       << " purity: "
				       << bt_serv->HitCollectionPurity(id, hits)
				       << " efficiency: "
				       << bt_serv->HitCollectionEfficiency(id, hits, hits, geo::k3D);


      setitr++;
    }

    return;

  } // end analyze

} // end namespace

namespace cheat{

  DEFINE_ART_MODULE(CheckBackTracking)

}
