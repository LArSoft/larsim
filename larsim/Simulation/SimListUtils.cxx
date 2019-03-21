////////////////////////////////////////////////////////////////////////
/// \file SimListUtils.cxx
//
/// \author  brebel@fnal.gov
/// 
/// this class is designed to hold methods that access the event handle
/// to make the various simulation lists, ie ParticleList, LArVoxelList, etc
////////////////////////////////////////////////////////////////////////

#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace sim{

  //----------------------------------------------------------------------
  SimListUtils::SimListUtils()
  {
  }
  
  //----------------------------------------------------------------------
  SimListUtils::~SimListUtils()
  {
  }
  
  //----------------------------------------------------------------------
  // moduleLabel is the label of the module that created the voxels you 
  // are putting into the list
  sim::LArVoxelList SimListUtils::GetLArVoxelList(const art::Event& evt, std::string moduleLabel)
  {
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // get the sim::SimChannels
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(moduleLabel, sccol);

    sim::LArVoxelList voxList;

    // loop over the voxels and put them into the list
    for(auto itr = sccol.begin(); itr != sccol.end(); ++itr){
      
      // get all sim::IDE associated with this channel
      const auto &idemap = (*itr)->TDCIDEMap();
      //std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mitr;

      // loop over all the sim::IDE values
      for(auto mitr = idemap.begin(); mitr != idemap.end(); mitr++){

	double time = (*mitr).first - detprop->TriggerOffset();
	time *= detprop->SamplingRate();
	
	// loop over the sim::IDE objects
	const std::vector<sim::IDE> &ide = (*mitr).second;
	for(size_t i = 0; i < ide.size(); ++i){

	  sim::LArVoxelID larVoxelID(ide[i].x, 
				     ide[i].y,
				     ide[i].z,
				     time);

	  // if energy is unassigned the TrackId is sim::kNoParticleId
	  voxList.Add(larVoxelID, ide[i].numElectrons/lgp->GeVToElectrons(), ide[i].trackID);

	  // set the voxel id for the just added LArVoxelData
	  (*voxList.find(larVoxelID)).second.SetVoxelID(larVoxelID);

	}// end loop over ide for this tdc
      }// end loop over map
    }// end loop over sim::SimChannels

    return voxList;
  }

  //----------------------------------------------------------------------
  // moduleLabel is the label of the module that created the pmthits you 
  // are putting into the list
  sim::SimPhotonsCollection SimListUtils::GetSimPhotonsCollection(const art::Event& evt, std::string moduleLabel)
  {
    /// get the voxels from the event handle
    art::Handle< std::vector<sim::SimPhotons> > pmtHandle;
    evt.getByLabel(moduleLabel, pmtHandle);
    const std::vector<sim::SimPhotons>& pmt(*pmtHandle);

    sim::SimPhotonsCollection pmtList;
    pmtList.clear();
    //std::cout << "Building SimPhotonsCollection" << std::endl;

    /// loop over the pmthits and put them into the list
    for(auto itr = pmt.begin(); itr != pmt.end(); ++itr){

      int ch = (*itr).OpChannel();
      /// make an entry in the list for this pmt id
      if(pmtList.find(ch) == pmtList.end()) {
	// Create new photon object
	sim::SimPhotons new_photons;
	new_photons.clear();
	new_photons.SetChannel(ch);
	new_photons.reserve((*itr).size());
	pmtList.insert(std::pair<int,sim::SimPhotons>(ch,new_photons));
      }
      
      /// add the photons to the entry
      for(auto pitr = (*itr).begin(); pitr != (*itr).end(); ++pitr)
	pmtList[ch].push_back(sim::OnePhoton((*pitr)));
    }
    
    return pmtList;
    //return std::move(pmtList);
  }


}//end namespace util
