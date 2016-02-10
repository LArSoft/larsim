////////////////////////////////////////////////////////////////////////
//
//  MCRecoEdep source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCRECOEDEP_CC
#define MCRECOEDEP_CC

#include "MCRecoEdep.h"

namespace sim {

  //const unsigned short MCEdepHit::kINVALID_USHORT = std::numeric_limits<unsigned short>::max();

  //const short MCEdep::kINVALID_SHORT = std::numeric_limits<short>::max();

  //##################################################################
  MCRecoEdep::MCRecoEdep(fhicl::ParameterSet const& pset)
  //##################################################################
  {
    _debug_mode = pset.get<bool>("DebugMode");
    _save_mchit = pset.get<bool>("SaveMCHit");
  }

  const std::vector<sim::MCEdep>& MCRecoEdep::GetEdepArrayAt(size_t edep_index) const
  {
    if(edep_index >= _mc_edeps.size()) 
      
      throw cet::exception(__FUNCTION__) << Form("Track ID %zu not found!",edep_index);

    return _mc_edeps.at(edep_index);
  }

  std::vector<sim::MCEdep>& MCRecoEdep::__GetEdepArray__(unsigned int track_id)
  { 
    if(ExistTrack(track_id)) return _mc_edeps.at((*_track_index.find(track_id)).second);

    _track_index.insert(std::pair<unsigned int,size_t>(track_id,_mc_edeps.size()));
    
    _mc_edeps.push_back(std::vector<sim::MCEdep>());
    
    return (*(_mc_edeps.rbegin()));
  }

  void MCRecoEdep::MakeMCEdep(const std::vector<sim::SimChannel>& schArray)
  {
    _mc_edeps.clear();
    _track_index.clear();

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detp;

    // Key map to identify a unique particle energy deposition point
    std::map<unsigned int, std::map<UniquePosition, int> > hit_index_m;

    if(_debug_mode) std::cout<<"Processing "<<schArray.size()<<" channels..."<<std::endl;
    // Loop over channels
    for(size_t i=0; i<schArray.size(); ++i) {

      // Get data to loop over
      auto const& sch = schArray[i];
      const std::map<unsigned short,std::vector<sim::IDE> > sch_map(sch.TDCIDEMap());
      // Channel
      UInt_t ch = sch.Channel();

      // Loop over ticks
      for(auto tdc_iter = sch_map.begin(); tdc_iter!=sch_map.end(); ++tdc_iter) {
	
	unsigned short hit_time = (*tdc_iter).first;

	// Loop over IDEs
	for(auto const ide : (*tdc_iter).second) {

	  int track_id = ide.trackID;
	  if(track_id < 0) track_id = track_id * (-1);
	  unsigned int real_track_id = track_id;

	  UniquePosition pos(ide.x, ide.y, ide.z);

	  int hit_index = -1;
	  
	  auto hit_index_track_iter = hit_index_m.find(real_track_id);
	  //std::cout<<"Inspecting: hit-time="<<hit_time<<" : PartID="<<real_track_id<<std::endl;
	  if(hit_index_track_iter == hit_index_m.end()) {
	    // create new entry here
	    //_track_index.insert(std::pair<unsigned int,size_t>(real_track_id,_mc_edeps.size()));
	    //_mc_edeps.push_back(std::vector<sim::MCEdep>());
	    int new_hit_index = this->__GetEdepArray__(real_track_id).size();
            hit_index_m[real_track_id].insert(std::make_pair(pos,new_hit_index));
	  }
	  else{
	    
	    auto hit_index_pos_iter = (*hit_index_track_iter).second.find(pos);
	    
	    if(hit_index_pos_iter == (*hit_index_track_iter).second.end()) {
	      int new_hit_index = this->__GetEdepArray__(real_track_id).size();
	      (*hit_index_track_iter).second.insert(std::make_pair(pos,new_hit_index));
	      
	    }
	    else
	      hit_index = (*hit_index_pos_iter).second;
	  }
	  //std::cout<<"Finished checking: hit-time="<<hit_time<<" : PartID="<<real_track_id<<std::endl;
	  if(hit_index < 0) {
	    
	    // This particle energy deposition is never recorded so far. Create a new Edep
	    MCEdep edep;
	    // Fill Edep
	    edep.pos = pos;
	    //float charge = ide.numElectrons * detp->ElectronsToADC();
	    double charge = ide.numElectrons;

	    for(auto const& pid : geom->PlaneIDs()){
	      edep.charge[pid] = 0;
	      edep.energy[pid] = 0;
	    }

	    auto pid = geom->ChannelToWire(ch)[0].planeID();
	    //edep.energy      = ide.energy;
	    edep.energy[pid] = ide.energy;
	    edep.charge[pid] = charge;
	    edep.pid         = pid;

	    // If configured to save MC hits, do so
	    if(_save_mchit) {
	      MCEdepHit mchit;
	      mchit.timeStart = hit_time;
	      mchit.timeMax = hit_time;
	      mchit.timeEnd = hit_time;
	      mchit.qSum = charge;
	      mchit.qMax = charge;
	      edep.mchits.insert(std::make_pair(ch,mchit));
	    }
	    //std::cout<<"Inserting: "<<this->__GetEdepArray__(real_track_id).size()<<std::flush;
	    this->__GetEdepArray__(real_track_id).push_back(edep);
	    //std::cout<<" now became "<<this->__GetEdepArray__(real_track_id).size()<<std::endl;
	  }else{
	    //std::cout<<"Existing ... "<<hit_index<<" (size = "<<this->__GetEdepArray__(real_track_id).size()<<std::endl;
	    // Append charge to the relevant edep (@ hit_index)
	    //float charge = ide.numElectrons * detp->ElectronsToADC();
	    double charge = ide.numElectrons;
	    MCEdep &edep = this->__GetEdepArray__(real_track_id).at(hit_index);

	    auto pid = geom->ChannelToWire(ch)[0].planeID();
	    edep.charge[pid] += charge;
	    //if(pid == edep.pid) edep.energy += ide.energy;
	    edep.energy[pid] += ide.energy;

	    // If configured to store hit, store
	    if(_save_mchit) {
	      
	      // Check if channel is already stored. If so, append. Else store new.
	      auto hit_iter = edep.mchits.find(ch);
	      if( hit_iter == edep.mchits.end() ) {
		MCEdepHit mchit;
		mchit.timeStart = hit_time;
		mchit.timeMax = hit_time;
		mchit.timeEnd = hit_time;
		mchit.qSum = charge;
		mchit.qMax = charge;
		edep.mchits.insert(std::pair<unsigned short,sim::MCEdepHit>(ch,mchit));

	      }else{

		if( (*hit_iter).second.timeStart > hit_time )
		  throw cet::exception(__FUNCTION__) << "Found hit before the already-registered hit's start time!";
		(*hit_iter).second.qSum += charge;
		(*hit_iter).second.timeEnd = hit_time;
		if( (*hit_iter).second.qMax < charge ) {
		  (*hit_iter).second.qMax = charge;
		  (*hit_iter).second.timeMax = hit_time;
		}
		
	      }
	    }
	  }
	  
	} // end looping over ides in this tick
	//std::cout<<"End looping over ides in tick: "<<hit_time<<std::endl;
      } // end looping over ticks in this channel
      //std::cout<<"End looping over ticks in channel: "<<ch<<std::endl;    
    }// end looping over channels

    if(_debug_mode) {

      std::cout<< Form("  Collected %zu particles' energy depositions...",_mc_edeps.size()) << std::endl;
      for(auto const& track_id_index : _track_index ) {

	auto track_id   = track_id_index.first;
	auto edep_index = track_id_index.second;

	//	std::cout<< Form("    Track ID: %d ... %zu Edep!", track_id, edep_index) << std::endl;

      }
      std::cout<<std::endl;
    }
    


  }



}

#endif
