////////////////////////////////////////////////////////////////////////
//
//  MCTrackRecoAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCTRACKRECOALG_CXX
#define MCTRACKRECOALG_CXX

#include "MCTrackRecoAlg.h"

namespace sim {

  //##################################################################
  MCTrackRecoAlg::MCTrackRecoAlg(fhicl::ParameterSet const& pset) 
  //##################################################################
  {
    fDebugMode = pset.get<bool>("DebugMode");
  }

  void MCTrackRecoAlg::Reconstruct(const MCRecoPart& part_v,
				   const MCRecoEdep& edep_v)
  {

    fMCTrack.clear();
    for(size_t i=0; i<part_v.size(); ++i) {

      auto const& mini_part = part_v[i];

      if( part_v._pdg_list.find(mini_part._pdgcode) == part_v._pdg_list.end() ) continue;

      ::sim::MCTrack mini_track;

      mini_track.Origin  ( mini_part._origin   );
      mini_track.PdgCode ( mini_part._pdgcode  );
      mini_track.TrackID ( mini_part._track_id );
      mini_track.Process ( mini_part._process  );
      mini_track.Start   ( MCStep( mini_part._start_vtx, mini_part._start_mom ) );
      mini_track.End     ( MCStep( mini_part._end_vtx,   mini_part._end_mom   ) );

      unsigned int mother_track   = part_v.MotherTrackID(i);
      unsigned int ancestor_track = part_v.AncestorTrackID(i);

      if(mother_track == kINVALID_UINT || ancestor_track == kINVALID_UINT)

	throw cet::exception(__FUNCTION__) << "LOGIC ERROR: mother/ancestor track ID is invalid!";

      MCMiniPart mother_part;
      MCMiniPart ancestor_part;

      unsigned int mother_index   = part_v.TrackToParticleIndex(mother_track);
      unsigned int ancestor_index = part_v.TrackToParticleIndex(ancestor_track);

      if(mother_index != kINVALID_UINT)   mother_part   = part_v[mother_index];
      else mother_part._track_id = mother_track;

      if(ancestor_index != kINVALID_UINT) ancestor_part = part_v[ancestor_index];
      else ancestor_part._track_id = ancestor_track;

      mini_track.MotherPdgCode ( mother_part._pdgcode  );
      mini_track.MotherTrackID ( mother_part._track_id );
      mini_track.MotherProcess ( mother_part._process  );
      mini_track.MotherStart   ( MCStep( mother_part._start_vtx, mother_part._start_mom ) );
      mini_track.MotherEnd     ( MCStep( mother_part._end_vtx,   mother_part._end_mom   ) );

      mini_track.AncestorPdgCode ( ancestor_part._pdgcode  );
      mini_track.AncestorTrackID ( ancestor_part._track_id );
      mini_track.AncestorProcess ( ancestor_part._process  );
      mini_track.AncestorStart   ( MCStep( ancestor_part._start_vtx, ancestor_part._start_mom ) );
      mini_track.AncestorEnd     ( MCStep( ancestor_part._end_vtx,   ancestor_part._end_mom   ) );


      // Fill trajectory points

      for(auto const& vtx_mom : mini_part._det_path)
	
	mini_track.push_back(MCStep(vtx_mom.first,vtx_mom.second));
      /*

      auto const& edep_index = edep_v.TrackToEdepIndex(mini_part._track_id);
      auto const& g4_step_v  = mini_part._det_path;
      if(edep_index>=0 && g4_step_v.size()) {

	// Get energy deposition vector
	auto const& edeps = edep_v.GetEdepArrayAt(edep_index);

	auto const& init_vtx = g4_step_v[0].first;
	auto const& init_mom = g4_step_v[0].second;

	// Initialize "previous step" holder 
	std::vector<double> previous_vtx(4,0);
	std::vector<double> previous_mom(4,0);
	for(size_t i=0; i<4; ++i){
	  previous_vtx[i] = init_vtx[i];
	  previous_mom[i] = init_mom[i];
	}
	double mass = sqrt(pow(init_mom[3],2) -
			   pow(init_mom[2],2) -
			   pow(init_mom[1],2) -
			   pow(init_mom[0],2));

	std::deque<size_t> edep_list;
	for(size_t i=0; i<edeps.size(); ++i) {
	  auto const& edep = edeps[i];
	  if(edep.energy>0 && part_v.InDetector(edep.x,edep.y,edep.z))
	    edep_list.push_back(i);
	}
	std::deque<size_t>::iterator iter, max_iter;

	double min_dist2 = ::sim::kINVALID_DOUBLE;

	TLorentzVector vtx,mom;
	mini_track.reserve(edep_list.size());
	while(edep_list.size()) {

	  min_dist2 = ::sim::kINVALID_DOUBLE;
	  for(iter = edep_list.begin(); iter != edep_list.end(); ++iter) {
	    
	    auto const& edep = edeps[(*iter)];
	    double dist2 = ( pow(edep.x - previous_vtx[0],2) +
			     pow(edep.y - previous_vtx[1],2) +
			     pow(edep.z - previous_vtx[2],2) );
	    if(dist2 < min_dist2) {
	      min_dist2 = dist2;
	      max_iter = iter;
	    }
	  }
	  
	  if(min_dist2 != ::sim::kINVALID_DOUBLE) {

	    min_dist2 = sqrt(min_dist2);
	    auto const& edep = edeps[(*max_iter)];
	    
	    vtx[0] = edep.x;
	    vtx[1] = edep.y;
	    vtx[2] = edep.z;
	    vtx[3] = (min_dist2/100. / 2.998e8)*1.e9 + previous_vtx[3];
	    
	    mom[3] = previous_mom[3] - edep.energy;
	    if(mom[3]<mass) throw cet::exception(__FUNCTION__) << "Negative energy!"<<previous_vtx[3]<<" - "<<edep.energy<< " < "<<mass;
	    double momentum_magnitude = sqrt(pow(mom[3],2)-pow(mass,2));
	    
	    mom[0] = vtx[0] - previous_vtx[0];
	    mom[1] = vtx[1] - previous_vtx[1];
	    mom[2] = vtx[2] - previous_vtx[2];
	    
	    double dir_magnitude = sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2));
	    if(dir_magnitude==0.) throw cet::exception(__FUNCTION__) << "Zero magnitude!" << edep.energy << " : " <<dir_magnitude;
	    mom[0] *= momentum_magnitude/dir_magnitude;
	    mom[1] *= momentum_magnitude/dir_magnitude;
	    mom[2] *= momentum_magnitude/dir_magnitude;
	    
	    mini_track.push_back(MCStep(vtx,mom));
	    
	    for(size_t i=0; i<4; ++i) {
	      previous_vtx[i] = vtx[i];
	      previous_mom[i] = mom[i];
	    }

	    edep_list.erase(max_iter);
	    
	  }else 
	    
	    break;
	  
	}
      }
      */
      fMCTrack.push_back(mini_track);
    }
    
    if(fDebugMode) {

      for(auto const& prof : fMCTrack) {
	
	std::cout
	  
	  << Form("  Track particle:      PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.PdgCode(),prof.TrackID(),
		  prof.Start().X(),prof.Start().Y(),prof.Start().Z(),prof.Start().T(),
		  prof.Start().Px(),prof.Start().Py(),prof.Start().Pz(),prof.Start().E())
	  << std::endl
	  << Form("    Mother particle:   PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.MotherPdgCode(),prof.MotherTrackID(),
		  prof.MotherStart().X(),prof.MotherStart().Y(),prof.MotherStart().Z(),prof.MotherStart().T(),
		  prof.MotherStart().Px(),prof.MotherStart().Py(),prof.MotherStart().Pz(),prof.MotherStart().E())
	  << std::endl
	  << Form("    Ancestor particle: PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.AncestorPdgCode(),prof.AncestorTrackID(),
		  prof.AncestorStart().X(),prof.AncestorStart().Y(),prof.AncestorStart().Z(),prof.AncestorStart().T(),
		  prof.AncestorStart().Px(),prof.AncestorStart().Py(),prof.AncestorStart().Pz(),prof.AncestorStart().E())
	  << std::endl
	  << Form("    ... with %zu trajectory points!",prof.size())
	  << std::endl;

	if(prof.size()) {
	  std::cout 
	    << Form("        Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		    prof[0].X(), prof[0].Y(), prof[0].Z(), prof[0].T(),
		    prof[0].Px(), prof[0].Py(), prof[0].Pz(), prof[0].E())
	    << std::endl
	    << Form("        End @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		    (*prof.rbegin()).X(), (*prof.rbegin()).Y(), (*prof.rbegin()).Z(), (*prof.rbegin()).T(),
		    (*prof.rbegin()).Px(), (*prof.rbegin()).Py(), (*prof.rbegin()).Pz(), (*prof.rbegin()).E())
	    << std::endl;
	}
      }

      std::cout<<std::endl<<std::endl;
    }
  }
}


#endif
