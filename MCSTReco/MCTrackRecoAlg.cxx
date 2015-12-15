////////////////////////////////////////////////////////////////////////
//
//  MCTrackRecoAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCTRACKRECOALG_CXX
#define MCTRACKRECOALG_CXX

#include "MCTrackRecoAlg.h"
#include <iostream>
#include <fstream>

namespace sim {

  //##################################################################
  MCTrackRecoAlg::MCTrackRecoAlg(fhicl::ParameterSet const& pset) 
  //##################################################################
  {
    fDebugMode = pset.get<bool>("DebugMode");
  }

  void MCTrackRecoAlg::Reconstruct(MCRecoPart& part_v,
				   MCRecoEdep& edep_v)
  {

    fMCTrack.clear();
    for(size_t i=0; i<part_v.size(); ++i) {

      auto const& mini_part = part_v[i];

      if( part_v._pdg_list.find(mini_part._pdgcode) == part_v._pdg_list.end() ) continue;

      ::sim::MCTrack mini_track;
      
      std::vector<double> dEdx; 
      std::vector<std::vector<double> > dQdx; 

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



      for(auto const& vtx_mom : mini_part._det_path){	
	mini_track.push_back(MCStep(vtx_mom.first,vtx_mom.second));
      }

      auto const& edep_index = edep_v.TrackToEdepIndex(mini_part._track_id);
      auto const& edeps = edep_v.GetEdepArrayAt(edep_index);      


      int n = 0;

      for(auto const& step_trk : mini_track){

        if( int(&step_trk - &mini_track[0])+1 == int(mini_track.size()) ){  //annoying way to check if this is last step
	  std::cout << "\t\t\t::: TRACK END!" << std::endl;
	  break;}
	
	
	auto const& nxt_step_trk = mini_track.at(int(&step_trk - &mini_track[0])+1);      	

	//Defining dEdx	
	
	//Find the distance step-to-step
	double dist = sqrt(pow(step_trk.X() - nxt_step_trk.X(),2) + 
			   pow(step_trk.Y() - nxt_step_trk.Y(),2) +
			   pow(step_trk.Z() - nxt_step_trk.Z(),2));


	//Make a plane at the step pointed at the next step
	
	//Need to define a plane through the first MCStep with a normal along the momentum vector of the step
	//The plane will be defined in the typical way:
	// a*x + b*y + c*z + d = 0
	// where, a = dir_x, b = dir_y, c = dir_z, d = - (a*x_0+b*y_0+c*z_0)
	// then the *signed* distance of any point (x_1, y_1, z_1) from this plane is: 
	// D = (a*x_1 + b*y_1 + c*z_1 + d )/sqrt( pow(a,2) + pow(b,2) + pow(c,2))  	        

	double a = 0, b = 0, c = 0, d = 0;
	a = nxt_step_trk.X() - step_trk.X();
	b = nxt_step_trk.Y() - step_trk.Y();
	c = nxt_step_trk.Z() - step_trk.Z();
	d = -1*(a*step_trk.X() + b*step_trk.Y() + c*step_trk.Z());	
	
	//Make a line connecting the two points 
	// Test point == x_0
	// First Step == x_1
	// Next step == x_2
	// A = x_1 - x_0
	// B = x_2 - x_1

	// B definition 
	TVector3 B(nxt_step_trk.Position().X() - step_trk.Position().X(), 
		   nxt_step_trk.Position().Y() - step_trk.Position().Y(),
		   nxt_step_trk.Position().Z() - step_trk.Position().Z());
	


	double step_dedx = 0; 
	std::vector<double> step_dqdx(3);
	step_dqdx.clear();

      
	for(auto const& edep : edeps){
	  // x_0 definition
	  TVector3 x_0(edep.pos._x, edep.pos._y, edep.pos._z);
	  // A definition 
	  TVector3 A(step_trk.Position().X() - x_0.X(), 
		     step_trk.Position().Y() - x_0.Y(), 
		     step_trk.Position().Z() - x_0.Z());
	  
	  // Distance from the line connecting x_1 and x_2
	  // 
	  //                 [A dot B]^2     [A dot B]^2
	  // d^2 = A^2 - 2* ____________ + ______________
	  //                    B^2             B^2
	  double LineDist = 0;
	 
	  if(B.Mag2() != 0){
	    LineDist = sqrt(A.Mag2() - 2*pow(A*B,2)/B.Mag2() + pow(A*B,2)/B.Mag2());
	  }
	  else{LineDist = 0;}
	  
	  //Planar Distance
	  if( (a*edep.pos._x + b*edep.pos._y + c*edep.pos._z + d)/sqrt( pow(a,2) + pow(b,2) + pow(c,2)) <= dist &&
	      (a*edep.pos._x + b*edep.pos._y + c*edep.pos._z + d)/sqrt( pow(a,2) + pow(b,2) + pow(c,2)) >= 0 && 
	      LineDist < 0.1){
	    
	    int npid = 0;
	    double engy = 0;
	    
	    for(auto const& pid_energy : edep.energy){
	      engy += pid_energy.second;
	      npid++;
	    }
	    engy /= npid;
	    step_dedx += engy;
	  }
	  
	  // Charge
	  auto q_iter = edep.charge.find(edep.pid);
	  if(q_iter != edep.charge.end())
	    step_dqdx[edep.pid.Plane] += (double)((*q_iter).second);	  
	  
	}

	step_dedx /= dist;
	step_dqdx[0] /= dist;
	step_dqdx[1] /= dist;
	step_dqdx[2] /= dist;	

	dEdx.push_back(step_dedx);
	dQdx.push_back(step_dqdx);

      }
      mini_track.dEdx(dEdx);
      mini_track.dQdx(dQdx);
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
