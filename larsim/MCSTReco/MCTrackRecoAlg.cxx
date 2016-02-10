////////////////////////////////////////////////////////////////////////
//
//  MCTrackRecoAlg source
//
//  dEdx and dQdx Estimates Added by Joseph Zennamo (jaz8600@fnal.gov)
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
      dQdx.resize(3);

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
      
      // No calorimetry for zero length tracks...
      // JZ : I think we should remove zero length MCTracks because I do not see their utility
      // JZ : Someone could make this a fcl parameter, I did not
      if(mini_track.size() == 0){
	fMCTrack.push_back(mini_track);
	continue;
      }
      
      auto const& edep_index = edep_v.TrackToEdepIndex(mini_part._track_id);
      if(edep_index < 0 ) continue;
      auto const& edeps = edep_v.GetEdepArrayAt(edep_index);      

      int n = 0;
      
      for(auto const& step_trk : mini_track){
	
        if( int(&step_trk - &mini_track[0])+1 == int(mini_track.size()) ){  //annoying way to check if this is last step
	  continue;}
	
	
	auto const& nxt_step_trk = mini_track.at(int(&step_trk - &mini_track[0])+1);      	

	//Defining the track step-by-step dEdx and dQdx
	
	//Find the distance between two adjacent MCSteps
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
	
	//Make a line connecting the two points and find the distance from that line 
	// 
	//                        [A dot B]^2     [A dot B]^2
	// distance^2 = A^2 - 2* ____________ + ______________
	//                            B^2             B^2	
	// Test point == x_0
	// First Step == x_1
	// Next step == x_2
	// A = x_1 - x_0
	// B = x_2 - x_1

	// 'B' definition 
	TVector3 B(nxt_step_trk.Position().X() - step_trk.Position().X(), 
		   nxt_step_trk.Position().Y() - step_trk.Position().Y(),
		   nxt_step_trk.Position().Z() - step_trk.Position().Z());
	

	//Initialize the step-by-step dEdx and dQdx containers
	double step_dedx = 0; 
	std::vector<double> step_dqdx;
	step_dqdx.resize(3);

	//Iterate through all the energy deposition points
	for(auto const& edep : edeps){
	  // 'x_0' definition
	  TVector3 x_0(edep.pos._x, edep.pos._y, edep.pos._z);
	  // 'A' definition 
	  TVector3 A(step_trk.Position().X() - x_0.X(), 
		     step_trk.Position().Y() - x_0.Y(), 
		     step_trk.Position().Z() - x_0.Z());
	  
	  // Distance from the line connecting x_1 and x_2
	  double LineDist = 0;
	 
	  if(B.Mag2() != 0){
	    LineDist = sqrt(A.Mag2() - 2*pow(A*B,2)/B.Mag2() + pow(A*B,2)/B.Mag2());	   
	  }
	  else{LineDist = 0;}
	  
	  //Planar Distance and Radial Line Distance Cuts 
	  // Add in a voxel before and after to account for MCSteps
	  // the line distance allows for 1mm GEANT multiple columb scattering correction, 
	  // small compared to average MCStep-to-MCStep distance
	  if( (a*edep.pos._x + b*edep.pos._y + c*edep.pos._z + d)/sqrt( pow(a,2) + pow(b,2) + pow(c,2)) <= dist + 0.03 &&
	      (a*edep.pos._x + b*edep.pos._y + c*edep.pos._z + d)/sqrt( pow(a,2) + pow(b,2) + pow(c,2)) >=    0 - 0.03 && 
	      LineDist < 0.1){

	    //dEdx Calculation 
	    int npid = 0;
	    double engy = 0;
	    
	    for(auto const& pid_energy : edep.energy){
	      engy += pid_energy.second;
	      npid++;
	    }

	    if(npid != 0){
	      engy /= npid;}
	    else{engy = 0;}
	    
	    step_dedx += engy;
	    
	    // dQdx Calculation
	    auto q_iter = edep.charge.find(edep.pid);
	    if(q_iter != edep.charge.end()){
	      step_dqdx[edep.pid.Plane] += (double)((*q_iter).second);	  	
	    }
	    
	  }
	  
	}
	
	// Normalize to the 3D distance between the MCSteps 

	//Disregard any energy deposition when 2 MCSteps are separated less than the voxel size
	if(dist > 0.03){
	  step_dedx /= dist;
	  step_dqdx[0] /= dist;
	  step_dqdx[1] /= dist;
	  step_dqdx[2] /= dist;	
	}
	else{
	  step_dedx = 0;
	  step_dqdx[0] = 0;
	  step_dqdx[1] = 0;
	  step_dqdx[2] = 0;	
	}
	
	// Build the vector(s) to add to data product
	dEdx.push_back(step_dedx);
	dQdx[0].push_back(step_dqdx[0]);
	dQdx[1].push_back(step_dqdx[1]);
	dQdx[2].push_back(step_dqdx[2]);
	
	

      }

      //Add calorimetry to the data product
      mini_track.dEdx(dEdx);
      mini_track.dQdx(dQdx);
      

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
