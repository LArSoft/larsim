////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCSHOWERRECOALG_CXX
#define MCSHOWERRECOALG_CXX

#include "MCShowerRecoAlg.h"

namespace sim {

  //##################################################################
  MCShowerRecoAlg::MCShowerRecoAlg(fhicl::ParameterSet const& pset) 
    : fPartAlg(pset.get< fhicl::ParameterSet >("MCShowerRecoPart"))
  //##################################################################
  {
    fDebugMode = pset.get<bool>("DebugMode");
    fMinShowerEnergy = pset.get<double>("MinShowerEnergy");
    fMinNumDaughters = pset.get<unsigned int>("MinNumDaughters");
  }

  void MCShowerRecoAlg::Reconstruct(const MCRecoPart& part_v,
				    const MCRecoEdep& edep_v)
  {
    
    art::ServiceHandle<geo::Geometry> geo;

    fPartAlg.ConstructShower(part_v);

    fMCShower.clear();

    // Get shower info from grouped particles
    const std::vector<unsigned int> shower_index_v = fPartAlg.ShowerMothers();
    fMCShower.reserve(shower_index_v.size());
    std::map<size_t,int> stored_mcs_index;
    for(size_t shower_index = 0; shower_index < shower_index_v.size(); ++shower_index) {

      unsigned int shower_candidate = shower_index_v.at(shower_index);
      auto const&  shower_part      = part_v.at(shower_candidate);

      unsigned int mother_track   = part_v.MotherTrackID(shower_candidate);
      unsigned int ancestor_track = part_v.AncestorTrackID(shower_candidate);

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

      auto mcs_index_iter = stored_mcs_index.insert(std::make_pair(shower_index,-1));

      double shower_g4_energy = shower_part._start_mom[3];

      if(fDebugMode)

	std::cout << "Found MCShower with mother energy: " << shower_g4_energy << " MeV";

      // Skip if mother energy is less than the enery threshold
      if(shower_g4_energy < fMinShowerEnergy) {
	if(fDebugMode)
	  std::cout << " ... below energy threshold: skipping!"<<std::endl;
	continue;
      }else if(shower_part._daughters.size() < fMinNumDaughters) {
	if(fDebugMode)
	  std::cout << " ... below # daughter particle count threshold: skipping!"<<std::endl;
	continue;
      }else if(fDebugMode) {
	std::cout << " ... condition matched. Storing this MCShower..."<<std::endl;
      }

      // Record this MCShower

      (*mcs_index_iter.first).second = fMCShower.size();

      if(fDebugMode)
	
	std::cout << " Storage index " << fMCShower.size() << " => Shower index " << shower_index
		  << std::endl;

      ::sim::MCShower shower_prof;

      shower_prof.Origin  ( shower_part._origin   );
      shower_prof.PdgCode ( shower_part._pdgcode  );
      shower_prof.TrackID ( shower_part._track_id );
      shower_prof.Process ( shower_part._process  );

      shower_prof.MotherPdgCode ( mother_part._pdgcode  );
      shower_prof.MotherTrackID ( mother_part._track_id );
      shower_prof.MotherProcess ( mother_part._process  );

      shower_prof.AncestorPdgCode ( ancestor_part._pdgcode  );
      shower_prof.AncestorTrackID ( ancestor_part._track_id );
      shower_prof.AncestorProcess ( ancestor_part._process  );

      shower_prof.Start         ( MCStep ( shower_part._start_vtx, shower_part._start_mom ) );
      shower_prof.End           ( MCStep ( shower_part._end_vtx,   shower_part._end_mom   ) );
      shower_prof.MotherStart   ( MCStep ( mother_part._start_vtx, mother_part._start_mom ) );
      shower_prof.MotherEnd     ( MCStep ( mother_part._end_vtx,   mother_part._end_mom   ) );
      shower_prof.AncestorStart ( MCStep ( mother_part._start_vtx, mother_part._start_mom ) );
      shower_prof.AncestorEnd   ( MCStep ( mother_part._end_vtx,   mother_part._end_mom   ) );

      std::vector<unsigned int> daughter_track_id;
      daughter_track_id.reserve( fPartAlg.ShowerDaughters(shower_index).size() );

      for(auto const& index : fPartAlg.ShowerDaughters(shower_index))

	daughter_track_id.push_back( part_v.at(index)._track_id );

      shower_prof.DaughterTrackID(daughter_track_id);

      fMCShower.push_back(shower_prof);
    }

    if(fDebugMode)
      std::cout << " Found " << fMCShower.size() << " MCShowers. Now computing DetProfile position..." << std::endl;
    
    // Next, loop over deposited energy and find the DetProfile vtx point
    std::vector<TLorentzVector>   mcs_daughter_vtx_v (fMCShower.size(),TLorentzVector());
    std::vector<double>           mcs_min_dist_v     (fMCShower.size(),kINVALID_DOUBLE);

    std::vector<std::vector<double> > shower_g4_dir(fMCShower.size(),std::vector<double>(3,0));
    for(size_t shower_index = 0; shower_index < fMCShower.size(); ++shower_index) {
    
      auto& g4_dir = shower_g4_dir[shower_index];
      auto const& g4_mom = fMCShower[shower_index].Start().Momentum();
      double magnitude = 0;
      for(size_t i=0; i<3; ++i) {
	g4_dir[i] = g4_mom[i];
	magnitude += pow(g4_mom[i],2);
      }
      magnitude = sqrt(magnitude);
      if(magnitude > 1.e-10)
	for(auto& v : g4_dir) v /= magnitude;
      
    }

    std::map<unsigned int,size_t> edep_index_map     (edep_v.TrackIndexMap());
    for(auto track_iter = edep_index_map.begin();
	track_iter != edep_index_map.end();
	++track_iter) {

      unsigned int part_track_id = (*track_iter).first;
      unsigned int edep_index    = (*track_iter).second;

      unsigned int part_index = part_v.TrackToParticleIndex(part_track_id);
      if(part_index == kINVALID_UINT) continue;

      int shower_index = fPartAlg.ShowerIndex(part_index);
      if(shower_index < 0 || shower_index == kINVALID_INT) continue;

      auto mcs_index_iter = stored_mcs_index.find(shower_index);
      if(mcs_index_iter == stored_mcs_index.end()) {
	std::cerr<<"Logic error: invalid index of stored MCShower!"<<std::endl;
	throw std::exception();
      }

      int stored_mcs_index = (*mcs_index_iter).second;
      if((*mcs_index_iter).second < 0) continue;

      auto const& shower_vtx = fMCShower[stored_mcs_index].Start().Position();
      auto const& shower_dir = shower_g4_dir[stored_mcs_index];
      auto& daughter_vtx = mcs_daughter_vtx_v[stored_mcs_index];

      for(auto const& edep : edep_v.GetEdepArrayAt((size_t)edep_index)) {

	std::vector<double> shower_dep_dir(3,0);
	shower_dep_dir[0] = edep.x - shower_vtx[0];
	shower_dep_dir[1] = edep.y - shower_vtx[1];
	shower_dep_dir[2] = edep.z - shower_vtx[2];

	double dist = sqrt( pow(shower_dep_dir[0],2) + pow(shower_dep_dir[1],2) + pow(shower_dep_dir[2],2) );
	for(auto& v : shower_dep_dir) v /= dist;

	double angle = acos( shower_dep_dir[0] * shower_dir[0] +
			     shower_dep_dir[1] * shower_dir[1] +
			     shower_dep_dir[2] * shower_dir[2] ) / TMath::Pi() * 180.;

	if(dist < mcs_min_dist_v[stored_mcs_index] && angle < 10) {
	  
	  mcs_min_dist_v[stored_mcs_index] = dist;
	  daughter_vtx[0] = edep.x;
	  daughter_vtx[1] = edep.y;
	  daughter_vtx[2] = edep.z;
	  daughter_vtx[3] = (dist/100. / 2.998e8)*1.e9 + shower_vtx[3];
	  
	}
      }
    }

    if(fDebugMode)
      std::cout << " Found " << fMCShower.size() << " MCShowers. Now computing DetProfile momentum..." << std::endl;

    // Next, loop over deposited energy and group them into showers
    //TLorentzVector mom_default(0,0,0,0);
    std::vector<TLorentzVector> mcs_daughter_mom_v   ( fMCShower.size(), TLorentzVector() );
    std::vector<std::vector<double> > plane_charge_v ( fMCShower.size(), std::vector<double>(geo->Nplanes(),0) );
    for(auto track_iter = edep_index_map.begin();
	track_iter != edep_index_map.end();
	++track_iter) {

      unsigned int part_track_id = (*track_iter).first;
      unsigned int edep_index    = (*track_iter).second;

      unsigned int part_index = part_v.TrackToParticleIndex(part_track_id);
      if(part_index == kINVALID_UINT) continue;

      int shower_index = fPartAlg.ShowerIndex(part_index);
      if(shower_index < 0 || shower_index == kINVALID_INT) continue;

      auto mcs_index_iter = stored_mcs_index.find(shower_index);
      if(mcs_index_iter == stored_mcs_index.end()) {
	std::cerr<<"Logic error: invalid index of stored MCShower!"<<std::endl;
	throw std::exception();
      }

      int stored_mcs_index = (*mcs_index_iter).second;
      if((*mcs_index_iter).second < 0) continue;

      auto const& daughter_vtx = mcs_daughter_vtx_v.at(stored_mcs_index);
      auto& daughter_mom       = mcs_daughter_mom_v.at(stored_mcs_index);
      auto& plane_charge       = plane_charge_v.at(stored_mcs_index);

      std::vector<double> vtx(3,0);
      std::vector<double> mom(3,0);
      //std::cout<<"Edep array @ "<<edep_index << " ... " << edep_v.GetEdepArrayAt((size_t)edep_index).size()<<" entries..."<<std::endl;
      for(auto const& edep : edep_v.GetEdepArrayAt((size_t)edep_index)) {

	vtx.at(0) = edep.x;
	vtx.at(1) = edep.y;
	vtx.at(2) = edep.z;

	// Compute unit vector to this energy deposition
	mom.at(0) = vtx.at(0) - daughter_vtx[0];
	mom.at(1) = vtx.at(1) - daughter_vtx[1];
	mom.at(2) = vtx.at(2) - daughter_vtx[2];

	// Weight by energy (momentum)
	double magnitude = sqrt(pow(mom.at(0),2) + pow(mom.at(1),2) + pow(mom.at(2),2));
	if(magnitude>1.e-10) {
	  mom.at(0) = mom.at(0) * (edep.energy/magnitude);
	  mom.at(1) = mom.at(1) * (edep.energy/magnitude);
	  mom.at(2) = mom.at(2) * (edep.energy/magnitude);
	  daughter_mom[0] += mom.at(0);
	  daughter_mom[1] += mom.at(1);
	  daughter_mom[2] += mom.at(2);
	}

	daughter_mom[3] += edep.energy;

	// Charge
	plane_charge.at(geo::kU) += edep.qU;
	plane_charge.at(geo::kV) += edep.qV;
	plane_charge.at(geo::kW) += edep.qW;

      }
    }
    
    if(fDebugMode)
      std::cout << " Found " << fMCShower.size() << " MCShowers. Now storing..." << std::endl;
    
    // Store plane charge & daughter momentum
    for(size_t mcs_index=0; mcs_index<fMCShower.size(); ++mcs_index) {
      
      auto& daughter_vtx = mcs_daughter_vtx_v[mcs_index];
      auto& daughter_mom = mcs_daughter_mom_v[mcs_index];
      
      // Correct for energy deposition normalization
      double magnitude = sqrt(pow(daughter_mom[0],2)+pow(daughter_mom[1],2)+pow(daughter_mom[2],2));

      if(magnitude>1.e-10) {
	daughter_mom[0] *= (daughter_mom[3]/magnitude);
	daughter_mom[1] *= (daughter_mom[3]/magnitude);
	daughter_mom[2] *= (daughter_mom[3]/magnitude);
      }else
	for(size_t i=0; i<4; ++i) daughter_mom[i]=0;
      
      fMCShower.at(mcs_index).DetProfile( MCStep( daughter_vtx, daughter_mom ) );
      fMCShower.at(mcs_index).Charge(plane_charge_v[mcs_index]);
    }
    
    if(fDebugMode) {
      
      for(auto const& prof : fMCShower) {
	
	std::cout
	  
	  << Form("  Shower particle:     PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.PdgCode(), prof.TrackID(),
		  prof.Start().X(),prof.Start().Y(),prof.Start().Z(),prof.Start().T(),
		  prof.Start().Px(),prof.Start().Py(),prof.Start().Pz(),prof.Start().E())
	  << std::endl
	  << Form("    Mother particle:   PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.MotherPdgCode(), prof.MotherTrackID(),
		  prof.MotherStart().X(),prof.MotherStart().Y(),prof.MotherStart().Z(),prof.MotherStart().T(),
		  prof.MotherStart().Px(),prof.MotherStart().Py(),prof.MotherStart().Pz(),prof.MotherStart().E())
	  << std::endl
	  << Form("    Ancestor particle: PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.AncestorPdgCode(), prof.AncestorTrackID(),
		  prof.AncestorStart().X(),prof.AncestorStart().Y(),prof.AncestorStart().Z(),prof.AncestorStart().T(),
		  prof.AncestorStart().Px(),prof.AncestorStart().Py(),prof.AncestorStart().Pz(),prof.AncestorStart().E())
	  << std::endl
	  << Form("    ... with %zu daughters: Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.DaughterTrackID().size(),
		  prof.DetProfile().X(),prof.DetProfile().Y(),prof.DetProfile().Z(),prof.DetProfile().T(),
		  prof.DetProfile().Px(),prof.DetProfile().Py(),prof.DetProfile().Pz(),prof.DetProfile().E())
	  << std::endl		  
	  << "    Charge per plane: ";
	
	for(size_t i=0; i<geo->Nplanes(); ++i) {
	  
	  std::cout << " | Plane " << i << std::flush;
	  std::cout << " ... Q = " << prof.Charge(i) << std::flush;
	  
	}
	std::cout<<std::endl<<std::endl;
      }
    }
  }
}


#endif
