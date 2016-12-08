////////////////////////////////////////////////////////////////////////
// $Id: PhotonBackTracker_service.cc,v 1.3 2011/12/13 05:57:02 bckhouse Exp $
//
//
// \file: PhotonBackTracker_service.cc
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#include <map>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace cheat{

  //----------------------------------------------------------------------
  PhotonBackTracker::PhotonBackTracker(const fhicl::ParameterSet& pset,
			   art::ActivityRegistry& reg)
  {
    reconfigure(pset);

    reg.sPreProcessEvent.watch(this, &PhotonBackTracker::Rebuild);
  }
  
  //----------------------------------------------------------------------
  PhotonBackTracker::~PhotonBackTracker()
  {
  }
  
  //----------------------------------------------------------------------
  void PhotonBackTracker::reconfigure(const fhicl::ParameterSet& pset)
  {
    fG4ModuleLabel          = pset.get<std::string>("G4ModuleLabel",            "largeant");
    fMinOpHitEnergyFraction = pset.get<double     >("MinimumOpHitEnergyFraction", 0.1);
    fDelay                  = pset.get< double > ("Delay");
    std::cout<<"Delay set to be :"<<fDelay<<"\n";
  }

  //----------------------------------------------------------------------
  void PhotonBackTracker::Rebuild(const art::Event& evt)
  {
    // do nothing if this is data
    if(evt.isRealData()) return;

    // get the particles from the event
    art::Handle<std::vector<simb::MCParticle>> pHandle;
    evt.getByLabel(fG4ModuleLabel, pHandle);

    // first check to see if we got called to early, that is the particles are
    // and sim channels are not made yet in a MC production job
    // if that is the case, we'll take care of it later
    if(pHandle.failedToGet()){
      mf::LogWarning("PhotonBackTracker") << "failed to get handle to simb::MCParticle from "
				    << fG4ModuleLabel
				    << ", return";
      return;
    }

    // Clear out anything remaining from previous calls to Rebuild
    fParticleList.clear();
    fMCTruthList .clear();
    cOpDetBacktrackerRecords .clear();
    //fVoxelList   .clear();

    art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);

    if( fo.isValid() ){
      for(size_t p = 0; p < pHandle->size(); ++p){
	
	simb::MCParticle *part = new simb::MCParticle(pHandle->at(p));
	fParticleList.Add(part);
	
	// get the simb::MCTruth associated to this sim::ParticleList
	try{
	  art::Ptr<simb::MCTruth> mct = fo.at(p);
	  if(fMCTruthList.size() < 1) fMCTruthList.push_back(mct);
	  else{
	    // check that we are not adding a simb::MCTruth twice to the collection
	    // we know that all the particles for a given simb::MCTruth are put into the
	    // collection of particles at the same time, so we can just check that the 
	    // current art::Ptr has a different id than the last one put 
	    if(!(mct == fMCTruthList.back())) fMCTruthList.push_back(mct);
	  }
	  // fill the track id to mctruth index map
	  fTrackIDToMCTruthIndex[pHandle->at(p).TrackId()] = fMCTruthList.size() - 1;
	}
	catch(cet::exception &ex){
	  mf::LogWarning("PhotonBackTracker") << "unable to find MCTruth from ParticleList "
					<< "created in " << fG4ModuleLabel << " " 
					<< "any attempt to get the MCTruth objects from "
					<< "the backtracker will fail\n"
					<< "message from caught exception:\n" << ex;
	}	
      }// end loop over particles to get MCTruthList  
    }// end if fo.isValid()

    // grab the sim::OpDetBacktrackerRecords for this event
    evt.getView(fG4ModuleLabel, cOpDetBacktrackerRecords);

    // grab the voxel list for this event
    //fVoxelList = sim::SimListUtils::GetLArVoxelList(evt, fG4ModuleLabel);

    fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);

    LOG_DEBUG("PhotonBackTracker") << "PhotonBackTracker has " << cOpDetBacktrackerRecords.size()
			     << " sim::OpDetBacktrackerRecords and " << GetSetOfTrackIDs().size()
			     << " tracks.  The particles are:\n"
			     << fParticleList
			     << "\n the MCTruth information is\n";
    for(size_t mc = 0; mc < fMCTruthList.size(); ++mc)
      LOG_DEBUG("PhotonBackTracker") << *(fMCTruthList.at(mc).get());
    
    return;
  }

  //----------------------------------------------------------------------
  const simb::MCParticle* PhotonBackTracker::TrackIDToParticle(int const& id) const
  {
    sim::ParticleList::const_iterator part_it = fParticleList.find(id);

    if(part_it == fParticleList.end()){
      mf::LogWarning("PhotonBackTracker") << "can't find particle with track id "
				    << id << " in sim::ParticleList"
				    << " returning null pointer";
      return 0;
    }

    return part_it->second;
  }

  //----------------------------------------------------------------------
  const simb::MCParticle* PhotonBackTracker::TrackIDToMotherParticle(int const& id) const
  {
    // get the mother id from the particle navigator
    // the EveId was adopted in the Rebuild method
 
    return this->TrackIDToParticle(fParticleList.EveId(abs(id)));
  }

  //----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& PhotonBackTracker::TrackIDToMCTruth(int const& id) const
  {
    // find the entry in the MCTruth collection for this track id
    size_t mct = fTrackIDToMCTruthIndex.find(abs(id))->second;

    if(/* mct < 0 || */ mct > fMCTruthList.size() ) 
      throw cet::exception("PhotonBackTracker") << "attempting to find MCTruth index for "
					  << "out of range value: " << mct
					  << "/" << fMCTruthList.size() << "\n";

    return fMCTruthList[mct];
  }

  //----------------------------------------------------------------------
  std::vector<sim::SDP> PhotonBackTracker::TrackIDToSimSDP(int const& id) const
  {
    std::vector<sim::SDP> sdps;

    // loop over all sim::OpDetBacktrackerRecords and fill a vector
    // of sim::SDP objects for the given track id
    for(size_t sc = 0; sc < cOpDetBacktrackerRecords.size(); ++sc){
      const auto & pdTimeSDPmap = cOpDetBacktrackerRecords[sc]->timePDclockSDPsMap();
      
      // loop over the SDPMAP      
      for(auto mapitr = pdTimeSDPmap.begin(); mapitr != pdTimeSDPmap.end(); mapitr++){
	
	// loop over the vector of SDP objects.
	const std::vector<sim::SDP>& sdpvec = (*mapitr).second;
	for(size_t iv = 0; iv < sdpvec.size(); ++iv){ 
	  if( abs(sdpvec[iv].trackID) == id) sdps.push_back(sdpvec[iv]);
	}

      } // end loop over map from sim::OpDetBacktrackerRecord
    } // end loop over sim::OpDetBacktrackerRecords
      
    return sdps;
  }

  //----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& PhotonBackTracker::ParticleToMCTruth(const simb::MCParticle* p) const
  {
    return this->TrackIDToMCTruth(p->TrackId());
  }

  //----------------------------------------------------------------------
  std::vector<const simb::MCParticle*> PhotonBackTracker::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
  {
    std::vector<const simb::MCParticle*> ret;
    
    // sim::ParticleList::value_type is a pair (track ID, particle pointer)
    for (const sim::ParticleList::value_type& TrackIDpair: fParticleList) {
      if( TrackIDToMCTruth(TrackIDpair.first) == mct )
        ret.push_back(TrackIDpair.second);
    }

    return ret;
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToTrackID(art::Ptr<recob::OpHit> const& opHit)
  {
    std::vector<sim::TrackSDP> trackSDPs;
    const double pTime = opHit->PeakTime();
    const double pWidth= opHit->Width();
    const double start = (pTime-pWidth)*1000-fDelay;
    const double end   = (pTime+pWidth)*1000-fDelay;

	  this->ChannelToTrackSDPs(trackSDPs, opHit->OpChannel(), start, end);

    return trackSDPs;
  }
  
  //----------------------------------------------------------------------
  const std::vector<std::vector<art::Ptr<recob::OpHit>>> PhotonBackTracker::TrackIDsToOpHits(std::vector<art::Ptr<recob::OpHit>> const& allOpHits, 
										   std::vector<int> const& tkIDs)
  {
    // returns a subset of the opHits in the allOpHits collection that are matched
    // to MC particles listed in tkIDs
    
    // temporary vector of TrackIDs and Ptrs to opHits so only one
    // loop through the (possibly large) allOpHits collection is needed
    std::vector<std::pair<int, art::Ptr<recob::OpHit>>> opHitList;
    std::vector<sim::TrackSDP> tids;
    for(auto itr = allOpHits.begin(); itr != allOpHits.end(); ++itr) {
      tids.clear();
      art::Ptr<recob::OpHit> const& opHit = *itr;
      const double pTime = opHit->PeakTime(), pWidth= opHit->Width();
      const double start = (pTime-pWidth)*1000.0-fDelay, end = (pTime+pWidth)*1000.0-fDelay;
      this->ChannelToTrackSDPs(tids, opHit->OpChannel(), start, end);
      for(auto itid = tids.begin(); itid != tids.end(); ++itid) {
        for(auto itkid = tkIDs.begin(); itkid != tkIDs.end(); ++itkid) {
          if(itid->trackID == *itkid) {
            if(itid->energyFrac > fMinOpHitEnergyFraction) 
              opHitList.push_back(std::make_pair(*itkid, opHit));
          }
        } // itkid 
      } // itid
    } // itr
    
    // now build the truOpHits vector that will be returned to the caller
    std::vector<std::vector<art::Ptr<recob::OpHit>>> truOpHits;
    // temporary vector containing opHits assigned to one MC particle
    std::vector<art::Ptr<recob::OpHit>> tmpOpHits;
    for(auto itkid = tkIDs.begin(); itkid != tkIDs.end(); ++itkid) {
      tmpOpHits.clear();
      for(auto itr = opHitList.begin(); itr != opHitList.end(); ++itr) {
        if(*itkid == (*itr).first) tmpOpHits.push_back((*itr).second);
      }
      truOpHits.push_back(tmpOpHits);
    }
    
    return truOpHits;
  }

  //----------------------------------------------------------------------
  // plist is assumed to have adopted the appropriate EveIdCalculator prior to 
  // having been passed to this method. It is likely that the EmEveIdCalculator is
  // the one you always want to use
  std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToEveID(art::Ptr<recob::OpHit> const& opHit)
  {
    std::vector<sim::TrackSDP> trackSDPs = this->OpHitToTrackID(opHit);

    // make a map of evd ID values and fraction of energy represented by
    // that eve id in this opHit
    std::map<int, float> eveIDtoEfrac;
    
    double totalE = 0.;
    for(size_t t = 0; t < trackSDPs.size(); ++t){
      eveIDtoEfrac[fParticleList.EveId( trackSDPs[t].trackID )] += trackSDPs[t].energy;
      totalE += trackSDPs[t].energy;
    }
    
    // now fill the eveSDPs vector from the map
    std::vector<sim::TrackSDP> eveSDPs;
    eveSDPs.reserve(eveIDtoEfrac.size());
    for(auto itr = eveIDtoEfrac.begin(); itr != eveIDtoEfrac.end(); itr++){
      sim::TrackSDP temp;
      temp.trackID    = (*itr).first;
      temp.energyFrac = (*itr).second/totalE;
      temp.energy     = (*itr).second;
      eveSDPs.push_back(std::move(temp));
    }

    return eveSDPs;
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfEveIDs()
  {
    std::set<int> eveIDs;

    sim::ParticleList::const_iterator plitr = fParticleList.begin();
    while(plitr != fParticleList.end() ){
      int eveID = fParticleList.EveId((*plitr).first);
      // look to see if this eveID is already in the set
      if( eveIDs.find(eveID) == eveIDs.end() ) eveIDs.insert(eveID);
      plitr++;
    }

    return eveIDs;
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfTrackIDs()
  {
    // fParticleList::value_type is a pair (track, particle pointer)
    std::set<int> trackIDs;
    for (const sim::ParticleList::value_type& pl: fParticleList)
      trackIDs.insert(pl.first);

    return trackIDs;
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfEveIDs(std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    std::set<int> eveIDs;

    std::vector< art::Ptr<recob::OpHit> >::const_iterator itr = opHits.begin();
    while(itr != opHits.end() ){
      
      // get the eve ids corresponding to this opHit
      const std::vector<sim::TrackSDP> sdps = OpHitToEveID(*itr);
      
      // loop over the ides and extract the track ids
      for(size_t i = 0; i < sdps.size(); ++i) eveIDs.insert(sdps[i].trackID);

      itr++;
    }

    return eveIDs;
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfTrackIDs(std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    std::set<int> trackIDs;

    std::vector< art::Ptr<recob::OpHit> >::const_iterator itr = opHits.begin();
    while(itr != opHits.end() ){
      
      std::vector<sim::TrackSDP> trackSDPs;
     
      // get the track ids corresponding to this opHit
      const double pTime = (*itr)->PeakTime();
      const double pWidth= (*itr)->Width();
      const double start = (pTime-pWidth)*1000.0-fDelay;
      const double end   = (pTime+pWidth)*1000.0-fDelay;
     // const double start = (*itr)->PeakTimeMinusRMS();
     // const double end   = (*itr)->PeakTimePlusRMS();

      this->ChannelToTrackSDPs(trackSDPs, (*itr)->OpChannel(), start, end);
      
      // loop over the ides and extract the track ids
      for(size_t i = 0; i < trackSDPs.size(); ++i) {
	trackIDs.insert(trackSDPs[i].trackID);
      }
   
      itr++;
    }

    return trackIDs;
  }

  //----------------------------------------------------------------------
  double PhotonBackTracker::OpHitCollectionPurity(std::set<int>                              trackIDs, 
					  std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection
    // if the EveID shows up in the input list of trackIDs, then it counts
    float total   = 1.*opHits.size();;
    float desired = 0.;

    // don't have to check the view in the opHits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < opHits.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackID(opHit);

      // don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
	if(trackIDs.find(opHitTrackIDs[e].trackID) != trackIDs.end()){
	  desired += 1.;
	  break;
	}
      }

    }// end loop over opHits

    double purity = 0.;
    if(total > 0) purity = desired/total;

    return purity;
  }

  //----------------------------------------------------------------------
  double PhotonBackTracker::OpHitChargeCollectionPurity(std::set<int>                              trackIDs, 
						std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection
    // if the EveID shows up in the input list of trackIDs, then it counts
    float total   = 0;
    float desired = 0.;

    // don't have to check the view in the opHits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < opHits.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackID(opHit);
       
      total+=opHit->Area(); // sum up the charge in the cluster

      // don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
	if(trackIDs.find(opHitTrackIDs[e].trackID) != trackIDs.end()){
	  desired += opHit->Area();
	  break;
	}
      }

    }// end loop over opHits

    double purity = 0.;
    if(total > 0) purity = desired/total;

    return purity;
  }


  //----------------------------------------------------------------------
  double PhotonBackTracker::OpHitCollectionEfficiency(std::set<int>                              trackIDs, 
					      std::vector< art::Ptr<recob::OpHit> > const& opHits,
					      std::vector< art::Ptr<recob::OpHit> > const& allOpHits,
					      geo::View_t                         const& view)
  {
     throw cet::exception("PhotonBackTracker")<<"This function is not supported. OpHits do not have type View.\n";
  }

  double PhotonBackTracker::OpHitCollectionEfficiency(std::set<int>                              trackIDs, 
					      std::vector< art::Ptr<recob::OpHit> > const& opHits,
					      std::vector< art::Ptr<recob::OpHit> > const& allOpHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection
    // and the energy associated with the desired trackID
    float desired = 0.;
    float total   = 0.;

    // don't have to check the view in the opHits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < opHits.size(); ++h){

      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackID(opHit);

      // don't worry about opHits where the energy fraction for the chosen
      // trackID is < 0.1
      // also don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
	if(trackIDs.find(opHitTrackIDs[e].trackID) != trackIDs.end() &&
	   opHitTrackIDs[e].energyFrac             >= fMinOpHitEnergyFraction){
	  desired += 1.;
	  break;
	}
      }
    }// end loop over opHits

    // now figure out how many opHits in the whole collection are associated with this id
    for(size_t h = 0; h < allOpHits.size(); ++h){

      art::Ptr<recob::OpHit> opHit = allOpHits[h];

      // check that we are looking at the appropriate view here
      // in the case of 3D objects we take all opHits
      //if(opHit->View() != view && view != geo::k3D ) continue;

      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackID(opHit);

      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
	// don't worry about opHits where the energy fraction for the chosen
	// trackID is < 0.1
	// also don't double count if this opHit has more than one of the
	// desired track IDs associated with it
	if(trackIDs.find(opHitTrackIDs[e].trackID) != trackIDs.end() &&
	   opHitTrackIDs[e].energyFrac             >= fMinOpHitEnergyFraction){
	  total += 1.;
	  break;
	}
      }

    }// end loop over all opHits
    
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;

    return efficiency;
  }

  //----------------------------------------------------------------------
  double PhotonBackTracker::OpHitChargeCollectionEfficiency(std::set<int>                              trackIDs, 
						    std::vector< art::Ptr<recob::OpHit> > const& opHits,
						    std::vector< art::Ptr<recob::OpHit> > const& allOpHits,
						    geo::View_t                         const& view)
  {
     throw cet::exception("PhotonBackTracker")<<"This function is not supported. OpHits do not have type View.\n";
  }
  double PhotonBackTracker::OpHitChargeCollectionEfficiency(std::set<int>                              trackIDs, 
						    std::vector< art::Ptr<recob::OpHit> > const& opHits,
						    std::vector< art::Ptr<recob::OpHit> > const& allOpHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection
    // and the energy associated with the desired trackID
    float desired = 0.;
    float total   = 0.;

    // don't have to check the view in the opHits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < opHits.size(); ++h){

      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackID(opHit);

      // don't worry about opHits where the energy fraction for the chosen
      // trackID is < 0.1
      // also don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
	if(trackIDs.find(opHitTrackIDs[e].trackID) != trackIDs.end() &&
	   opHitTrackIDs[e].energyFrac             >= fMinOpHitEnergyFraction){
	  desired += opHit->Area();
	  break;
	}
      }
    }// end loop over opHits

    // now figure out how many opHits in the whole collection are associated with this id
    for(size_t h = 0; h < allOpHits.size(); ++h){

      art::Ptr<recob::OpHit> opHit = allOpHits[h];

      // check that we are looking at the appropriate view here
      // in the case of 3D objects we take all opHits
      //if(opHit->View() != view && view != geo::k3D ) continue;

      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackID(opHit);

      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
	// don't worry about opHits where the energy fraction for the chosen
	// trackID is < 0.1
	// also don't double count if this opHit has more than one of the
	// desired track IDs associated with it
	if(trackIDs.find(opHitTrackIDs[e].trackID) != trackIDs.end() &&
	   opHitTrackIDs[e].energyFrac             >= fMinOpHitEnergyFraction){
	  total += opHit->Area();
	  break;
	}
      }

    }// end loop over all opHits
    
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;

    return efficiency;
  }



  //----------------------------------------------------------------------
  const sim::OpDetBacktrackerRecord* PhotonBackTracker::FindOpDetBacktrackerRecord(int opDetNum) const
  {
    const sim::OpDetBacktrackerRecord* opDet = 0;

    for(size_t sc = 0; sc < cOpDetBacktrackerRecords.size(); ++sc){
      if(cOpDetBacktrackerRecords[sc]->OpDetNum() == opDetNum) opDet = cOpDetBacktrackerRecords[sc];
    }

    if(!opDet)
      throw cet::exception("PhotonBackTracker") << "No sim::OpDetBacktrackerRecord corresponding "
					  << "to opDetNum: " << opDetNum << "\n";

    return opDet;
  }

  //----------------------------------------------------------------------
  void PhotonBackTracker::ChannelToTrackSDPs(std::vector<sim::TrackSDP>&   trackSDPs,
				     int channel,
				     const double opHit_start_time,
				     const double opHit_end_time)
  {
    trackSDPs.clear();

    double totalE = 0.;

    try{
      const sim::OpDetBacktrackerRecord* schannel = this->FindOpDetBacktrackerRecord(channel);
      
      // loop over the photons in the channel and grab those that are in time 
      // with the identified opHit start and stop times
      const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
      int start_tdc = ts->OpticalG4Time2TDC( opHit_start_time );
      int end_tdc   = ts->OpticalG4Time2TDC( opHit_end_time   );
      if(start_tdc<0) start_tdc = 0;
      if(end_tdc<0) end_tdc = 0;
      std::vector<sim::SDP> simSDPs = schannel->TrackIDsAndEnergies(start_tdc, end_tdc);
      
      // first get the total energy represented by all track ids for 
      // this channel and range of tdc values
      for(size_t e = 0; e < simSDPs.size(); ++e)
	totalE += simSDPs[e].energy;
      
      
      // protect against a divide by zero below
      if(totalE < 1.e-5) totalE = 1.;
      
      // loop over the entries in the map and fill the input vectors
      
      for(size_t e = 0; e < simSDPs.size(); ++e){
	
	if(simSDPs[e].trackID == sim::NoParticleId) continue;
	
	sim::TrackSDP info;
	info.trackID    = simSDPs[e].trackID;
	info.energyFrac = simSDPs[e].energy/totalE;
	info.energy     = simSDPs[e].energy;
	
	trackSDPs.push_back(info);
	
      }
    }// end try
    catch(cet::exception e){
      mf::LogWarning("PhotonBackTracker") << "caught exception \n"
				    << e;
    }

    return;
  }

  //----------------------------------------------------------------------
  void PhotonBackTracker::OpHitToSimSDPs(recob::OpHit const& opHit,
                                 std::vector<sim::SDP>&      sdps) const
  {
    // Get services.
    const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    
    double fPeakTime       = opHit.PeakTime();
    double fWidth          = opHit.Width();
    sim::OpDetBacktrackerRecord::timePDclock_t start_time = ((fPeakTime-fWidth)*1000.0)-fDelay;
    sim::OpDetBacktrackerRecord::timePDclock_t end_time = ((fPeakTime+fWidth)*1000.0)-fDelay;
    
    sdps = FindOpDetBacktrackerRecord( geom->OpDetFromOpChannel(opHit.OpChannel()) )->TrackIDsAndEnergies(start_time, end_time);

  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTracker::SimSDPsToXYZ(std::vector<sim::SDP> const& sdps)
  {
    std::vector<double> xyz(3, -999.);

    double x = 0.;
    double y = 0.;
    double z = 0.;
    double w = 0.;

    // loop over photons.

    for(auto const& sdp : sdps) {

      double weight = sdp.numPhotons;
      
      w += weight;
      x += weight * sdp.x;
      y += weight * sdp.y;
      z += weight * sdp.z;

    }// end loop over sim::SDPs
	
    // if the sum of the weights is still 0, then return
    // the obviously stupid default values
    if(w < 1.e-5)
      throw cet::exception("PhotonBackTracker") << "No sim::SDPs providing non-zero number of photons"
					  << " can't determine originating location from truth\n";

    xyz[0] = x/w;
    xyz[1] = y/w;
    xyz[2] = z/w;

    return xyz;
  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTracker::OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit)
  {
    std::vector<sim::SDP> sdps;
    OpHitToSimSDPs(opHit, sdps);
    return SimSDPsToXYZ(sdps);
  }

} // namespace

namespace cheat{
  DEFINE_ART_SERVICE(PhotonBackTracker)
}
