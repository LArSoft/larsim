////////////////////////////////////////////////////////////////////////
// $Id: BackTracker_service.cc,v 1.3 2011/12/13 05:57:02 bckhouse Exp $
//
//
// \file: BackTracker_service.cc
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
#include "MCCheater/BackTracker.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "Simulation/EmEveIdCalculator.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"

namespace cheat{

  //----------------------------------------------------------------------
  BackTracker::BackTracker(const fhicl::ParameterSet& pset,
			   art::ActivityRegistry& reg)
  {
    reconfigure(pset);

    reg.sPreProcessEvent.watch(this, &BackTracker::Rebuild);
  }
  
  //----------------------------------------------------------------------
  BackTracker::~BackTracker()
  {
  }
  
  //----------------------------------------------------------------------
  void BackTracker::reconfigure(const fhicl::ParameterSet& pset)
  {
    fG4ModuleLabel        = pset.get<std::string>("G4ModuleLabel",            "largeant");
    fMinHitEnergyFraction = pset.get<double     >("MinimumHitEnergyFraction", 0.1);
  }

  //----------------------------------------------------------------------
  void BackTracker::Rebuild(const art::Event& evt)
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
      mf::LogWarning("BackTracker") << "failed to get handle to simb::MCParticle from "
				    << fG4ModuleLabel
				    << ", return";
      return;
    }

    // Clear out anything remaining from previous calls to Rebuild
    fTrackIDs    .clear();
    fParticleList.clear();
    fMCTruthList .clear();
    fSimChannels .clear();
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
	  mf::LogWarning("BackTracker") << "unable to find MCTruth from ParticleList "
					<< "created in " << fG4ModuleLabel << " " 
					<< "any attempt to get the MCTruth objects from "
					<< "the backtracker will fail\n"
					<< "message from caught exception:\n" << ex;
	}	
      }// end loop over particles to get MCTruthList  
    }// end if fo.isValid()


    // fill the set of track ids from the ParticleNavigator
    for(size_t t = 0; t < fParticleList.size(); ++t){
      fTrackIDs.insert(fParticleList.TrackId(t));
    }

    // grab the sim::SimChannels for this event
    evt.getView(fG4ModuleLabel, fSimChannels);

    // grab the voxel list for this event
    //fVoxelList = sim::SimListUtils::GetLArVoxelList(evt, fG4ModuleLabel);

    fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);

    LOG_DEBUG("BackTracker") << "BackTracker has " << fSimChannels.size()
			     << " sim::SimChannels and " << fTrackIDs.size()
			     << " tracks.  The particles are:\n"
			     << fParticleList
			     << "\n the MCTruth information is\n";
    for(size_t mc = 0; mc < fMCTruthList.size(); ++mc)
      LOG_DEBUG("BackTracker") << *(fMCTruthList.at(mc).get());
    
    return;
  }

  //----------------------------------------------------------------------
  const simb::MCParticle* BackTracker::TrackIDToParticle(int const& id) const
  {
    sim::ParticleList::const_iterator part_it = fParticleList.find(id);

    if(part_it == fParticleList.end()){
      mf::LogWarning("BackTracker") << "can't find particle with track id "
				    << id << " in sim::ParticleList"
				    << " returning null pointer";
      return 0;
    }

    return part_it->second;
  }

  //----------------------------------------------------------------------
  const simb::MCParticle* BackTracker::TrackIDToMotherParticle(int const& id) const
  {
    // get the mother id from the particle navigator
    // the EveId was adopted in the Rebuild method
 
    return this->TrackIDToParticle(fParticleList.EveId(abs(id)));
  }

  //----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& BackTracker::TrackIDToMCTruth(int const& id) const
  {
    // find the entry in the MCTruth collection for this track id
    size_t mct = fTrackIDToMCTruthIndex.find(abs(id))->second;

    if(mct < 0 || mct > fMCTruthList.size() ) 
      throw cet::exception("BackTracker") << "attempting to find MCTruth index for "
					  << "out of range value: " << mct
					  << "/" << fMCTruthList.size();

    return fMCTruthList[mct];
  }

  //----------------------------------------------------------------------
  const std::vector<sim::IDE> BackTracker::TrackIDToSimIDE(int const& id) const
  {
    std::vector<sim::IDE> ides;

    // loop over all sim::SimChannels and fill a vector
    // of sim::IDE objects for the given track id
    for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
      const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = fSimChannels[sc]->TDCIDEMap();
      
      // loop over the IDEMAP      
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
	
	// loop over the vector of IDE objects.
	const std::vector<sim::IDE> idevec = (*mapitr).second;
	for(size_t iv = 0; iv < idevec.size(); ++iv){ 
	  if( abs(idevec[iv].trackID) == id) ides.push_back(idevec[iv]);
	}

      } // end loop over map from sim::SimChannel
    } // end loop over sim::SimChannels
      
    return ides;
  }

  //----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& BackTracker::ParticleToMCTruth(const simb::MCParticle* p) const
  {
    return this->TrackIDToMCTruth(p->TrackId());
  }

  //----------------------------------------------------------------------
  std::vector<const simb::MCParticle*> BackTracker::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
  {
    std::vector<const simb::MCParticle*> ret;

    for(auto itr = fTrackIDs.begin(); itr != fTrackIDs.end(); itr++){
      if( this->TrackIDToMCTruth(*itr) == mct ) ret.push_back(fParticleList.find(*itr)->second);
    }

    return ret;
  }

  //----------------------------------------------------------------------
  const std::vector<TrackIDE> BackTracker::HitToTrackID(art::Ptr<recob::Hit> const& hit)
  {
    std::vector<TrackIDE> trackIDEs;
    
    double start = hit->StartTime();
    double end   = hit->EndTime();
	
    this->ChannelToTrackID(trackIDEs, hit->Channel(), start, end);

    return trackIDEs;
  }

  //----------------------------------------------------------------------
  // plist is assumed to have adopted the appropriate EveIdCalculator prior to 
  // having been passed to this method. It is likely that the EmEveIdCalculator is
  // the one you always want to use
  const std::vector<TrackIDE> BackTracker::HitToEveID(art::Ptr<recob::Hit> const& hit)
  {
    std::vector<TrackIDE> eveides;
    std::vector<TrackIDE> trackides = this->HitToTrackID(hit);

    // make a map of evd ID values and fraction of energy represented by
    // that eve id in this hit
    std::map<int, float> eveToE;
    
    double totalE = 0.;
    for(size_t t = 0; t < trackides.size(); ++t){
      eveToE[fParticleList.EveId( trackides[t].trackID )] += trackides[t].energy;
      totalE += trackides[t].energy;
    }
    
    // now fill the eveides vector from the map
    for(auto itr = eveToE.begin(); itr != eveToE.end(); itr++){
      TrackIDE temp;
      temp.trackID    = (*itr).first;
      temp.energyFrac = (*itr).second/totalE;
      temp.energy     = (*itr).second;
      eveides.push_back(temp);
    }

    return eveides;
  }

  //----------------------------------------------------------------------
  std::set<int> BackTracker::GetSetOfEveIDs()
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
  std::set<int> BackTracker::GetSetOfTrackIDs()
  {
    std::set<int> trackIDs;

    sim::ParticleList::const_iterator plitr = fParticleList.begin();
    while(plitr != fParticleList.end() ){
      int trackID = (*plitr).first;
      // look to see if this eveID is already in the set
      if( trackIDs.find(trackID) == trackIDs.end() ) trackIDs.insert(trackID);
      plitr++;
    }

    return trackIDs;
  }

  //----------------------------------------------------------------------
  std::set<int> BackTracker::GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits)
  {
    std::set<int> eveIDs;

    std::vector< art::Ptr<recob::Hit> >::const_iterator itr = hits.begin();
    while(itr != hits.end() ){
      
      // get the eve ids corresponding to this hit
      const std::vector<TrackIDE> ides = HitToEveID(*itr);
      
      // loop over the ides and extract the track ids
      for(size_t i = 0; i < ides.size(); ++i) eveIDs.insert(ides[i].trackID);

      itr++;
    }

    return eveIDs;
  }

  //----------------------------------------------------------------------
  std::set<int> BackTracker::GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits)
  {
    std::set<int> trackIDs;

    std::vector< art::Ptr<recob::Hit> >::const_iterator itr = hits.begin();
    while(itr != hits.end() ){
      
      std::vector<TrackIDE> trackIDEs;
     
      // get the track ids corresponding to this hit
      double start = (*itr)->StartTime();
      double end   = (*itr)->EndTime();

      this->ChannelToTrackID(trackIDEs, (*itr)->Channel(), start, end);
      
      // loop over the ides and extract the track ids
      for(size_t i = 0; i < trackIDEs.size(); ++i) {
	trackIDs.insert(trackIDEs[i].trackID);
      }
   
      itr++;
    }

    return trackIDs;
  }

  //----------------------------------------------------------------------
  double BackTracker::HitCollectionPurity(std::set<int>                              trackIDs, 
					  std::vector< art::Ptr<recob::Hit> > const& hits)
  {
    // get the list of EveIDs that correspond to the hits in this collection
    // if the EveID shows up in the input list of trackIDs, then it counts
    float total   = 1.*hits.size();;
    float desired = 0.;

    // don't have to check the view in the hits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < hits.size(); ++h){
      art::Ptr<recob::Hit> hit = hits[h];
      std::vector<TrackIDE> hitTrackIDs = this->HitToTrackID(hit);

      // don't double count if this hit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < hitTrackIDs.size(); ++e){
	if(trackIDs.find(hitTrackIDs[e].trackID) != trackIDs.end()){
	  desired += 1.;
	  break;
	}
      }

    }// end loop over hits

    double purity = 0.;
    if(total > 0) purity = desired/total;

    return purity;
  }

  //----------------------------------------------------------------------
  double BackTracker::HitCollectionEfficiency(std::set<int>                              trackIDs, 
					      std::vector< art::Ptr<recob::Hit> > const& hits,
					      std::vector< art::Ptr<recob::Hit> > const& allhits,
					      geo::View_t                         const& view)
  {
    // get the list of EveIDs that correspond to the hits in this collection
    // and the energy associated with the desired trackID
    float desired = 0.;
    float total   = 0.;

    // don't have to check the view in the hits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < hits.size(); ++h){

      art::Ptr<recob::Hit> hit = hits[h];
      std::vector<TrackIDE> hitTrackIDs = this->HitToTrackID(hit);

      // don't worry about hits where the energy fraction for the chosen
      // trackID is < 0.1
      // also don't double count if this hit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < hitTrackIDs.size(); ++e){
	if(trackIDs.find(hitTrackIDs[e].trackID) != trackIDs.end() &&
	   hitTrackIDs[e].energyFrac             >= fMinHitEnergyFraction){
	  desired += 1.;
	  break;
	}
      }
    }// end loop over hits

    // now figure out how many hits in the whole collection are associated with this id
    for(size_t h = 0; h < allhits.size(); ++h){

      art::Ptr<recob::Hit> hit = allhits[h];

      // check that we are looking at the appropriate view here
      // in the case of 3D objects we take all hits
      if(hit->View() != view && view != geo::k3D ) continue;

      std::vector<TrackIDE> hitTrackIDs = this->HitToTrackID(hit);

      for(size_t e = 0; e < hitTrackIDs.size(); ++e){
	// don't worry about hits where the energy fraction for the chosen
	// trackID is < 0.1
	// also don't double count if this hit has more than one of the
	// desired track IDs associated with it
	if(trackIDs.find(hitTrackIDs[e].trackID) != trackIDs.end() &&
	   hitTrackIDs[e].energyFrac             >= fMinHitEnergyFraction){
	  total += 1.;
	  break;
	}
      }

    }// end loop over all hits
    
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;

    return efficiency;
  }

  //----------------------------------------------------------------------
  const sim::SimChannel* BackTracker::FindSimChannel(uint32_t channel) const
  {
    const sim::SimChannel* chan = 0;

    for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
      if(fSimChannels[sc]->Channel() == channel) chan = fSimChannels[sc];
    }

    if(!chan)
      throw cet::exception("BackTracker") << "No sim::SimChannel corresponding "
					  << "to channel: " << channel << "\n";

    return chan;
  }

  //----------------------------------------------------------------------
  void BackTracker::ChannelToTrackID(std::vector<TrackIDE>& trackIDEs,
				     uint32_t               channel,
				     double                 startTime,
				     double                 endTime)
  {
    trackIDEs.clear();

    double totalE = 0.;

    // Get services.
    art::ServiceHandle<util::DetectorProperties> detprop;

    try{
      const sim::SimChannel* schannel = this->FindSimChannel(channel);
      
      // loop over the electrons in the channel and grab those that are in time 
      // with the identified hit start and stop times

      unsigned int start_tdc = detprop->ConvertTicksToTDC(startTime);
      unsigned int end_tdc = detprop->ConvertTicksToTDC(endTime);
      std::vector<sim::IDE> simides = schannel->TrackIDsAndEnergies(start_tdc, end_tdc);
      
      // first get the total energy represented by all track ids for 
      // this channel and range of tdc values
      for(size_t e = 0; e < simides.size(); ++e)
	totalE += simides[e].energy;
      
      
      // protect against a divide by zero below
      if(totalE < 1.e-5) totalE = 1.;
      
      // loop over the entries in the map and fill the input vectors
      
      for(size_t e = 0; e < simides.size(); ++e){
	
	if(simides[e].trackID == sim::NoParticleId) continue;
	
	TrackIDE info;
	info.trackID    = simides[e].trackID;
	info.energyFrac = simides[e].energy/totalE;
	info.energy     = simides[e].energy;
	
	trackIDEs.push_back(info);
	
      }
    }// end try
    catch(cet::exception e){
      mf::LogWarning("BackTracker") << "caught exception \n"
				    << e;
    }

    return;
  }

  //----------------------------------------------------------------------
  void BackTracker::HitToSimIDEs(art::Ptr<recob::Hit> const& hit,
				 std::vector<sim::IDE>&      ides)
  {
    // Get services.
    art::ServiceHandle<util::DetectorProperties> detprop;

    unsigned int start_tdc = detprop->ConvertTicksToTDC(hit->StartTime());
    unsigned int end_tdc   = detprop->ConvertTicksToTDC(hit->EndTime());

    ides = this->FindSimChannel(hit->Channel())->TrackIDsAndEnergies(start_tdc, end_tdc);
  }

  //----------------------------------------------------------------------
  std::vector<double> BackTracker::SimIDEsToXYZ(std::vector<sim::IDE> const& ides)
  {
    std::vector<double> xyz(3, -999.);

    double x = 0.;
    double y = 0.;
    double z = 0.;
    double w = 0.;

    // loop over electrons.

    for(auto const& ide : ides) {

      double weight = ide.numElectrons;
      
      w += weight;
      x += weight * ide.x;
      y += weight * ide.y;
      z += weight * ide.z;

    }// end loop over sim::IDEs
	
    // if the sum of the weights is still 0, then return
    // the obviously stupid default values
    if(w < 1.e-5)
      throw cet::exception("BackTracker") << "No sim::IDEs providing non-zero number of electrons"
					  << " can't determine originating location from truth";

    xyz[0] = x/w;
    xyz[1] = y/w;
    xyz[2] = z/w;

    return xyz;
  }

  //----------------------------------------------------------------------
  std::vector<double> BackTracker::HitToXYZ(art::Ptr<recob::Hit> const& hit)
  {
    std::vector<sim::IDE> ides;
    HitToSimIDEs(hit, ides);
    return SimIDEsToXYZ(ides);
  }

  //----------------------------------------------------------------------
  std::vector<double> BackTracker::SpacePointToXYZ(art::Ptr<recob::SpacePoint>         const& spt,
						   art::Event                          const& evt,
						   std::string                         const& label)
  {
    // Get hits that make up this space point.
    art::PtrVector<recob::SpacePoint> spv;
    spv.push_back(spt);
    art::FindManyP<recob::Hit> fmh(spv, evt, label);
    std::vector< art::Ptr<recob::Hit> > hitv = fmh.at(0);

    // make a PtrVector
    art::PtrVector<recob::Hit> hits;
    for(size_t h = 0; h < hitv.size(); ++h) hits.push_back(hitv[h]);

    return this->SpacePointHitsToXYZ(hits);
  }

  //----------------------------------------------------------------------
  std::vector<double> BackTracker::SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits)
  {
    // Get services.
    art::ServiceHandle<geo::Geometry> geom;

    // Result vector.
    std::vector<double> xyz(3, -999.);

    // Count the hits in each plane.
    std::vector< std::vector< std::vector<int> > > numhits( geom->Ncryostats() );

    // Weight of each hit
    std::vector< std::vector< std::vector<double> > > hitWeight( geom->Ncryostats() );

    // Position of each hit   
    std::vector< std::vector< std::vector< std::vector<double> > > > hitPos( geom->Ncryostats() );
   
    for(size_t c = 0; c < numhits.size(); ++c){
      numhits[c].resize( geom->NTPC(c) );
      hitWeight[c].resize( geom->NTPC(c) );
      hitPos[c].resize( geom->NTPC(c) );
      for(size_t t = 0; t < numhits[c].size(); ++t){
	numhits[c][t].resize( geom->Nplanes(t, c) );
	hitWeight[c][t].resize( geom->Nplanes(t, c) );
	hitPos[c][t].resize( geom->Nplanes(t, c) );
      }
    }

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;

      // use the HitToXYZ and Geometry::PositionToTPC 
      // to figure out which drift volume the hit originates from
      std::vector<double> hitOrigin = this->HitToXYZ(*ihit);
      unsigned int cstat = 0;
      unsigned int tpc   = 0;
      const double worldLoc[3] = {hitOrigin[0], hitOrigin[1], hitOrigin[2]};
      geom->PositionToTPC(worldLoc, tpc, cstat);

      if(hit.WireID().Cryostat == cstat && hit.WireID().TPC == tpc){
	++numhits[cstat][tpc][hit.WireID().Plane];
	hitWeight[cstat][tpc][hit.WireID().Plane] = hit.Charge();
	hitPos[cstat][tpc][hit.WireID().Plane] = hitOrigin;
      }
	
    }

    // loop over the vectors we made and find the average position for the hits 
    // in the future we might take a weighted average
    int nhits = 0;
    xyz[0] = 0.;
    xyz[1] = 0.;
    xyz[2] = 0.;
    for(size_t c = 0; c < numhits.size(); ++c){
      for(size_t t = 0; t < numhits[c].size(); ++t){
	for(size_t p = 0; p < numhits[c][t].size(); ++p){

	  if(numhits[c][t][p] == 1) {
	    ++nhits;
	    xyz[0] += hitPos[c][t][p][0];
	    xyz[1] += hitPos[c][t][p][1];
	    xyz[2] += hitPos[c][t][p][2];
	  }

	} // end loop over planes
      } // end loop over tpcs
    } // end loop over cryostats

    // get the average position
    if(nhits < 1) 
      throw cet::exception("BackTracker") << "No hits to determine originating location from truth";

    
    xyz[0] /= nhits;
    xyz[1] /= nhits;
    xyz[2] /= nhits;

    // Done.
    return xyz;
  }

} // namespace

namespace cheat{
  DEFINE_ART_SERVICE(BackTracker)
}
