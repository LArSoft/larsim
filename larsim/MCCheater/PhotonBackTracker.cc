//
////////////////////////////////////////////////////////////////////
//
// \file PhotonBackTracker.cc
// \brief The functions needed for the PhotonBackTracker class needed by the PhotonBackTrackerService in order to connect truth information with reconstruction.
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
///////////////////////////////////////////////////////////////////
//
//TODO: Impliment alternate backtracking scheme developed by T. Usher
//TODO: OpChanToOpDetSDPs (Expanded Clone of OpDetNumToOpDetSDPs
//
///////////////////////////////////////////////////////////////////


//Includes
#include "larsim/MCCheater/PhotonBackTracker.h"

//CPP
#include <map>
//Framework

//LArSoft
#include "lardataobj/Simulation/sim.h"
#include "larsim/Simulation/SimListUtils.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

namespace cheat{

  //----------------------------------------------------------------
  PhotonBackTracker::PhotonBackTracker(fhiclConfig const& config,
      const cheat::ParticleInventory* partInv,
      const geo::GeometryCore*        geom)//,
//      const detinfo::DetectorClocks*  detClock)
    :fPartInv  (partInv),
    fGeom      (geom),
//    fDetClocks (detClock),
    fDelay     (config.Delay()),
    fG4ModuleLabel(config.G4ModuleLabel()),
    fOpHitLabel(config.OpHitLabel()),
    fOpFlashLabel(config.OpFlashLabel()),
    //fWavLabel(config.WavLabel()),
    fMinOpHitEnergyFraction(config.MinOpHitEnergyFraction())
  {}

  //----------------------------------------------------------------
  PhotonBackTracker::PhotonBackTracker( fhicl::ParameterSet const& pSet,
      const cheat::ParticleInventory* partInv,
      const geo::GeometryCore* geom)//,
//      const detinfo::DetectorClocks* detClock)
    :fPartInv (partInv),
    fGeom (geom),
//    fDetClocks(detClock),
    fDelay(pSet.get<double>("Delay")),
    fG4ModuleLabel(pSet.get<art::InputTag>("G4ModuleLabel", "largeant")),
    fOpHitLabel(pSet.get<art::InputTag>("OpHitLabel", "ophit")),
    fOpFlashLabel(pSet.get<art::InputTag>("OpFlashLabel", "opflash")),
    fMinOpHitEnergyFraction(pSet.get<double>("MinimumOpHitEnergyFraction", 0.1))
  {}


  //----------------------------------------------------------------
  const double PhotonBackTracker::GetDelay(){ return this->fDelay;}

  //----------------------------------------------------------------
  void PhotonBackTracker::ClearEvent(){
    priv_OpDetBTRs.clear();
  }

  //----------------------------------------------------------------
  const bool PhotonBackTracker::BTRsReady()
  {
    return !( priv_OpDetBTRs.empty() ) ;
  }

  //----------------------------------------------------------------
  const bool PhotonBackTracker::OpFlashToOpHitsReady()
  {
    return !( fOpFlashToOpHits.empty() ) ;
  }

  //----------------------------------------------------------------
  const std::vector< art::Ptr< sim::OpDetBacktrackerRecord >>& PhotonBackTracker::OpDetBTRs()
  {
    return priv_OpDetBTRs;
  }

  //----------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTracker::TrackIdToSimSDPs_Ps(int const& id)
  {
    std::vector< const sim::SDP* > sdp_Ps;
    for(size_t odet=0; odet<priv_OpDetBTRs.size(); ++odet){
      const auto & pdTimeSDPmap = priv_OpDetBTRs[odet]->timePDclockSDPsMap();
      for(auto mapitr = pdTimeSDPmap.begin(); mapitr != pdTimeSDPmap. end(); mapitr++){
        std::vector<sim::SDP> const& sdpvec = (*mapitr).second;
        for(size_t iv = 0; iv < sdpvec.size(); ++iv){
      //    const sim::SDP* const sdp_P = &sdpvec[iv];
          if( abs(sdpvec[iv].trackID) == id) sdp_Ps.push_back(&(sdpvec[iv]));
        }
      } // end loop over map from sim::OpDetBacktrackerRecord


    }// end loop over sim::OpDetBacktrackerRecords
    return sdp_Ps;
  }

  //----------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTracker::TrackIdToSimSDPs_Ps(int const& id, geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTracker")
      <<"PhotonBackTracker is not equiped to handle geo::Views.";
  }

  //----------------------------------------------------------------
  const art::Ptr< sim::OpDetBacktrackerRecord > PhotonBackTracker::FindOpDetBTR(int const& opDetNum) const
  {
    art::Ptr< sim::OpDetBacktrackerRecord > opDet;
    for(size_t sc = 0; sc < priv_OpDetBTRs.size(); ++sc){
      //This could become a bug. What if it occurs twice (shouldn't happen in correct records, but still, no error handeling included for the situation
      if(priv_OpDetBTRs.at(sc)->OpDetNum() == opDetNum) opDet = priv_OpDetBTRs.at(sc);
    }
    if(!opDet)
    {
      throw cet::exception("PhotonBackTracker2") << "No sim:: OpDetBacktrackerRecord corresponding "
        << "to opDetNum: " << opDetNum << "\n";
    }
    return opDet;
  }

  //----------------------------------------------------------------
  const std::vector < sim::TrackSDP > PhotonBackTracker::OpDetToTrackSDPs( int const& OpDetNum,
      double const& opHit_start_time, double const& opHit_end_time)  const
  {
    std::vector< sim::TrackSDP > tSDPs;
    double totalE=0;
    try{
      const art::Ptr< sim::OpDetBacktrackerRecord > opDetBTR =
        this->FindOpDetBTR(OpDetNum);
      // ( fGeom->OpDetFromOpChannel(channel) );
      std::vector<sim::SDP> simSDPs =
        opDetBTR->TrackIDsAndEnergies(opHit_start_time, opHit_end_time);
      for(size_t e = 0; e < simSDPs.size(); ++e)
        totalE += simSDPs[e].energy;
      if(totalE < 1.e-5) totalE = 1.;
      for(size_t e = 0; e < simSDPs.size(); ++e){
        if(simSDPs[e].trackID == sim::NoParticleId) continue;
        sim::TrackSDP info;
        info.trackID = simSDPs[e].trackID;
        info.energyFrac = simSDPs[e].energy/totalE;
        info.energy = simSDPs[e].energy;
        tSDPs.push_back(info);
      }
    }
    catch(cet::exception e){
      mf::LogWarning("PhotonBackTracker")<<"Exception caught\n"
        <<e;
    }
    return tSDPs;
  }

  //----------------------------------------------------------------
  const std::vector< sim::TrackSDP> PhotonBackTracker::OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P)  const
  {
    //auto opHit = *opHit_P;
    auto OpDetNum =  fGeom->OpDetFromOpChannel(opHit_P->OpChannel()) ;
    std::vector<sim::TrackSDP> trackSDPs;
    const double pTime = opHit_P->PeakTime();
    const double pWidth= opHit_P->Width();
    const double start = (pTime-pWidth)*1000-fDelay;
    const double end = (pTime+pWidth)*1000-fDelay;

    //this->OpDetToTrackSDPs(trackSDPs, opHit_P->OpChannel(), start, end);


    //return trackSDPs;
    //return this->OpDetToTrackSDPs( opHit_P->OpChannel(), start, end);
    return this->OpDetToTrackSDPs( OpDetNum, start, end);

  }

  //----------------------------------------------------------------
  const std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToTrackSDPs(recob::OpHit const& opHit) const
  {
    auto OpDetNum =  fGeom->OpDetFromOpChannel(opHit.OpChannel()) ;
    std::vector<sim::TrackSDP> trackSDPs;
    const double pTime = opHit.PeakTime();
    const double pWidth= opHit.Width();
    const double start = (pTime-pWidth)*1000-fDelay;
    const double end = (pTime+pWidth)*1000-fDelay;


    return this->OpDetToTrackSDPs( OpDetNum, start, end);

  }

  //----------------------------------------------------------------
  const std::vector < int > PhotonBackTracker::OpHitToTrackIds(recob::OpHit const& opHit) const
  {
    std::vector< int > retVec;
    for( auto const trackSDP : this->OpHitToTrackSDPs(opHit) ){
      retVec.push_back( trackSDP.trackID);
    }
    return retVec;
  }

  //----------------------------------------------------------------
  const std::vector < int > PhotonBackTracker::OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit) const
  {
    return this->OpHitToTrackIds(*opHit);
  }

  //----------------------------------------------------------------
  const std::vector < int > PhotonBackTracker::OpHitToEveTrackIds(recob::OpHit const& opHit)
  {/*NEW*/ /*COMPLETE*/
    std::vector< int > retVec;
    for( auto const trackSDP : this->OpHitToEveTrackSDPs(opHit) ){
      retVec.push_back( trackSDP.trackID);
    }
    return retVec;
  }

  //----------------------------------------------------------------
  const std::vector < int > PhotonBackTracker::OpHitToEveTrackIds(art::Ptr<recob::OpHit> const& opHit_P)
  {/*NEW*/ /*COMPLETE*/
    return this->OpHitToEveTrackIds(*opHit_P);
  }

  //----------------------------------------------------------------
  const std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToEveTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P)  const
  {
    return this->OpHitToEveTrackSDPs(*opHit_P);
  }

  //----------------------------------------------------------------
  const std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToEveTrackSDPs(recob::OpHit const& opHit) const
  {
    std::vector<sim::TrackSDP> trackSDPs = this->OpHitToTrackSDPs(opHit);

    // make a map of evd ID values and fraction of energy represented by
    // // that eve id in this opHit
    std::map<int, float> eveIDtoEfrac;

    double totalE = 0.;
    for(size_t t = 0; t < trackSDPs.size(); ++t){
      eveIDtoEfrac[(fPartInv->ParticleList()).EveId( trackSDPs[t].trackID )] += trackSDPs[t].energy;
      totalE += trackSDPs[t].energy;
    }

    // now fill the eveSDPs vector from the map
    std::vector<sim::TrackSDP> eveSDPs;
    eveSDPs.reserve(eveIDtoEfrac.size());
    for(auto itr = eveIDtoEfrac.begin(); itr != eveIDtoEfrac.end(); itr++){
      sim::TrackSDP temp;
      temp.trackID = (*itr).first;
      temp.energyFrac = (*itr).second/totalE;
      temp.energy = (*itr).second;
      eveSDPs.push_back(std::move(temp));
    }
    return eveSDPs;
  }

  //----------------------------------------------------------------
  //TODO: Make a copy of this function that uses an allOpHits list.
  const std::vector<art::Ptr<recob::OpHit>> PhotonBackTracker::TrackIdToOpHits_Ps( int const& tkId, std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
  {
    //One would think we would want to have this function defined, and call this function in the std::vector<tkids> to opHits, but that would require more loops (and a higher overhead.) Instead, to provide this, we will just call the existing std::vector<tkids>ToOpHits with an input of 1.
    std::vector<int> tkidFake(1, tkId);
    //std::vector<art::Ptr<recob::OpHit>> retVec = (this->TrackIdsToOpHits_Ps(tkidFake, hitsIn)).at(0);
    // return (this->TrackIdsToOpHits_Ps(tkidFake, hitsIn));
    const std::vector<art::Ptr<recob::OpHit>> out = (this->TrackIdsToOpHits_Ps(tkidFake, hitsIn)).at(0);
    return out;
  }

  //----------------------------------------------------------------
  const std::vector<std::vector<art::Ptr<recob::OpHit>>> PhotonBackTracker::TrackIdsToOpHits_Ps(std::vector<int> const& tkIds, std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
  {
    std::vector<std::pair<int, art::Ptr<recob::OpHit>>> opHitList;
    for(auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) {
      art::Ptr<recob::OpHit> const& opHit = *itr;
      auto OpDetNum = fGeom->OpDetFromOpChannel(opHit->OpChannel());
      const double pTime = opHit->PeakTime(), pWidth= opHit->Width();
      const double start = (pTime-pWidth)*1000.0-fDelay, end = (pTime+ pWidth)*1000.0-fDelay;
      std::vector<sim::TrackSDP> tids = this->OpDetToTrackSDPs( OpDetNum, start, end);
      //std::vector<sim::TrackSDP> tids = this->OpDetToTrackSDPs( opHit->OpChannel(), start, end);
      for(auto itid = tids.begin(); itid != tids.end(); ++itid) {
        for(auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
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
    for(auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
      tmpOpHits.clear();
      for(auto itr = opHitList.begin(); itr != opHitList.end(); ++itr) {
        if(*itkid == (*itr).first) tmpOpHits.push_back((*itr).second);
      }
      truOpHits.push_back(tmpOpHits);
    }
    return truOpHits;
  }

  //----------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTracker::OpHitToSimSDPs_Ps(recob::OpHit const& opHit) const
  {
    std::vector<const sim::SDP*> retVec;
    double fPeakTime = opHit.PeakTime();
    double fWidth = opHit.Width();
    sim::OpDetBacktrackerRecord::timePDclock_t start_time = ((fPeakTime- fWidth)*1000.0)-fDelay;
    sim::OpDetBacktrackerRecord::timePDclock_t end_time = ((fPeakTime+ fWidth)*1000.0)-fDelay;
    if(start_time > end_time){throw;}

    //BUG!!!fGeom->OpDetFromOpChannel(channel)
    const std::vector<std::pair<double, std::vector<sim::SDP>> >& timeSDPMap
      = (this->FindOpDetBTR(fGeom->OpDetFromOpChannel(opHit.OpChannel()) ))->timePDclockSDPsMap(); //Not guranteed to be sorted.
    //const std::vector<std::pair<double, std::vector<sim::SDP>> >& timeSDPMap = (this->FindOpDetBTR(opHit.OpChannel()))->timePDclockSDPsMap(); //Not guranteed to be sorted.

    std::vector<const std::pair<double, std::vector<sim::SDP>>*> timePDclockSDPMap_SortedPointers;
    for ( auto& pair : timeSDPMap ){ timePDclockSDPMap_SortedPointers.push_back(&pair); }
    auto pairSort = [](auto& a, auto& b) { return a->first < b->first ; } ;
    if( !std::is_sorted( timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), pairSort)){
      std::sort(timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), pairSort);
    }

    //This section is a hack to make comparisons work right.
    std::vector<sim::SDP> dummyVec;
    std::pair<double, std::vector<sim::SDP>> start_timePair = std::make_pair(start_time, dummyVec);
    std::pair<double, std::vector<sim::SDP>> end_timePair = std::make_pair(end_time, dummyVec);
    auto start_timePair_P = &start_timePair;
    auto end_timePair_P = &end_timePair;
    //First interesting iterator.
    auto mapFirst = std::lower_bound(timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), start_timePair_P, pairSort);
    //Last interesting iterator.
    auto mapLast = std::upper_bound(mapFirst, timePDclockSDPMap_SortedPointers.end(), end_timePair_P, pairSort);

    for( auto& mapitr = mapFirst; mapitr != mapLast; ++mapitr )
      for( auto& sdp : (*mapitr)->second)
        retVec.push_back(&sdp);

    return retVec;

    //sdps = FindOpDetBTR( geom->OpDetFromOpChannel(opHit. OpChannel()) )->TrackIDsAndEnergies(start_time, end_time);
    // return (this->FindOpDetBTR( fGeom->OpDetFromOpChannel(opHit.OpChannel()) ))->TrackIDsAndEnergies(start_time, end_time);

  }


  //----------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTracker::OpHitToSimSDPs_Ps(art::Ptr<recob::OpHit> const& opHit_P) const
  {
    return this->OpHitToSimSDPs_Ps(*opHit_P);
  }

  //----------------------------------------------------------------
  const std::vector< double> PhotonBackTracker::SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&
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
    //If the sum of the weights is still zero, then fail to return a value.
    //A hit with no contributing photons does't make sense.
    if(w < 1.e-5)
      throw cet::exception("PhotonBackTracker") << "No sim::SDPs providing non-zero number of photons"
        << " can't determine originating location from truth\n";
    xyz[0] = x/w;
    xyz[1] = y/w;
    xyz[2] = z/w;
    return xyz;
  }

  //----------------------------------------------------------------
  const std::vector< double> PhotonBackTracker::SimSDPsToXYZ(std::vector< const sim::SDP*> const& sdps_Ps) const&
  {
    std::vector<double> xyz(3, -999.);
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double w = 0.;
    // loop over photons.
    for(const sim::SDP* sdp_P : sdps_Ps) {
      auto& sdp = *sdp_P;
      double weight = sdp.numPhotons;
      w += weight;
      x += weight * sdp.x;
      y += weight * sdp.y;
      z += weight * sdp.z;
    }// end loop over sim::SDPs
    //If the sum of the weights is still zero, then fail to return a value.
    //A hit with no contributing photons does't make sense.
    if(w < 1.e-5)
      throw cet::exception("PhotonBackTracker") << "No sim::SDPs providing non-zero number of photons"
        << " can't determine originating location from truth\n";
    xyz[0] = x/w;
    xyz[1] = y/w;
    xyz[2] = z/w;
    return xyz;
  }

  //----------------------------------------------------------------
  const std::vector< double> PhotonBackTracker::OpHitToXYZ( recob::OpHit const& opHit)
  {
    return SimSDPsToXYZ(this->OpHitToSimSDPs_Ps(opHit));
  }

  //----------------------------------------------------------------
  const std::vector< double> PhotonBackTracker::OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit)
  {
    return SimSDPsToXYZ(this->OpHitToSimSDPs_Ps(*opHit));
  }

  //----------------------------------------------------------------
  //const std::vector< const sim::SDP* > PhotonBackTracker::OpHitToSimSDPs_Ps(recob::OpHit const& opHit)
  //  const std::vector< const sim::SDP* > PhotonBackTracker::OpHitsToSimSDPs_Ps(const std::vector< art::Ptr < recob::OpHit > >& opHits_Ps)
  const std::vector< const sim::SDP* > PhotonBackTracker::OpHitsToSimSDPs_Ps( std::vector< art::Ptr < recob::OpHit > > const& opHits_Ps) const
  {
    std::vector < const sim::SDP* > sdps_Ps;
    for ( auto opHit_P : opHits_Ps ){
      std::vector < const sim::SDP* > to_add_sdps_Ps = this->OpHitToSimSDPs_Ps(opHit_P);
      sdps_Ps.insert( sdps_Ps.end(), to_add_sdps_Ps.begin(), to_add_sdps_Ps.end() );
    }
    return sdps_Ps;
  }

  //----------------------------------------------------------------
  const std::vector< double > PhotonBackTracker::OpHitsToXYZ( std::vector < art::Ptr < recob::OpHit > > const& opHits_Ps) const
  {
    const std::vector<const sim::SDP*> SDPs_Ps = OpHitsToSimSDPs_Ps(opHits_Ps);
    return this->SimSDPsToXYZ(SDPs_Ps);
  }

  //----------------------------------------------------------------
  const std::unordered_set<const sim::SDP*> PhotonBackTracker::OpHitToEveSimSDPs_Ps(recob::OpHit const& opHit_P)
  { /*NEW*/ /*COMPLETE*/
    const std::vector < int > ids = this->OpHitToEveTrackIds(opHit_P);
    std::unordered_set <const sim::SDP* > sdps;
    for( auto const& id : ids ){
      std::vector<const sim::SDP* > tmp_sdps = TrackIdToSimSDPs_Ps(id);
      for( const sim::SDP* tmp_sdp : tmp_sdps ){
        sdps.insert(tmp_sdp); //emplace not needed here.
      }
    }
    return sdps;
  }

  //----------------------------------------------------------------
  const std::unordered_set<const sim::SDP*> PhotonBackTracker::OpHitToEveSimSDPs_Ps(art::Ptr<recob::OpHit>&  opHit)
  { /*NEW*/ /*COMPLETE*/
    const std::vector < int > ids = this->OpHitToEveTrackIds(opHit);
    std::unordered_set <const sim::SDP* > sdps;
    for( auto const& id : ids ){
      std::vector<const sim::SDP* > tmp_sdps = TrackIdToSimSDPs_Ps(id);
      for( const sim::SDP* tmp_sdp : tmp_sdps ){
        sdps.insert(tmp_sdp); //emplace not needed here.
      }
    }
    return sdps;
  }

  //----------------------------------------------------------------
  const std::set<int> PhotonBackTracker::GetSetOfEveIds() const
  {
    //std::set<int> out = fPartInv->GetSetOfEveIds();
    return fPartInv->GetSetOfEveIds();
    //return out;
  }

  //----------------------------------------------------------------
  const std::set<int> PhotonBackTracker::GetSetOfTrackIds() const
  {
    return fPartInv->GetSetOfTrackIds();
  }

  //----------------------------------------------------------------
  const std::set<int> PhotonBackTracker::GetSetOfEveIds(std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps) const
  {
    std::set<int> eveIds;
    for(auto const& opHit_P : opHits_Ps){
      const std::vector<sim::TrackSDP> sdps = this->OpHitToEveTrackSDPs(opHit_P);
      for(auto const& sdp : sdps){eveIds.insert(sdp.trackID);}//end sdps
    }//End for hits
    return eveIds;
  }

  //----------------------------------------------------------------
  const std::set<int> PhotonBackTracker::GetSetOfEveIds(std::vector< recob::OpHit> const& opHits) const
  { /*NEW*/ /*COMPLETE*/
    std::set<int> eveIds;
    for(auto const& opHit : opHits){
      const std::vector<sim::TrackSDP> sdps = this->OpHitToEveTrackSDPs(opHit);
      for(auto const& sdp : sdps){eveIds.insert(sdp.trackID);}//end sdps
    }//End for hits
    return eveIds;
  }

  //----------------------------------------------------------------
  const std::set<int> PhotonBackTracker::GetSetOfTrackIds(std::vector< art::Ptr<recob::OpHit> > const& opHits) const
  {
    std::set<int> tids;
    for( auto const& opHit : opHits){
      for(auto const& sdp : this->OpHitToTrackSDPs(opHit)) {
        tids.insert(sdp.trackID);
      }//End for TrackSDPs
    }//End for hits
    return tids;
  }

  //----------------------------------------------------------------
  const std::set<int> PhotonBackTracker::GetSetOfTrackIds(std::vector< recob::OpHit > const& opHits) const
  { /*NEW*/ /*COMPLETE*/
    std::set<int> tids;
    for( auto const& opHit : opHits){
      for(auto const& sdp : this->OpHitToTrackSDPs(opHit)) {
        tids.insert(sdp.trackID);
      }//End for TrackSDPs
    }//End for hits
    return tids;
  }

  //----------------------------------------------------------------
  const double PhotonBackTracker::OpHitCollectionPurity(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection
    // if the EveID shows up in the input list of tkIds, then it counts
    float total = 1.*opHits.size();;
    float desired = 0.;
    for(size_t h = 0; h < opHits.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackSDPs = this->OpHitToTrackSDPs(opHit);
      // don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackSDPs.size(); ++e){
        if(tkIds.find(opHitTrackSDPs[e].trackID) != tkIds.end()){
          desired += 1.;
          break;
        }
      }
    }// end loop over opHits
    double purity = 0.;
    if(total > 0) purity = desired/total;
    return purity;
  }

  //----------------------------------------------------------------
  const double PhotonBackTracker::OpHitLightCollectionPurity(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection
    // if the EveID shows up in the input list of tkIds, then it counts
    float total = 0;
    float desired = 0.;
    // don't have to check the view in the opHits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < opHits.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackSDPs(opHit);
      total+=opHit->Area(); // sum up the charge in the cluster
      // don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
        if(tkIds.find(opHitTrackIDs[e].trackID) != tkIds.end()){
          desired += opHit->Area();
          break;
        }
      }
    }// end loop over opHits
    double purity = 0.;
    if(total > 0) purity = desired/total;
    return purity;
  }

  //----------------------------------------------------------------
  const double PhotonBackTracker::OpHitCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn,
      geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTracker")<<"This function is not supported. OpHits do not have type View.\n";
  }

  //----------------------------------------------------------------
  const double PhotonBackTracker::OpHitCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn)
  {
    float desired = 0.;
    float total = 0.;
    for(size_t h = 0; h < opHits.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackSDPs = this->OpHitToTrackSDPs(opHit);
      // also don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackSDPs.size(); ++e){
        if(tkIds.find(opHitTrackSDPs[e].trackID) != tkIds.end() &&
            opHitTrackSDPs[e].energyFrac >= fMinOpHitEnergyFraction){
          desired += 1.;
          break;
        }
      }
    }// end loop over opHits
    // now figure out how many opHits in the whole collection are associated with this id
    for(size_t h = 0; h < opHitsIn.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHitsIn[h];
      std::vector<sim::TrackSDP> opHitTrackSDPs = this->OpHitToTrackSDPs(opHit);
      for(size_t e = 0; e < opHitTrackSDPs.size(); ++e){
        // don't worry about opHits where the energy fraction for the chosen
        // trackID is < 0.1
        // also don't double count if this opHit has more than one of the
        // desired track IDs associated with it
        if(tkIds.find(opHitTrackSDPs[e].trackID) != tkIds.end() &&
            opHitTrackSDPs[e].energyFrac >= fMinOpHitEnergyFraction){
          total += 1.;
          break;
        }
      }
    }// end loop over all opHits
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;
    return efficiency;
  }

  //----------------------------------------------------------------
  const double PhotonBackTracker::OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn,
      geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTracker")<<"This function is not supported. OpHits do not have type View.\n";
  }

  //----------------------------------------------------------------
  const double PhotonBackTracker::OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn)
  {
    float desired = 0.;
    float total = 0.;

    // don't have to check the view in the opHits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(size_t h = 0; h < opHits.size(); ++h){

      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackSDPs(opHit);

      // don't worry about opHits where the energy fraction for the chosen
      // trackID is < 0.1
      // also don't double count if this opHit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
        if(tkIds.find(opHitTrackIDs[e].trackID) != tkIds.end() &&
            opHitTrackIDs[e].energyFrac >= fMinOpHitEnergyFraction){
          desired += opHit->Area();
          break;
        }
      }
    }// end loop over opHits
    for(size_t h = 0; h < opHitsIn.size(); ++h){
      art::Ptr<recob::OpHit> opHit = opHitsIn[h];
      // check that we are looking at the appropriate view here
      // in the case of 3D objects we take all opHits
      //if(opHit->View() != view && view != geo::k3D ) continue;
      std::vector<sim::TrackSDP> opHitTrackIDs = this->OpHitToTrackSDPs(opHit);
      for(size_t e = 0; e < opHitTrackIDs.size(); ++e){
        // don't worry about opHits where the energy fraction for the chosen
        // trackID is < 0.1
        // also don't double count if this opHit has more than one of the
        // desired track IDs associated with it
        if(tkIds.find(opHitTrackIDs[e].trackID) != tkIds.end() &&
            opHitTrackIDs[e].energyFrac >= fMinOpHitEnergyFraction){
          total += opHit->Area();
          break;
        }
      }
    }// end loop over all opHits
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;
    return efficiency;
  }
  //--------------------------------------------------
//  const std::vector< const recob::OpHit* > PhotonBackTracker::OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P) const
    //const std::vector<art::Ptr<recob::OpHit>> PhotonBackTracker::OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P, Evt const& evt) const
//  {//There is not "non-pointer" version of this because the art::Ptr is needed to look up the assn. One could loop the Ptrs and dereference them, but I will not encourage the behavior by building the tool to do it.
//
    //      art::FindManyP< recob::OpHit > fmoh(std::vector<art::Ptr<recob::OpFlash>>({flash_P}), evt, fOpHitLabel.label());
    //      std::vector<art::Ptr<recob::OpHit>> const& hits_Ps = fmoh.at(0);
//    std::vector<const recob::OpHit*> const& hits_Ps = fOpFlashToOpHits.at(flash_P);
//    std::vector<art::Ptr<recob::OpHit>> const& hits_Ps = fOpFlashToOpHits.at(flash_P);
//    return hits_Ps;
//
//  }

  //--------------------------------------------------
//  const std::vector<double> PhotonBackTracker::OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P) const
//  {
//    const std::vector< const recob::OpHit *> opHits_Ps = this->OpFlashToOpHits_Ps(flash_P);
//    const std::vector<double> retVec = this->OpHitsToXYZ(opHits_Ps);
//    const std::vector<double> retVec(0.0,3);
//    //This feature temporarily disabled.
//    return retVec;
//  }

  //--------------------------------------------------
//  const std::set<int> PhotonBackTracker::OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P) const{
    /* Temporarily disabled.
    std::vector<art::Ptr<recob::OpHit> > opHits_Ps = this->OpFlashToOpHits_Ps(flash_P);
    std::set<int> ids;
    for( auto& opHit_P : opHits_Ps){
      for( const int& id : this->OpHitToTrackIds(opHit_P) ){
        ids.insert( id) ;
      } // end for ids
    }// end for opHits
    */
//    std::set<int> ids;
//    return ids;
//  }// end OpFlashToTrackIds

  //----------------------------------------------------- /*NEW*/
  //----------------------------------------------------- /*NEW*/
  //----------------------------------------------------- /*NEW*/
  //----------------------------------------------------- /*NEW*/
  //----------------------------------------------------- /*NEW*/
  //std::vector<sim::TrackSDP> OpFlashToTrackSDPs(art::Ptr<recob::OpFlash> flash_P);
  //----------------------------------------------------- /*NEW*/
  //std::vector<sim::TrackSDP> OpFlashToTrackSDPs(recob::OpFlash flash);
  //----------------------------------------------------- /*NEW*/
  //std::vector<sim::TrackSDP> OpFlashToEveTrackSDPs(recob::OpFlash flash);
  //----------------------------------------------------- /*NEW*/
  //std::vector<sim::TrackSDP> OpFlashToEveTrackSDPs(art::Ptr<recob::OpFlash> flash_P);
  //----------------------------------------------------- /*NEW*/
  //std::vector<sim::SDP*> OpFlashToSimSDPs_Ps(recob::OpFlash flash);
  //----------------------------------------------------- /*NEW*/
  //std::vector<sim::SDP*> OpFlashToSimSDPs_Ps(art::Ptr<recob::OpFlash> flash_P);




  //----------------------------------------------------------------------
} // namespace

