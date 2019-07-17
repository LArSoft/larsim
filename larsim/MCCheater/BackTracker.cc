////////////////////////////////////////////////////////////////////
//
//
// \file BackTracker.cc
// \brief The functions needed for the BackTracker class needed by the BackTracker service in order to connect truth information with reconstruction.
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
//
///////////////////////////////////////////////////////////////////

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "BackTracker.h"
#include "lardataobj/Simulation/sim.h"

#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"


namespace cheat{

  //-----------------------------------------------------------------------
  BackTracker::BackTracker(const fhiclConfig& config, const cheat::ParticleInventory* partInv, const geo::GeometryCore* geom, const detinfo::DetectorClocks* detClock )
    :fPartInv  (partInv),
     fGeom     (geom),
     fDetClocks(detClock),
     fG4ModuleLabel(config.G4ModuleLabel()),
     fHitLabel(config.DefaultHitModuleLabel()),
     fMinHitEnergyFraction(config.MinHitEnergyFraction()),
     fOverrideRealData(config.OverrideRealData())
  {
  }

  //-----------------------------------------------------------------------
  BackTracker::BackTracker(const fhicl::ParameterSet& pSet, const cheat::ParticleInventory* partInv, const geo::GeometryCore* geom, const detinfo::DetectorClocks* detClock)
    :fPartInv  (partInv),
     fGeom     (geom),
     fDetClocks(detClock),
     fG4ModuleLabel       (pSet.get<art::InputTag>("G4ModuleLabel", "largeant")),
     fHitLabel            (pSet.get<art::InputTag>("DefaultHitModuleLabel", "hitfd")),
     fMinHitEnergyFraction(pSet.get<double>       ("MinHitEnergyFraction", 0.010)),
     fOverrideRealData    (pSet.get<bool>         ("OverrideRealData",false))
  {
  }

  //-----------------------------------------------------------------------
  void BackTracker::ClearEvent(){
    fSimChannels.clear();
//    fAllHitList.clear();
  }

  //-----------------------------------------------------------------------
  std::vector< const sim::IDE* > BackTracker::TrackIdToSimIDEs_Ps( int const& id ) const{
    std::vector< const sim::IDE* > ideps;
    for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
      const auto & tdcidemap = fSimChannels[sc]->TDCIDEMap(); //This returns a reference.
      // loop over the IDEMAP
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
        // loop over the vector of IDE objects.
        const std::vector<sim::IDE>& idevec = (*mapitr).second; //note, mapitr.second returns the actual data from the map, not a copy
        for(size_t iv = 0; iv < idevec.size(); ++iv){
          //const sim::IDE* const idep = &idevec[iv];
          //if( abs(idevec[iv].trackID) == id) continue;
          //ideps.push_back(idep);
          if( abs(idevec[iv].trackID) == id) ideps.push_back(&(idevec[iv]));
        }//end for index in idevec
      } // end loop over map from sim::SimChannel
    } // end loop over sim::SimChannels
    return ideps;
  }


  //-----------------------------------------------------------------------
  std::vector<const sim::IDE* >   BackTracker::TrackIdToSimIDEs_Ps (int const& id, const geo::View_t view) const
  {
    std::vector<const sim::IDE*> ide_Ps;
    for(const art::Ptr<sim::SimChannel> sc : fSimChannels){
      if (fGeom->View(sc->Channel()) != view) continue;

      // loop over the IDEMAP
      for(const auto & item : sc->TDCIDEMap()){

        // loop over the vector of IDE objects.
        for(const sim::IDE & ide : item.second){
          if (abs(ide.trackID) == id) ide_Ps.push_back(&ide);
        }
      } // end loop over map from sim::SimChannel
    } // end loop over sim::SimChannels

    return ide_Ps;
  }


  //-----------------------------------------------------------------------
  art::Ptr<sim::SimChannel> BackTracker::FindSimChannel(raw::ChannelID_t channel) const{
    art::Ptr<sim::SimChannel> chan;
    auto ilb = std::lower_bound(fSimChannels.begin(),fSimChannels.end(),channel,[](art::Ptr<sim::SimChannel> a, raw::    ChannelID_t channel) {return(a->Channel()<channel);});
    if (ilb != fSimChannels.end())
      if ( (*ilb)->Channel() == channel) {chan = *ilb;}
    if(!chan)
      throw cet::exception("BackTracker") << "No sim::SimChannel corresponding "
        << "to channel: " << channel << "\n";
    return chan;
  }


  //-----------------------------------------------------------------------
  std::vector< sim::TrackIDE > BackTracker::ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time, const double hit_end_time) const{
    std::vector< sim::TrackIDE > trackIDEs;
    double totalE=0.;
    try{
      art::Ptr<sim::SimChannel> schannel = this->FindSimChannel(channel);

      // loop over the electrons in the channel and grab those that are in time
      // with the identified hit start and stop times
      int start_tdc = fDetClocks->TPCTick2TDC( hit_start_time );
      int end_tdc   = fDetClocks->TPCTick2TDC( hit_end_time   );
      if(start_tdc<0) start_tdc = 0;
      if(end_tdc<0) end_tdc = 0;
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

        sim::TrackIDE info;
        info.trackID    = simides[e].trackID;
        info.energyFrac = simides[e].energy/totalE;
        info.energy     = simides[e].energy;
        info.numElectrons = simides[e].numElectrons;

        trackIDEs.push_back(info);

      }
    }// end try
    catch(cet::exception e){
      mf::LogWarning("BackTracker") << "caught exception \n"
        << e;
    }

    return trackIDEs;

  }


  //-----------------------------------------------------------------------
  std::vector< sim::TrackIDE > BackTracker::HitToTrackIDEs( recob::Hit const& hit) const {
    std::vector<  sim::TrackIDE > trackIDEs;
    const double start = hit.PeakTimeMinusRMS();
    const double end   = hit.PeakTimePlusRMS();
    trackIDEs = this->ChannelToTrackIDEs(hit.Channel(), start, end);
    return trackIDEs;
  }


  //-----------------------------------------------------------------------
  std::vector< int > BackTracker::HitToTrackIds(recob::Hit const& hit) const {
    std::vector< int > retVec;
    for(auto const trackIDE : this->HitToTrackIDEs( hit ) ){
      retVec.push_back( trackIDE.trackID );
    }
    return retVec;
  }

  //These don't exist in the event. They are generated on the fly.
  //----------------------------------------------------------------------------
  std::vector<sim::TrackIDE> BackTracker::HitToEveTrackIDEs(recob::Hit const& hit) const{
    std::vector<sim::TrackIDE> eveIDEs;
    std::vector<sim::TrackIDE> trackIDEs = this->HitToTrackIDEs(hit);
    std::map<int, std::pair<double,double>> eveToEMap;
    double totalE=0.0;
    for (const auto& trackIDE : trackIDEs){
      auto check = eveToEMap.emplace(
          fPartInv->TrackIdToEveTrackId(trackIDE.trackID),
          std::make_pair(trackIDE.energy,trackIDE.numElectrons));
      if(check.second==false)
      {
        check.first->second.first +=trackIDE.energy;
        check.first->second.second +=trackIDE.numElectrons;
      }
      //      eveToEMap[fPartInv->TrackIdToEveTrackId(trackIDE.trackID)].first
      //      += trackIDE.energy;
      totalE+=trackIDE.energy;
    }//End for trackIDEs
    eveIDEs.reserve(eveToEMap.size());
    for(const auto& eveToE : eveToEMap){
      sim::TrackIDE eveTrackIDE_tmp;

      eveTrackIDE_tmp.trackID       = eveToE.first;
      eveTrackIDE_tmp.energy        = eveToE.second.first;
      eveTrackIDE_tmp.energyFrac    = (eveTrackIDE_tmp.energy)/(totalE);
      eveTrackIDE_tmp.numElectrons  = eveToE.second.second;

      eveIDEs.push_back(eveTrackIDE_tmp);
    }//END eveToEMap loop
    return eveIDEs;
  }

  //-----------------------------------------------------------------------
  std::vector < art::Ptr< recob::Hit > > BackTracker::TrackIdToHits_Ps( const int& tkId, std::vector<art::Ptr<recob::Hit>> const& hitsIn ) const{
    // returns a subset of the hits in the hitsIn collection that are matched
    // to the given track

    // temporary vector of TrackIds and Ptrs to hits so only one
    // loop through the (possibly large) hitsIn collection is needed
    std::vector< art::Ptr<recob::Hit>> hitList;
    std::vector<sim::TrackIDE> trackIDE;
    for(auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) {
      trackIDE.clear();
      art::Ptr<recob::Hit> const& hit = *itr;
      trackIDE = this->ChannelToTrackIDEs( hit->Channel(), hit->PeakTimeMinusRMS(), hit->PeakTimePlusRMS());
      for(auto itr_trakIDE = trackIDE.begin(); itr_trakIDE != trackIDE.end(); ++itr_trakIDE) {
        if(itr_trakIDE->trackID == tkId && itr_trakIDE->energyFrac > fMinHitEnergyFraction)
          hitList.push_back( hit);
      } // itr_trakIDE
    } // itr
    return hitList;
  }


  //-----------------------------------------------------------------------
  //This function could clearly be made by calling TrackIdToHits for each trackId, but that would be significantly slower because we would loop through all hits many times.
  std::vector< std::vector< art::Ptr<recob::Hit> > > BackTracker::TrackIdsToHits_Ps( std::vector<int> const& tkIds, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const{
    // returns a subset of the hits in the hitsIn collection that are matched
    // to MC particles listed in tkIds

    // temporary vector of TrackIds and Ptrs to hits so only one
    // loop through the (possibly large) hitsIn collection is needed
    std::vector<std::pair<int, art::Ptr<recob::Hit>>> hitList;
    std::vector<sim::TrackIDE> tids;
    for(auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) {
      tids.clear();
      art::Ptr<recob::Hit> const& hit = *itr;
      tids=this->ChannelToTrackIDEs( hit->Channel(), hit->PeakTimeMinusRMS(), hit->PeakTimePlusRMS());
      for(auto itid = tids.begin(); itid != tids.end(); ++itid) {
        for(auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
          if(itid->trackID == *itkid) {
            if(itid->energyFrac > fMinHitEnergyFraction)
              hitList.push_back(std::make_pair(*itkid, hit));
          }
        } // itkid
      } // itid
    } // itr

    // now build the truHits vector that will be returned to the caller
    std::vector<std::vector<art::Ptr<recob::Hit>>> truHits;
    // temporary vector containing hits assigned to one MC particle
    std::vector<art::Ptr<recob::Hit>> tmpHits;
    for(auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
      tmpHits.clear();
      for(auto itr = hitList.begin(); itr != hitList.end(); ++itr) {
        if(*itkid == (*itr).first) tmpHits.push_back((*itr).second);
      }
      truHits.push_back(tmpHits);
    }
    return truHits;
  }

  //-----------------------------------------------------------------------
  //Cannot be returned as a reference, as these IDEs do not exist in the event. They are constructed on the fly.
  std::vector<sim::IDE> BackTracker::HitToAvgSimIDEs (recob::Hit const& hit) const{
    // Get services.

    int start_tdc = fDetClocks->TPCTick2TDC( hit.PeakTimeMinusRMS() );
    int end_tdc   = fDetClocks->TPCTick2TDC( hit.PeakTimePlusRMS()   );
    if(start_tdc<0) start_tdc = 0;
    if(end_tdc<0) end_tdc = 0;

    return (this->FindSimChannel(hit.Channel()))->TrackIDsAndEnergies(start_tdc, end_tdc);
  }

  //-----------------------------------------------------------------------
  std::vector< const sim::IDE* > BackTracker::HitToSimIDEs_Ps (recob::Hit const& hit) const{
    std::vector< const sim::IDE* > retVec;
    //    const auto start_tdc = hit.PeakTimeMinusRMS();
    //    const auto end_tdc = hit.PeakTimePlusRMS();
    int start_tdc = fDetClocks->TPCTick2TDC( hit.PeakTimeMinusRMS() );
    int end_tdc   = fDetClocks->TPCTick2TDC( hit.PeakTimePlusRMS()   );
    if(start_tdc<0) start_tdc = 0;
    if(end_tdc<0) end_tdc = 0;

    if(start_tdc > end_tdc){throw;}

    /*    auto sc = this->FindSimChannel(hit.Channel());
          auto tdcidemap = sc->TDCIDEMap();
          auto& tdcidemapref = sc->TDCIDEMap();*/

    const std::vector< std::pair<unsigned short, std::vector<sim::IDE>> >& tdcIDEMap = (this->FindSimChannel(hit.Channel()))->TDCIDEMap(); //This in fact does not return a map. It returns a vector... with no guarantee that it is sorted...
    std::vector<const std::pair<unsigned short, std::vector<sim::IDE>>*> tdcIDEMap_SortedPointers;
    for ( auto& pair : tdcIDEMap ){
      tdcIDEMap_SortedPointers.push_back(&pair);
    }

    //Sort a vector of pointers to the TDCIDEs
    /*
       auto pairSort = [](auto& a, auto& b) { return a.first < b.first ; } ;
       if( !std::is_sorted( tdcIDEMap.begin(), tdcIDEMap.end(), pairSort)) {
       std::cout<<"AlreadySorted\n";
    //   std::sort (tdcIDEMap.begin(), tdcIDEMap.end(), pairSort);
    }else{
    throw;
    }*/

    //This is a bunch of extra steps, due to needing a vector we can sort, and needing those items in the sorted vector to be the items from the sim channels (so a pointer to the IDEs inside the sim channels can be made). The work around is to make a vector of pointers to IDEs inside the TDCIDEMap (which is a constant reference to the fTDCIDEs in the SimChannel.)
    auto pairSort = [](auto& a, auto& b) { return a->first < b->first ; } ;
    if( !std::is_sorted( tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), pairSort)) {
      std::sort (tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), pairSort);
    }


    std::vector<sim::IDE> dummyVec; //I need something to stick in a pair to compare pair<tdcVal, IDE>. This is an otherwise useless "hack".
    std::pair<double, std::vector<sim::IDE>> start_tdcPair = std::make_pair(start_tdc,dummyVec); //This pair is a "hack" to make my comparison work for lower and upper bound.
    std::pair<double, std::vector<sim::IDE>> end_tdcPair = std::make_pair(end_tdc,dummyVec);
    //const std::vector<std::pair<unsigned short, std::vector<sim::IDE> > >::iterator  //iterator to the first interesting IDE
    auto start_tdcPair_P = &start_tdcPair;
    auto end_tdcPair_P = &end_tdcPair;
    auto
      mapFirst = std::lower_bound(tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), start_tdcPair_P, pairSort);
    //mapFirst = std::lower_bound(tdcIDEMap.begin(), tdcIDEMap.end(), start_tdcPair, pairSort);
    //const std::vector<std::pair<unsigned short, std::vector<sim::IDE> > >::iterator  //iterator to just after the last interesting IDE
    auto
      mapLast  = std::upper_bound(tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), end_tdcPair_P, pairSort);
    for( auto& mapitr = mapFirst; mapitr != mapLast; ++mapitr ){
      for( auto& ide : (*mapitr)->second){
        retVec.push_back(&ide);
        //std::cout<<"Dumping Selected IDEs from tdcIDEMap:\nTrackID: "<<ide.trackID<<"\nnumElectrons: "<<ide.numElectrons<<"\nenergy: "<<ide.energy<<"\nXYZ: "<<ide.x<<" "<<ide.y<<" "<<ide.z<<"\n";
      } //Add all interesting IDEs to the retVec
    }
    /*    for( auto ide : retVec){
          std::cout<<"Dumping Selected IDEs from retVec:\nTrackID: "<<ide->trackID<<"\nnumElectrons: "<<ide->numElectrons<<"\nenergy: "<<ide->energy<<"\nXYZ: "<<ide->x<<" "<<ide->y<<" "<<ide->z<<"\n";
          }*/
    return retVec;
  }

  //------------------------------------------------------------------------------
  std::vector<double> BackTracker::SimIDEsToXYZ(std::vector<sim::IDE> const& ides) const
  {
    std::vector<double> xyz(3,0.0);
    double w = 0.0;
    for( auto const& ide : ides){
      double weight=ide.numElectrons;
      w      += weight;
      xyz[0] += (weight*ide.x);
      xyz[1] += (weight*ide.y);
      xyz[2] += (weight*ide.z);
    }
    if(w<1.e-5)
      throw cet::exception("BackTracker")<<"No sim::IDEs providing non-zero number of electrons"
        << " can't determine originating location from truth\n";
    xyz[0] = xyz[0]/w;
    xyz[1] = xyz[1]/w;
    xyz[2] = xyz[2]/w;
    return xyz;
  }

  //-------------------------------------------------------------------------------
  std::vector<double> BackTracker::SimIDEsToXYZ( std::vector< const sim::IDE* > const& ide_Ps) const{
    std::vector<sim::IDE> ides;
    for(auto ide_P : ide_Ps ){ ides.push_back(*ide_P);}
    return this->SimIDEsToXYZ(ides);
  }


  //--------------------------------------------------------------------------------
  std::vector<double> BackTracker::HitToXYZ(const recob::Hit& hit) const{
    std::vector<const sim::IDE*> ide_Ps = this->HitToSimIDEs_Ps(hit);
    return this->SimIDEsToXYZ(ide_Ps);
  }

  //-----------------------------------------------------------------------------------
  double BackTracker::HitCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits) const {
    int desired =0;
    for( const auto& hit : hits ){
      std::vector<sim::TrackIDE> hitTrackIDEs=this->HitToTrackIDEs(hit);
      for(const auto& tIDE : hitTrackIDEs){
        if(trackIds.find(tIDE.trackID)!=trackIds.end()){
          ++desired;
          break;
        }//End if TID Found
      }//END for trackIDE in TrackIDEs
    }//End for hit in hits
    if(hits.size()>0){return double(double(desired)/double(hits.size()));}
    return 0;
  }

  //-----------------------------------------------------------------------------------
  double BackTracker::HitChargeCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const&     hits) const {
    double totalCharge=0.,desired=0.;
    for(const auto& hit : hits){
      totalCharge+=hit->Integral();
      std::vector<sim::TrackIDE> trackIDEs= this->HitToTrackIDEs(hit);
      for(const auto& trackIDE : trackIDEs){
        if(trackIds.find(trackIDE.trackID)!=trackIds.end()){
          desired+=hit->Integral();
          break;
        }//End if trackId in trackIds.
      }//End for trackIDE in trackIDEs
    }//End for Hit in Hits
    if(totalCharge>0.0){return (desired/totalCharge);}
    return 0.0;
  }

  //-----------------------------------------------------------------------------------
  double BackTracker::HitCollectionEfficiency( std::set<int> const& trackIds,
      std::vector< art::Ptr<recob::Hit> > const& hits,
      std::vector< art::Ptr<recob::Hit> > const& allHits,
      geo::View_t const& view) const {

    int desired=0,total=0;

    for( const auto& hit : hits){
      std::vector<sim::TrackIDE> hitTrackIDEs = this->HitToTrackIDEs(hit);
      for( const auto& trackIDE : hitTrackIDEs){
        if( trackIds.find(trackIDE.trackID)!=trackIds.end() && trackIDE.energyFrac >= fMinHitEnergyFraction){
          ++desired;
          break;
        }//End if trackID in trackIds.
      }//end for trackIDE in TrackIDEs
    }//end for hit in hits

    for( const auto& hit : allHits){
      if(hit->View()!=view && view != geo::k3D) {continue;}//End if hit.view = view or view = geo::k3D
      std::vector<sim::TrackIDE> hitTrackIDEs = this->HitToTrackIDEs(hit);
      for( const auto& hitIDE : hitTrackIDEs){
        if(trackIds.find(hitIDE.trackID)!=trackIds.end() && hitIDE.energyFrac>=fMinHitEnergyFraction){
          ++total;
          break;
        }
      }//END for all IDEs in HitTrackIDEs.
    }//end for hit in allHits.
    if(total >= 0){ return double(double(desired)/double(total));}
    return 0.;
  }

  //-----------------------------------------------------------------------------------
  double BackTracker::HitChargeCollectionEfficiency(std::set<int> trackIds,
      std::vector< art::Ptr<recob::Hit> > const& hits,
      std::vector< art::Ptr<recob::Hit> > const& allHits,
      geo::View_t const& view) const{
    double desired=0.,total=0.;
    for( const auto& hit : hits){
      std::vector<sim::TrackIDE> hitTrackIDEs = this->HitToTrackIDEs(hit);
      for(const auto& hitIDE : hitTrackIDEs){
        if(trackIds.find(hitIDE.trackID) != trackIds.end() && hitIDE.energyFrac >= fMinHitEnergyFraction){
          desired+=hit->Integral();
          break;
        }//end if hit id matches and energy sufficient.
      }//End for IDE in HitTrackIDEs.
    }//End for hit in hits.

    for( const auto& hit : allHits ){
      if(hit->View() != view && view != geo::k3D){ continue; }
      std::vector<sim::TrackIDE> hitTrackIDEs = this->HitToTrackIDEs(hit);
      for( const auto& hitIDE : hitTrackIDEs ){
        if(trackIds.find(hitIDE.trackID) != trackIds.end() && hitIDE.energyFrac >= fMinHitEnergyFraction){
          total += hit->Integral();
          break;
        }//end if hit matches
      }//end for ide in ides
    }//End for hit in allHits

    if(total>0.) {return desired/total;}
    return 0.;

  }


  //-----------------------------------------------------------------------------------
  std::set<int> BackTracker::GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const{
    std::set<int> tids;
    for( const auto& hit : hits){
      const double start = hit->PeakTimeMinusRMS();
      const double end   = hit->PeakTimePlusRMS();
      std::vector<sim::TrackIDE> trackIDEs = this->ChannelToTrackIDEs(hit->Channel(), start, end);
      for(const auto& ide : trackIDEs) {
        tids.insert(ide.trackID);
      }//End for TrackIDEs
    }//End for hits
    return tids;
  }//End GetSetOfTrackIds

  //-----------------------------------------------------------------------------------
  std::set<int> BackTracker::GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const {
    std::set<int> eveIds;
    for(const auto& hit : hits){
      const std::vector<sim::TrackIDE> ides = this->HitToEveTrackIDEs(hit);
      for(const auto& ide : ides){eveIds.insert(ide.trackID);}//end ides
    }//End for hits
    return eveIds;
  }


  //This function definitely needs a new implimentation. There must be abetter way than so many loops.
  std::vector< double> BackTracker::SpacePointHitsToWeightedXYZ(std::vector<art::Ptr<recob::Hit>> const& hits) const {
    std::vector<double> xyz(3,-99999.9);
    std::vector< std::vector<std::vector<int>>> numHits (fGeom->Ncryostats());
    std::vector< std::vector<std::vector<double>>> hitWeight (fGeom->Ncryostats());
    std::vector< std::vector<std::vector<std::vector<double>>>> hitPos (fGeom->Ncryostats());
    //Do we need to resize everything...
    for(size_t c = 0; c < numHits.size(); ++c){
      numHits[c].resize( fGeom->NTPC(c) );
      hitWeight[c].resize( fGeom->NTPC(c) );
      hitPos[c].resize( fGeom->NTPC(c) );
      for(size_t t = 0; t < numHits[c].size(); ++t){
        numHits[c][t].resize( fGeom->Nplanes(t, c) );
        hitWeight[c][t].resize( fGeom->Nplanes(t, c) );
        hitPos[c][t].resize( fGeom->Nplanes(t, c) );
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
      fGeom->PositionToTPC(worldLoc, tpc, cstat);

      if(hit.WireID().Cryostat == cstat && hit.WireID().TPC == tpc){
        ++numHits[cstat][tpc][hit.WireID().Plane];
        hitWeight[cstat][tpc][hit.WireID().Plane] = hit.Integral();
        hitPos[cstat][tpc][hit.WireID().Plane] = hitOrigin;
      }

    }

    // loop over the vectors we made and find the average position for the hits
    // in the future we might take a weighted average
    int nhits = 0;
    xyz[0] = 0.;
    xyz[1] = 0.;
    xyz[2] = 0.;
    for(size_t c = 0; c < numHits.size(); ++c){
      for(size_t t = 0; t < numHits[c].size(); ++t){
        for(size_t p = 0; p < numHits[c][t].size(); ++p){

          if(numHits[c][t][p] == 1) {
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
      throw cet::exception("BackTracker") << "No hits to determine originating location from truth\n";


    xyz[0] /= nhits;
    xyz[1] /= nhits;
    xyz[2] /= nhits;

    // Done.
    return xyz;

  }


}//End namespace cheat
