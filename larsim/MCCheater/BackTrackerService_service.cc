////////////////////////////////////////////////////////////////////////////////////////
//
// \file BackTrackerService_service.cc
// \brief A service for backtracking reconstruction information to its truth information
//
// \author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
//
////////////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "larsim/MCCheater/BackTrackerService.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace cheat{

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService( const fhicl::ParameterSet& pSet,  art::ActivityRegistry& reg)
    :BackTracker(
        pSet.get<fhicl::ParameterSet>("BackTracker"),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>()
        )
  {
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
  }

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg)
    :BackTracker(
        config.BackTrackerTable(),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>()
        )
  {
    //This line is only to be included until art/LArSoft discussion about lazy rebuild can be completed.
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
  }


  ////////////////////////////////////////////////
  //////////Event Rebuild Implimentation//////////
  /// This section contains the implimentation ///
  /// of all Prep templates from the service   ///
  /// provider.                                ///
  ////////////////////////////////////////////////

  //-----Temp rebuild function.
  void BackTrackerService::Rebuild( const art::Event& evt ){
    this->priv_PrepEvent(evt, art::ScheduleContext::invalid());
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepEvent( const art::Event& evt, art::ScheduleContext ){
    fEvt=&evt;
    BackTracker::ClearEvent();
    if( ! this->priv_CanRun(evt) ){ return; }
    this->priv_PrepSimChannels(evt);
    fEvt=nullptr; //don't save the pointer because it will be useless after this anyways. I want to make sure calls at the wrong time crash.
  }

  //---------------------------------------------------------------------
  bool BackTrackerService::priv_CanRun(const art::Event& evt){
    return BackTracker::CanRun(evt);
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepFailed(){
    throw cet::exception("BackTrackerService")
      <<"The called Prep Function failed. This is most likely due "
      <<"to a BackTracker function being used on a file that is real data.";
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepSimChannels(const art::Event& evt){
    if( !this->priv_CanRun(evt)) { this->priv_PrepFailed(); }
    if( this->priv_SimChannelsReady()){ return; }
    try{BackTracker::PrepSimChannels(evt);}
    catch(...){ mf::LogWarning("BackTrackerService") << "Rebuild failed to get the SimChannels. This is expected when running on a generation or simulation step.";}
  }

  //---------------------------------------------------------------------
  /*  void BackTrackerService::priv_PrepAllHitList(){
      if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
      if( this->priv_AllHitListReady() ){ return; }
      BackTracker::PrepAllHitList(*fEvt);
      }
      */
  /////////////////////////////////////////////
  // End of the Event Rebuild Implimentation //
  /////////////////////////////////////////////

  ///////////////////////////////////////////////
  /////BackTracking Functions Implimentation/////
  // This section contains the implimentation  //
  // of all BackTrackerService end user        //
  // functionality                             //
  ///////////////////////////////////////////////

  //---------------------------------------------------------------------
  const std::vector < art::Ptr < sim::SimChannel > >& BackTrackerService::SimChannels() const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::SimChannels();
  }

  //---------------------------------------------------------------------
  std::vector < const sim::IDE* > BackTrackerService::TrackIdToSimIDEs_Ps(int const& id) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::TrackIdToSimIDEs_Ps(id);
  }

  //---------------------------------------------------------------------
  std::vector < const sim::IDE* > BackTrackerService::TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::TrackIdToSimIDEs_Ps(id, view);
  }

  //---------------------------------------------------------------------
  art::Ptr < sim::SimChannel > BackTrackerService::FindSimChannel( raw::ChannelID_t channel ) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::FindSimChannel(channel);
  }

  //---------------------------------------------------------------------
  std::vector < sim::TrackIDE > BackTrackerService::ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time,    const double hit_end_time) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::ChannelToTrackIDEs(channel, hit_start_time, hit_end_time);
  }

  //---------------------------------------------------------------------
  std::vector < sim::TrackIDE >  BackTrackerService::HitToTrackIDEs(recob::Hit const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  std::vector < sim::TrackIDE >  BackTrackerService::HitToTrackIDEs(art::Ptr < recob::Hit > const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  std::vector < int >  BackTrackerService::HitToTrackIds(recob::Hit const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToTrackIds(hit);
  }

  //---------------------------------------------------------------------
  std::vector < sim::TrackIDE >  BackTrackerService::HitToEveTrackIDEs(recob::Hit const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToEveTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  std::vector < sim::TrackIDE >  BackTrackerService::HitToEveTrackIDEs(art::Ptr < recob::Hit > const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToEveTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  std::vector < art::Ptr < recob::Hit > > BackTrackerService::TrackIdToHits_Ps( const int& tkId, std::vector < art::Ptr < recob::Hit > > const& hitsIn ) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::TrackIdToHits_Ps( tkId, hitsIn);
  }

  //---------------------------------------------------------------------
  /*  const std::vector < art::Ptr < recob::Hit > > BackTrackerService::TrackIdToHits_Ps( const int& tkId ) {
  //Removed until Lazy Rebuild works
  //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
  //    if(!this->priv_AllHitListReady()){this->priv_PrepAllHitList();}
  return BackTracker::TrackIdToHits_Ps( tkId);
  }*/ //I can't support this functino and caching all the hits without lazy implimentation

  //---------------------------------------------------------------------
  std::vector < std::vector < art::Ptr < recob::Hit > > > BackTrackerService::TrackIdsToHits_Ps( std::vector < int > const&      tkIds, std:: vector < art::Ptr < recob::Hit > > const& hitsIn ) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::TrackIdsToHits_Ps( tkIds, hitsIn);
  }

  //---------------------------------------------------------------------
  /*  const std::vector < std::vector < art::Ptr < recob::Hit > > > BackTrackerService::TrackIdsToHits_Ps( std::vector < int > const&      tkIds ) {
  //Removed until Lazy Rebuild works
  //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
  //    if(!this->priv_AllHitListReady()){this->priv_PrepAllHitList();}
  return BackTracker::TrackIdsToHits_Ps( tkIds);
  }*/

  //---------------------------------------------------------------------
  std::vector< sim::IDE > BackTrackerService::HitToAvgSimIDEs ( recob::Hit const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToAvgSimIDEs(hit);
  }

  //---------------------------------------------------------------------
  std::vector< sim::IDE > BackTrackerService::HitToAvgSimIDEs ( art::Ptr<recob::Hit> hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToAvgSimIDEs(hit);
  }

  //---------------------------------------------------------------------
  std::vector< const sim::IDE* > BackTrackerService::HitToSimIDEs_Ps (recob::Hit const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToSimIDEs_Ps(hit);
  }

  //---------------------------------------------------------------------
  std::vector< const sim::IDE* > BackTrackerService::HitToSimIDEs_Ps (art::Ptr< recob::Hit > const& hit) const {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToSimIDEs_Ps(hit);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::SimIDEsToXYZ( std::vector< sim::IDE > const& ides) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::SimIDEsToXYZ(ides);
  }


  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::SimIDEsToXYZ( std::vector< const sim::IDE* > const& ide_Ps) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::SimIDEsToXYZ(ide_Ps);
  }



  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::HitToXYZ(const recob::Hit& hit) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitToXYZ(hit);
  }


  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::HitToXYZ(art::Ptr<recob::Hit> const& hit) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return HitToXYZ(*hit);
  }



  //---------------------------------------------------------------------
  double BackTrackerService::HitCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const&      hits) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitCollectionPurity(trackIds, hits);
  }


  //---------------------------------------------------------------------
  double BackTrackerService::HitChargeCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> >       const&     hits) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitChargeCollectionPurity(trackIds, hits);
  }



  //---------------------------------------------------------------------
  double BackTrackerService::HitCollectionEfficiency( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const&  hits, std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitCollectionEfficiency(trackIds, hits, allhits, view);
  }



  //---------------------------------------------------------------------
  double BackTrackerService::HitChargeCollectionEfficiency( std::set<int> trackIds, std::vector< art::Ptr<recob::Hit> >          const&        hits,        std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::HitChargeCollectionEfficiency(trackIds, hits, allhits, view);
  }



  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfTrackIds() const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::GetSetOfTrackIds();
  }


  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfEveIds() const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::GetSetOfEveIds();
  }



  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::GetSetOfTrackIds(hits);
  }


  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const
  {
    //Removed until Lazy Rebuild works
    //if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return BackTracker::GetSetOfEveIds(hits);
  }

  //---------------------------------------------------------------------
  std::vector< double> BackTrackerService::SpacePointHitsToWeightedXYZ(std::vector<art::Ptr<recob::Hit>> const& hits) const {
    //if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    return BackTracker::SpacePointHitsToWeightedXYZ(hits);
  }

  //---------------------------------------------------------------------
  std::vector< art::Ptr< recob::Hit > > BackTrackerService::SpacePointToHits_Ps(art::Ptr<recob::SpacePoint> const&   spt) const {
    //Note, this function is goofy in that it doesn't use the hits prepared by the module. That should likely be changed to make things more uniform.
    //if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    if( ! fEvt ){
      throw cet::exception("BackTrackerService")
        <<"This function is not yet implimented pending the implimentation of backtracker lazy loading.";
    }
    return BackTracker::SpacePointToHits_Ps(spt, *fEvt);
  }

  //---------------------------------------------------------------------
  std::vector< double > BackTrackerService::SpacePointToXYZ( art::Ptr< recob::SpacePoint > const& spt) const {
    //if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    if( ! fEvt ){
      throw cet::exception("BackTrackerService")
        <<"This function is not yet implimented pending the implimentation of backtracker lazy loading.";
    }
    return BackTracker::SpacePointToXYZ(spt, *fEvt);
  }

  DEFINE_ART_SERVICE(BackTrackerService)


}//end namespace cheat
