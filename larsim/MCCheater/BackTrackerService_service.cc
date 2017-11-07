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

#include "BackTrackerService.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace cheat{

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService( const fhicl::ParameterSet& pSet,  art::ActivityRegistry& reg)
    :fBackTracker(
        pSet.get<fhicl::ParameterSet>("BackTracker"), 
        lar::providerFrom<cheat::ParticleInventoryService>(), 
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>())
  {
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
  }

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg)
    :fBackTracker(
        config.BackTrackerTable(),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>())
  {
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
  }


  ////////////////////////////////////////////////
  //////////Event Rebuild Implimentation//////////
  /// This section contains the implimentation ///
  /// of all Prep templates from the service   ///
  /// provider.                                ///
  ////////////////////////////////////////////////

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepEvent( const art::Event& evt ){
    fEvt=&evt;
    fBackTracker.ClearEvent();
  }

  //---------------------------------------------------------------------
  bool BackTrackerService::priv_CanRun(const art::Event& evt){
    return fBackTracker.CanRun(evt);
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepFailed(){
    throw cet::exception("BackTrackerService")
      <<"The called Prep Function failed. This is most likely due "
      <<"to a BackTracker function being used on a file that is real data.";
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepSimChannels(){
    if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    if( this->priv_SimChannelsReady()){ return; }
    fBackTracker.PrepSimChannels(*fEvt);
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepAllHitList(){
    if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    if( this->priv_AllHitListReady() ){ return; }
    fBackTracker.PrepAllHitList(*fEvt);
  }

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
  const std::vector < art::Ptr < sim::SimChannel > >& BackTrackerService::SimChannels() {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.SimChannels(); 
  }
      
  //---------------------------------------------------------------------
  const std::vector < const sim::IDE* > BackTrackerService::TrackIdToSimIDEs_Ps(int const& id) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.TrackIdToSimIDEs_Ps(id);
  }

  //---------------------------------------------------------------------
  const std::vector < const sim::IDE* > BackTrackerService::TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.TrackIdToSimIDEs_Ps(id, view);
  }

  //---------------------------------------------------------------------
  const art::Ptr < sim::SimChannel > BackTrackerService::FindSimChannel( raw::ChannelID_t channel ){
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.FindSimChannel(channel);
  }

  //---------------------------------------------------------------------
  const std::vector < sim::TrackIDE > BackTrackerService::ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time,    const double hit_end_time){
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.ChannelToTrackIDEs(channel, hit_start_time, hit_end_time);
  }

  //---------------------------------------------------------------------
  const std::vector < sim::TrackIDE >  BackTrackerService::HitToTrackIDEs(recob::Hit const& hit){
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  const std::vector < sim::TrackIDE >  BackTrackerService::HitToTrackIDEs(art::Ptr < recob::Hit > const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  const std::vector < int >  BackTrackerService::HitToTrackIds(recob::Hit const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToTrackIds(hit);
  }

  //---------------------------------------------------------------------
  const std::vector < sim::TrackIDE >  BackTrackerService::HitToEveTrackIDEs(recob::Hit const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToEveTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  const std::vector < sim::TrackIDE >  BackTrackerService::HitToEveTrackIDEs(art::Ptr < recob::Hit > const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToEveTrackIDEs(hit);
  }

  //---------------------------------------------------------------------
  const std::vector < art::Ptr < recob::Hit > > BackTrackerService::TrackIdToHits_Ps( const int& tkId, std::vector < art::Ptr < recob::Hit > > const& hitsIn ) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.TrackIdToHits_Ps( tkId, hitsIn);
  }

  //---------------------------------------------------------------------
  const std::vector < art::Ptr < recob::Hit > > BackTrackerService::TrackIdToHits_Ps( const int& tkId ) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    if(!this->priv_AllHitListReady()){this->priv_PrepAllHitList();}
    return fBackTracker.TrackIdToHits_Ps( tkId);
  }

  //---------------------------------------------------------------------
  const std::vector < std::vector < art::Ptr < recob::Hit > > > BackTrackerService::TrackIdsToHits_Ps( std::vector < int > const&      tkIds, std:: vector < art::Ptr < recob::Hit > > const& hitsIn ) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.TrackIdsToHits_Ps( tkIds, hitsIn);
  }

  //---------------------------------------------------------------------
  const std::vector < std::vector < art::Ptr < recob::Hit > > > BackTrackerService::TrackIdsToHits_Ps( std::vector < int > const&      tkIds ) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    if(!this->priv_AllHitListReady()){this->priv_PrepAllHitList();}
    return fBackTracker.TrackIdsToHits_Ps( tkIds);
  }

  //---------------------------------------------------------------------
  const std::vector< sim::IDE > BackTrackerService::HitToAvgSimIDEs ( recob::Hit const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToAvgSimIDEs(hit);
  }

  //---------------------------------------------------------------------
  const std::vector< sim::IDE > BackTrackerService::HitToAvgSimIDEs ( art::Ptr<recob::Hit> hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToAvgSimIDEs(hit);
  }

  //---------------------------------------------------------------------
  const std::vector< const sim::IDE* > BackTrackerService::HitToSimIDEs_Ps (recob::Hit const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToSimIDEs_Ps(hit);
  }

  //---------------------------------------------------------------------
  const std::vector< const sim::IDE* > BackTrackerService::HitToSimIDEs_Ps (art::Ptr< recob::Hit > const& hit) {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToSimIDEs_Ps(hit);
  }

  //---------------------------------------------------------------------
  const std::vector<double> BackTrackerService::SimIDEsToXYZ( std::vector< sim::IDE > const& ides) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.SimIDEsToXYZ(ides);
  }


  //---------------------------------------------------------------------
  const std::vector<double> BackTrackerService::SimIDEsToXYZ( std::vector< const sim::IDE* > const& ide_Ps) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.SimIDEsToXYZ(ide_Ps);
  }



  //---------------------------------------------------------------------
  const std::vector<double> BackTrackerService::HitToXYZ(const recob::Hit& hit) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitToXYZ(hit);
  }


  //---------------------------------------------------------------------
  const std::vector<double> BackTrackerService::HitToXYZ(art::Ptr<recob::Hit> const& hit) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return HitToXYZ(hit);
  }



  //---------------------------------------------------------------------
  const double BackTrackerService::HitCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const&      hits) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitCollectionPurity(trackIds, hits);
  }


  //---------------------------------------------------------------------
  const double BackTrackerService::HitChargeCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> >       const&     hits) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitChargeCollectionPurity(trackIds, hits);
  }



  //---------------------------------------------------------------------
  const double BackTrackerService::HitCollectionEfficiency( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const&  hits, std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view)
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitCollectionEfficiency(trackIds, hits, allhits, view);
  }



  //---------------------------------------------------------------------
  const double BackTrackerService::HitChargeCollectionEfficiency( std::set<int> trackIds, std::vector< art::Ptr<recob::Hit> >          const&        hits,        std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.HitChargeCollectionEfficiency(trackIds, hits, allhits, view);
  }



  //---------------------------------------------------------------------
  const std::set<int> BackTrackerService::GetSetOfTrackIds() 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.GetSetOfTrackIds();
  }


  //---------------------------------------------------------------------
  const std::set<int> BackTrackerService::GetSetOfEveIds() 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.GetSetOfEveIds();
  }



  //---------------------------------------------------------------------
  const std::set<int> BackTrackerService::GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits ) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.GetSetOfTrackIds(hits);
  }


  //---------------------------------------------------------------------
  const std::set<int> BackTrackerService::GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits ) 
  {
    if(!this->priv_SimChannelsReady()){this->priv_PrepSimChannels();}
    return fBackTracker.GetSetOfEveIds(hits);
  }

  //---------------------------------------------------------------------
  const std::vector< double> BackTrackerService::SpacePointHitsToWeightedXYZ(std::vector<art::Ptr<recob::Hit>> const& hits) {
    if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    return fBackTracker.SpacePointHitsToWeightedXYZ(hits);
  }

  //---------------------------------------------------------------------
  const std::vector< art::Ptr< recob::Hit > > BackTrackerService::SpacePointToHits_Ps(art::Ptr<recob::SpacePoint> const&   spt) {
    //Note, this function is goofy in that it doesn't use the hits prepared by the module. That should likely be changed to make things more uniform.
    if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    return fBackTracker.SpacePointToHits_Ps(spt, *fEvt);
  }

  //---------------------------------------------------------------------
  const std::vector< double > BackTrackerService::SpacePointToXYZ( art::Ptr< recob::SpacePoint > const& spt) {
    if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    return fBackTracker.SpacePointToXYZ(spt, *fEvt);
  }

    DEFINE_ART_SERVICE(BackTrackerService)


}//end namespace cheat
