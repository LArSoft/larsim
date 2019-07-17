////////////////////////////////////////////////////////////////////////
//
// \file PhotonBackTrackerServiceService_service.cc
// \brief This is the LArSoft facing implimentation of the PhotonBackTracker.
// \author jason.stock@mines.sdsmt.edu
//
//
////////////////////////////////////////////////////////////////////////

#include "larsim/MCCheater/PhotonBackTrackerService.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"                     // for provid...
#include "larcore/Geometry/Geometry.h"                         // for Geometry
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"      // for View_t
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"      // for SDP (p...
#include "larsim/MCCheater/ParticleInventoryService.h"         // for Partic...
#include "larsim/MCCheater/PhotonBackTracker.h"                // for Photon...

namespace cheat{

  //----------------------------------------------------------------------
  PhotonBackTrackerService::PhotonBackTrackerService(fhicl::ParameterSet const& pSet, art::ActivityRegistry& reg)
    :PhotonBackTracker(
        pSet.get<fhicl::ParameterSet>("PhotonBackTracker"),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>()//,
//        lar::providerFrom<detinfo::DetectorClocksService>()
        )
  {
    reg.sPreProcessEvent.watch(this, &PhotonBackTrackerService::priv_PrepEvent);
  }

  //----------------------------------------------------------------------
  PhotonBackTrackerService::PhotonBackTrackerService(fhiclConfig const& config, art::ActivityRegistry& reg)
    :PhotonBackTracker(
        config.PhotonBackTrackerTable(),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>()//,
//        lar::providerFrom<detinfo::DetectorClocksService>()
        )
  {
    reg.sPreProcessEvent.watch(this, &PhotonBackTrackerService::priv_PrepEvent);
  }

  ////////////////////////////////////////////////
  //////////Event Rebuild Implimentation//////////
  /// This section contains the implimentation ///
  /// of all Prep templates from the service   ///
  /// provider.                                ///
  ////////////////////////////////////////////////


  //----------------------------------------------------------------------
  void PhotonBackTrackerService::Rebuild(art::Event const& evt)
  {
    this->priv_PrepEvent(evt, art::ScheduleContext::invalid());
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepEvent( art::Event const& evt, art::ScheduleContext)
  {
    PhotonBackTracker::ClearEvent();
    if( ! this->priv_CanRun(evt) ){ return; }
    this->priv_PrepOpDetBTRs(evt);
    this->priv_PrepOpFlashToOpHits(evt);
  }

  //----------------------------------------------------------------------
  bool PhotonBackTrackerService::priv_CanRun(art::Event const& evt){
    return PhotonBackTracker::CanRun(evt);
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepFailed(){
    throw cet::exception("PhotonBackTrackerService")
      <<"The called Prep Fucntion failed. This is most likely due "
      <<"to a PhotonBackTracker fucntion being used on a file that "
      <<"is real data.";
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepOpDetBTRs(art::Event const& evt){
    if( !this->priv_CanRun(evt) ) {this->priv_PrepFailed(); }
    if( this->priv_OpDetBTRsReady()){ return; }
    try{PhotonBackTracker::PrepOpDetBTRs(evt);}
    //catch(...){ mf::LogWarning("PhotonBackTrackerService")//This needs to go. Catch all should not be used.
    catch(cet::exception e){//This needs to go. Make it specific if there is a really an exception we would like to catch.
      mf::LogWarning("PhotonBackTrackerService")
        <<"Rebuild failed to get the OpDetBTRs. This is expected when "
        <<"running on a generation or simulation step.";}
  }

  void PhotonBackTrackerService::priv_PrepOpFlashToOpHits(art::Event const& evt){
    if( !this->priv_CanRun(evt) ) {this->priv_PrepFailed();}
    if( this->priv_OpFlashToOpHitsReady()){ return; }
    try{PhotonBackTracker::PrepOpFlashToOpHits(evt);}
    //catch(...){ //This needs to go. Catch all should not be used.
    catch(cet::exception e){//This needs to go. Make it specific if there is a really an exception we would like to catch.
      mf::LogWarning("PhotonBackTrackerService")
        <<"Rebuild failed to get the OpFlashToOpHits. This is expected when "
        <<"running on a generation or simulation stage.";
    }
  }

  /////////////////////////////////////////////
  // End of the Event Rebuild Implimentation //
  /////////////////////////////////////////////

  ///////////////////////////////////////////////
  /////BackTracking Functions Implimentation/////
  // This section contains the implimentation  //
  // of all PhotonBackTrackerService end user  //
  // functionality                             //
  ///////////////////////////////////////////////


  //----------------------------------------------------------------------
  const std::vector< art::Ptr< sim::OpDetBacktrackerRecord >>& PhotonBackTrackerService::OpDetBTRs()
  {
    return PhotonBackTracker::OpDetBTRs();
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::GetDelay(){ return PhotonBackTracker::GetDelay();}

  //----------------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTrackerService::TrackIdToSimSDPs_Ps(int const& id)
  {
    return PhotonBackTracker::TrackIdToSimSDPs_Ps(id);
  }

  //----------------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTrackerService::TrackIdToSimSDPs_Ps(int const& id,  geo::View_t const& view )
  {
    return PhotonBackTracker::TrackIdToSimSDPs_Ps(id, view);
  }

  //----------------------------------------------------------------------
  art::Ptr< sim::OpDetBacktrackerRecord > PhotonBackTrackerService::FindOpDetBTR(int const& opDetNum)
  {
    return PhotonBackTracker::FindOpDetBTR(opDetNum);
  }

  //----------------------------------------------------------------------
  const std::vector < sim::TrackSDP > PhotonBackTrackerService::OpDetToTrackSDPs( int const& OpDetNum,
      double const& opHit_start_time,
      double const& opHit_end_time)
  {
    return PhotonBackTracker::OpDetToTrackSDPs(OpDetNum, opHit_start_time, opHit_end_time);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P )
  {
    return PhotonBackTracker::OpHitToTrackSDPs(opHit_P );
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToTrackSDPs(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToTrackSDPs(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector < int > PhotonBackTrackerService::OpHitToTrackIds(recob::OpHit const& opHit) {
    return PhotonBackTracker::OpHitToTrackIds(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector < int > PhotonBackTrackerService::OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit_P) {
    return PhotonBackTracker::OpHitToTrackIds(opHit_P);
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  const std::vector < int > PhotonBackTrackerService::OpHitToEveTrackIds(recob::OpHit const& opHit) {
    return PhotonBackTracker::OpHitToEveTrackIds(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector < int > PhotonBackTrackerService::OpHitToEveTrackIds(art::Ptr<recob::OpHit> const& opHit_P) {
    return PhotonBackTracker::OpHitToEveTrackIds(opHit_P);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToEveTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P )
  {
    return PhotonBackTracker::OpHitToEveTrackSDPs(opHit_P);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToEveTrackSDPs(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToEveTrackSDPs(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector<art::Ptr<recob::OpHit>> PhotonBackTrackerService::TrackIdToOpHits_Ps(int const& tkId, std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
  {
    return PhotonBackTracker::TrackIdToOpHits_Ps(tkId, hitsIn);
  }

  //----------------------------------------------------------------------
  const std::vector<std::vector<art::Ptr<recob::OpHit>>> PhotonBackTrackerService::TrackIdsToOpHits_Ps(std::vector<int> const& tkIds, std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
  {
    return PhotonBackTracker::TrackIdsToOpHits_Ps(tkIds, hitsIn);
  }

  //----------------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTrackerService::OpHitToSimSDPs_Ps(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToSimSDPs_Ps(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector< const sim::SDP* > PhotonBackTrackerService::OpHitToSimSDPs_Ps(art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToSimSDPs_Ps(opHit_P);
  }

  //----------------------------------------------------------------------
  /*
     const std::vector< sim::SDP > PhotonBackTrackerService::OpHitToChannelWeightedSimSDPs(art::Ptr<recob::OpHit> const& opHit_P)
     {
     return PhotonBackTracker::OpHitToChannelWeightedSimSDPs(opHit_P);
     }*/

  //----------------------------------------------------------------------
  const std::unordered_set< const sim::SDP* > PhotonBackTrackerService::OpHitToEveSimSDPs_Ps(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToEveSimSDPs_Ps(opHit);
  }

  //----------------------------------------------------------------------
  const std::unordered_set< const sim::SDP* > PhotonBackTrackerService::OpHitToEveSimSDPs_Ps(art::Ptr<recob::OpHit> & opHit_P)
  {
    return PhotonBackTracker::OpHitToEveSimSDPs_Ps(opHit_P);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&
  {
    return PhotonBackTracker::SimSDPsToXYZ(sdps);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::SimSDPsToXYZ(std::vector<const sim::SDP*> const& sdps_Ps )
  {
    return PhotonBackTracker::SimSDPsToXYZ(sdps_Ps);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::OpHitToXYZ(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToXYZ(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToXYZ(opHit_P);
  }

  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfEveIds()
  {
    return PhotonBackTracker::GetSetOfEveIds();
  }

  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfTrackIds()
  {
    return PhotonBackTracker::GetSetOfTrackIds();
  }

  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfEveIds(std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps)
  {
    return PhotonBackTracker::GetSetOfEveIds(opHits_Ps);
  }
  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfEveIds(std::vector< recob::OpHit > const& opHits)
  {
    return PhotonBackTracker::GetSetOfEveIds(opHits);
  }

  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfTrackIds(std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps)
  {
    return PhotonBackTracker::GetSetOfTrackIds(opHits_Ps);
  }
  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfTrackIds(std::vector< recob::OpHit > const& opHits)
  {
    return PhotonBackTracker::GetSetOfTrackIds(opHits);
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitCollectionPurity(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps)
  {
    return PhotonBackTracker::OpHitCollectionPurity(tkIds, opHits_Ps);
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitLightCollectionPurity(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps)
  {
    return PhotonBackTracker::OpHitLightCollectionPurity(tkIds, opHits_Ps);
  }


  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn_Ps,
      geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTrackerService")<<"This function is not supported. OpHits do not have type View.\n";
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn_Ps)
  {
    return PhotonBackTracker::OpHitCollectionEfficiency(tkIds, opHits_Ps, opHitsIn_Ps);
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn_Ps,
      geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTrackerService")<<"This function is not supported. OpHits do not have type View.\n";
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
      std::vector< art::Ptr<recob::OpHit> > const& opHits_Ps,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn_Ps)
  {
    return PhotonBackTracker::OpHitLightCollectionEfficiency(tkIds, opHits_Ps, opHitsIn_Ps);
  }

  //----------------------------------------------------------------------
  const std::set<int> PhotonBackTrackerService::OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P ) const{
    return PhotonBackTracker::OpFlashToTrackIds( flash_P);
  }

  //----------------------------------------------------------------------
  const std::vector < art::Ptr< recob::OpHit > > PhotonBackTrackerService::OpFlashToOpHits_Ps ( art::Ptr < recob::OpFlash > & flash_P ){
    return PhotonBackTracker::OpFlashToOpHits_Ps( flash_P);
  }

  //----------------------------------------------------------------------
  const std::vector < double > PhotonBackTrackerService::OpFlashToXYZ ( art::Ptr < recob::OpFlash > & flash_P ){
    return PhotonBackTracker::OpFlashToXYZ( flash_P );
  }


  DEFINE_ART_SERVICE(PhotonBackTrackerService)
    //-------------------------------------------------------------------
  } // namespace
