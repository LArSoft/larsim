////////////////////////////////////////////////////////////////////////
//
// \file: PhotonBackTrackerServiceService_service.cc
//
//jason.stock@mines.sdsmt.edu
//Based on the BackTracker_service by Brian Rebel
//
////////////////////////////////////////////////////////////////////////

#include "larsim/MCCheater/PhotonBackTrackerService.h"

//CPP includes
#include <map>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/Simulation/SimListUtils.h"

//#include "lardata/Utilities/AssociationUtil.h"
//#include "lardataobj/Simulation/sim.h"

namespace cheat{

  //----------------------------------------------------------------------
  PhotonBackTrackerService::PhotonBackTrackerService(const fhicl::ParameterSet& pSet, art::ActivityRegistry& reg)
    :PhotonBackTracker(
        pSet.get<fhicl::ParameterSet>("PhotonBackTracker"),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>()
        )
  {
    reg.sPreProcessEvent.watch(this, &PhotonBackTrackerService::priv_PrepEvent);
  }

  //----------------------------------------------------------------------
  PhotonBackTrackerService::PhotonBackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg)
    :PhotonBackTracker(
        config.PhotonBackTrackerTable(),
        lar::providerFrom<cheat::ParticleInventoryService>(),
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>()
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
  void PhotonBackTrackerService::Rebuild(const art::Event& evt)
  {
    this->priv_PrepEvent(evt);
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepEvent( const art::Event& evt)
  {
    PhotonBackTracker::ClearEvent();
    if( ! this->priv_CanRun(evt) ){ return; }
    this->priv_PrepOpDetBTRs(evt);
  }

  //----------------------------------------------------------------------
  bool PhotonBackTrackerService::priv_CanRun(const art::Event& evt){
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
  void PhotonBackTrackerService::priv_PrepOpDetBTRs(const art::Event& evt){
    if( !this->priv_CanRun(evt) ) {this->priv_PrepFailed(); }
    if( this->priv_OpDetBTRsReady()){ return; }
    try{PhotonBackTracker::PrepOpDetBTRs(evt);}
    catch(...){ mf::LogWarning("PhotonBackTrackerService")
      <<"Rebuild failed to get the OpDetBTRs. This is expected when "
        <<"running on a generation or simulation step.";}
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
  const art::Ptr< sim::OpDetBacktrackerRecord > PhotonBackTrackerService::FindOpDetBTR(int const& opDetNum) 
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
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& opHit)
  {
    return PhotonBackTracker::OpHitToTrackSDPs(opHit);
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
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToEveTrackSDPs(art::Ptr<recob::OpHit> const& opHit)
  {
    return PhotonBackTracker::OpHitToEveTrackSDPs(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToEveTrackSDPs(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToEveTrackSDPs(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector<art::Ptr<recob::OpHit>> PhotonBackTrackerService::TrackIdToOpHits_Ps(const int& tkId, std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
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
  const std::vector< const sim::SDP* > PhotonBackTrackerService::OpHitToSimSDPs_Ps(art::Ptr<recob::OpHit> const& opHit) 
  {
    return PhotonBackTracker::OpHitToSimSDPs_Ps(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&
  {
    return PhotonBackTracker::SimSDPsToXYZ(sdps);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::SimSDPsToXYZ(std::vector<const sim::SDP*> const& sdps)
  {
    return PhotonBackTracker::SimSDPsToXYZ(sdps);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::OpHitToXYZ(const recob::OpHit& opHit)
  {
    return PhotonBackTracker::OpHitToXYZ(opHit);
  }

  //----------------------------------------------------------------------
  const std::vector< double> PhotonBackTrackerService::OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit)
  {
    return PhotonBackTracker::OpHitToXYZ(opHit);
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
  const std::set< int> PhotonBackTrackerService::GetSetOfEveIds(std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    return PhotonBackTracker::GetSetOfEveIds(opHits);
  }

  //----------------------------------------------------------------------
  const std::set< int> PhotonBackTrackerService::GetSetOfTrackIds(std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    return PhotonBackTracker::GetSetOfTrackIds(opHits);
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitCollectionPurity(std::set<int> const& tkIds, 
      std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    return PhotonBackTracker::OpHitCollectionPurity(tkIds, opHits);
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitLightCollectionPurity(std::set<int> const& tkIds, 
      std::vector< art::Ptr<recob::OpHit> > const& opHits)
  {
    return PhotonBackTracker::OpHitLightCollectionPurity(tkIds, opHits);
  }


  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitCollectionEfficiency(std::set<int> const& tkIds, 
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn,
      geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTrackerService")<<"This function is not supported. OpHits do not have type View.\n";
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitCollectionEfficiency(std::set<int> const& tkIds, 
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn)
  {
    return PhotonBackTracker::OpHitCollectionEfficiency(tkIds, opHits, opHitsIn);
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitLightCollectionEfficiency(std::set<int> const& tkIds, 
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn,
      geo::View_t const& view)
  {
    throw cet::exception("PhotonBackTrackerService")<<"This function is not supported. OpHits do not have type View.\n";
  }

  //----------------------------------------------------------------------
  const double PhotonBackTrackerService::OpHitChargeCollectionEfficiency(std::set<int> const& tkIds, 
      std::vector< art::Ptr<recob::OpHit> > const& opHits,
      std::vector< art::Ptr<recob::OpHit> > const& opHitsIn)
  {
    return PhotonBackTracker::OpHitLightCollectionEfficiency(tkIds, opHits, opHitsIn);
  }




  DEFINE_ART_SERVICE(PhotonBackTrackerService)
    //-------------------------------------------------------------------
} // namespace

