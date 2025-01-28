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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"             // for provid...
#include "larcore/Geometry/WireReadout.h"              // for Geometry
#include "larsim/MCCheater/ParticleInventoryService.h" // for Partic...

namespace cheat {

  //----------------------------------------------------------------------
  PhotonBackTrackerService::PhotonBackTrackerService(fhicl::ParameterSet const& pSet,
                                                     art::ActivityRegistry& reg)
    : PhotonBackTracker(pSet.get<fhicl::ParameterSet>("PhotonBackTracker"),
                        lar::providerFrom<cheat::ParticleInventoryService>(),
                        &art::ServiceHandle<geo::WireReadout const>()->Get())
  {
    reg.sPreProcessEvent.watch(this, &PhotonBackTrackerService::priv_PrepEvent);
  }

  //----------------------------------------------------------------------
  PhotonBackTrackerService::PhotonBackTrackerService(fhiclConfig const& config,
                                                     art::ActivityRegistry& reg)
    : PhotonBackTracker(config.PhotonBackTrackerTable(),
                        lar::providerFrom<cheat::ParticleInventoryService>(),
                        &art::ServiceHandle<geo::WireReadout const>()->Get())
  {
    reg.sPreProcessEvent.watch(this, &PhotonBackTrackerService::priv_PrepEvent);
  }

  //----------------------------------------------------------------------
  // Event Rebuild Implementation: This section contains the implementation
  // of all Prep templates from the service  provider.

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::Rebuild(art::Event const& evt)
  {
    priv_PrepEvent(evt, art::ScheduleContext::invalid());
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepEvent(art::Event const& evt, art::ScheduleContext)
  {
    PhotonBackTracker::ClearEvent();
    if (!priv_CanRun(evt)) { return; }
    priv_PrepOpDetBTRs(evt);
    priv_PrepOpFlashToOpHits(evt);
  }

  //----------------------------------------------------------------------
  bool PhotonBackTrackerService::priv_CanRun(art::Event const& evt)
  {
    return PhotonBackTracker::CanRun(evt);
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepFailed()
  {
    throw cet::exception("PhotonBackTrackerService")
      << "The called Prep Fucntion failed. This is most likely due "
      << "to a PhotonBackTracker fucntion being used on a file that "
      << "is real data.";
  }

  //----------------------------------------------------------------------
  void PhotonBackTrackerService::priv_PrepOpDetBTRs(art::Event const& evt)
  {
    if (!priv_CanRun(evt)) { priv_PrepFailed(); }
    if (priv_OpDetBTRsReady()) { return; }
    try {
      PhotonBackTracker::PrepOpDetBTRs(evt);
    }
    catch (cet::exception const&) {
      //This needs to go. Make it specific if there is a really an exception we would like to catch.
      mf::LogWarning("PhotonBackTrackerService")
        << "Rebuild failed to get the OpDetBTRs. This is expected when "
        << "running on a generation or simulation step.";
    }
  }

  void PhotonBackTrackerService::priv_PrepOpFlashToOpHits(art::Event const& evt)
  {
    if (!priv_CanRun(evt)) { priv_PrepFailed(); }
    if (priv_OpFlashToOpHitsReady()) { return; }
    try {
      PhotonBackTracker::PrepOpFlashToOpHits(evt);
    }
    catch (cet::exception const&) {
      // This needs to go. Make it specific if there is a really an exception we would like to catch.
      mf::LogWarning("PhotonBackTrackerService")
        << "Rebuild failed to get the OpFlashToOpHits. This is expected when "
        << "running on a generation or simulation stage.";
    }
  }

  //----------------------------------------------------------------------
  // BackTracking Function Implementation: This section contains the implimentation of all
  // PhotonBackTrackerService end user functionality

  //----------------------------------------------------------------------
  const std::vector<art::Ptr<sim::OpDetBacktrackerRecord>>& PhotonBackTrackerService::OpDetBTRs()
  {
    return PhotonBackTracker::OpDetBTRs();
  }

  //----------------------------------------------------------------------
  double PhotonBackTrackerService::GetDelay()
  {
    return PhotonBackTracker::GetDelay();
  }

  //----------------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTrackerService::TrackIdToSimSDPs_Ps(int const id)
  {
    return PhotonBackTracker::TrackIdToSimSDPs_Ps(id);
  }

  //----------------------------------------------------------------------
  art::Ptr<sim::OpDetBacktrackerRecord> PhotonBackTrackerService::FindOpDetBTR(int const opDetNum)
  {
    return PhotonBackTracker::FindOpDetBTR(opDetNum);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpDetToTrackSDPs(
    int const OpDetNum,
    double const opHit_start_time,
    double const opHit_end_time)
  {
    return PhotonBackTracker::OpDetToTrackSDPs(OpDetNum, opHit_start_time, opHit_end_time);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToTrackSDPs(
    art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToTrackSDPs(opHit_P);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToTrackSDPs(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToTrackSDPs(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<int> PhotonBackTrackerService::OpHitToTrackIds(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToTrackIds(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<int> PhotonBackTrackerService::OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToTrackIds(opHit_P);
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  std::vector<int> PhotonBackTrackerService::OpHitToEveTrackIds(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToEveTrackIds(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<int> PhotonBackTrackerService::OpHitToEveTrackIds(
    art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToEveTrackIds(opHit_P);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToEveTrackSDPs(
    art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToEveTrackSDPs(opHit_P);
  }

  //----------------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTrackerService::OpHitToEveTrackSDPs(
    recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToEveTrackSDPs(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<art::Ptr<recob::OpHit>> PhotonBackTrackerService::TrackIdToOpHits_Ps(
    int const tkId,
    std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
  {
    return PhotonBackTracker::TrackIdToOpHits_Ps(tkId, hitsIn);
  }

  //----------------------------------------------------------------------
  std::vector<std::vector<art::Ptr<recob::OpHit>>> PhotonBackTrackerService::TrackIdsToOpHits_Ps(
    std::vector<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& hitsIn)
  {
    return PhotonBackTracker::TrackIdsToOpHits_Ps(tkIds, hitsIn);
  }

  //----------------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTrackerService::OpHitToSimSDPs_Ps(
    recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToSimSDPs_Ps(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTrackerService::OpHitToSimSDPs_Ps(
    art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToSimSDPs_Ps(opHit_P);
  }

  //----------------------------------------------------------------------
  std::unordered_set<const sim::SDP*> PhotonBackTrackerService::OpHitToEveSimSDPs_Ps(
    recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToEveSimSDPs_Ps(opHit);
  }

  //----------------------------------------------------------------------
  std::unordered_set<const sim::SDP*> PhotonBackTrackerService::OpHitToEveSimSDPs_Ps(
    art::Ptr<recob::OpHit>& opHit_P)
  {
    return PhotonBackTracker::OpHitToEveSimSDPs_Ps(opHit_P);
  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTrackerService::SimSDPsToXYZ(
    std::vector<sim::SDP> const& sdps) const&
  {
    return PhotonBackTracker::SimSDPsToXYZ(sdps);
  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTrackerService::SimSDPsToXYZ(
    std::vector<const sim::SDP*> const& sdps_Ps)
  {
    return PhotonBackTracker::SimSDPsToXYZ(sdps_Ps);
  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTrackerService::OpHitToXYZ(recob::OpHit const& opHit)
  {
    return PhotonBackTracker::OpHitToXYZ(opHit);
  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTrackerService::OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit_P)
  {
    return PhotonBackTracker::OpHitToXYZ(opHit_P);
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::GetSetOfEveIds()
  {
    return PhotonBackTracker::GetSetOfEveIds();
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::GetSetOfTrackIds()
  {
    return PhotonBackTracker::GetSetOfTrackIds();
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::GetSetOfEveIds(
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps)
  {
    return PhotonBackTracker::GetSetOfEveIds(opHits_Ps);
  }
  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::GetSetOfEveIds(std::vector<recob::OpHit> const& opHits)
  {
    return PhotonBackTracker::GetSetOfEveIds(opHits);
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::GetSetOfTrackIds(
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps)
  {
    return PhotonBackTracker::GetSetOfTrackIds(opHits_Ps);
  }
  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::GetSetOfTrackIds(std::vector<recob::OpHit> const& opHits)
  {
    return PhotonBackTracker::GetSetOfTrackIds(opHits);
  }

  //----------------------------------------------------------------------
  double PhotonBackTrackerService::OpHitCollectionPurity(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps)
  {
    return PhotonBackTracker::OpHitCollectionPurity(tkIds, opHits_Ps);
  }

  //----------------------------------------------------------------------
  double PhotonBackTrackerService::OpHitLightCollectionPurity(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps)
  {
    return PhotonBackTracker::OpHitLightCollectionPurity(tkIds, opHits_Ps);
  }

  //----------------------------------------------------------------------
  double PhotonBackTrackerService::OpHitCollectionEfficiency(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps,
    std::vector<art::Ptr<recob::OpHit>> const& opHitsIn_Ps)
  {
    return PhotonBackTracker::OpHitCollectionEfficiency(tkIds, opHits_Ps, opHitsIn_Ps);
  }

  //----------------------------------------------------------------------
  double PhotonBackTrackerService::OpHitLightCollectionEfficiency(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps,
    std::vector<art::Ptr<recob::OpHit>> const& opHitsIn_Ps)
  {
    return PhotonBackTracker::OpHitLightCollectionEfficiency(tkIds, opHits_Ps, opHitsIn_Ps);
  }

  //----------------------------------------------------------------------
  std::set<int> PhotonBackTrackerService::OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P) const
  {
    return PhotonBackTracker::OpFlashToTrackIds(flash_P);
  }

  //----------------------------------------------------------------------
  std::vector<art::Ptr<recob::OpHit>> PhotonBackTrackerService::OpFlashToOpHits_Ps(
    art::Ptr<recob::OpFlash>& flash_P)
  {
    return PhotonBackTracker::OpFlashToOpHits_Ps(flash_P);
  }

  //----------------------------------------------------------------------
  std::vector<double> PhotonBackTrackerService::OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P)
  {
    return PhotonBackTracker::OpFlashToXYZ(flash_P);
  }

} // namespace
