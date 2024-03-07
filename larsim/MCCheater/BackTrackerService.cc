////////////////////////////////////////////////////////////////////////////////////////
//
// \file BackTrackerService_service.cc
// \brief A service for backtracking reconstruction information to its truth
// information
//
// \author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
//
////////////////////////////////////////////////////////////////////////////////////////

#include "larsim/MCCheater/BackTrackerService.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cheat {

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService(const fhicl::ParameterSet& pSet,
                                         art::ActivityRegistry& reg)
    : BackTracker(pSet.get<fhicl::ParameterSet>("BackTracker"),
                  lar::providerFrom<cheat::ParticleInventoryService>(),
                  lar::providerFrom<geo::Geometry>(),
                  &art::ServiceHandle<geo::WireReadout const>()->Get())
  {
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
  }

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg)
    : BackTracker(config.BackTrackerTable(),
                  lar::providerFrom<cheat::ParticleInventoryService>(),
                  lar::providerFrom<geo::Geometry>(),
                  &art::ServiceHandle<geo::WireReadout const>()->Get())
  {
    // This line is only to be included until art/LArSoft discussion about lazy rebuild
    // can be completed.
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
  }

  //---------------------------------------------------------------------
  // Event Rebuild Implimentation: This section contains the implimentation of all Prep
  // templates from the service provider.

  //-----Temp rebuild function.
  void BackTrackerService::Rebuild(const art::Event& evt)
  {
    this->priv_PrepEvent(evt, art::ScheduleContext::invalid());
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepEvent(const art::Event& evt, art::ScheduleContext)
  {
    fEvt = &evt;
    BackTracker::ClearEvent();
    if (!this->priv_CanRun(evt)) { return; }
    this->priv_PrepSimChannels(evt);
    fEvt = nullptr; // don't save the pointer because it will be useless after this
                    // anyways. I want to make sure calls at the wrong time crash.
  }

  //---------------------------------------------------------------------
  bool BackTrackerService::priv_CanRun(const art::Event& evt)
  {
    return BackTracker::CanRun(evt);
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepFailed()
  {
    throw cet::exception("BackTrackerService")
      << "The called Prep Function failed. This is most likely due "
      << "to a BackTracker function being used on a file that is real data.";
  }

  //---------------------------------------------------------------------
  void BackTrackerService::priv_PrepSimChannels(const art::Event& evt)
  {
    if (!this->priv_CanRun(evt)) { this->priv_PrepFailed(); }
    if (this->priv_SimChannelsReady()) { return; }
    try {
      BackTracker::PrepSimChannels(evt);
    }
    catch (...) {
      mf::LogWarning("BackTrackerService")
        << "Rebuild failed to get the SimChannels. This is expected when running on a generation "
           "or simulation step.";
    }
  }

  //---------------------------------------------------------------------
  // BackTracking Functions Implementation: This section contains the implimentation of
  // all BackTrackerService end user functionality

  //---------------------------------------------------------------------
  const std::vector<art::Ptr<sim::SimChannel>>& BackTrackerService::SimChannels() const
  {
    return BackTracker::SimChannels();
  }

  //---------------------------------------------------------------------
  std::vector<const sim::IDE*> BackTrackerService::TrackIdToSimIDEs_Ps(int const& id) const
  {
    return BackTracker::TrackIdToSimIDEs_Ps(id);
  }

  //---------------------------------------------------------------------
  std::vector<const sim::IDE*> BackTrackerService::TrackIdToSimIDEs_Ps(int const& id,
                                                                       const geo::View_t view) const
  {
    return BackTracker::TrackIdToSimIDEs_Ps(id, view);
  }

  //---------------------------------------------------------------------
  art::Ptr<sim::SimChannel> BackTrackerService::FindSimChannel(raw::ChannelID_t channel) const
  {
    return BackTracker::FindSimChannel(channel);
  }

  //---------------------------------------------------------------------
  std::vector<sim::TrackIDE> BackTrackerService::ChannelToTrackIDEs(
    detinfo::DetectorClocksData const& clockData,
    raw::ChannelID_t channel,
    const double hit_start_time,
    const double hit_end_time) const
  {
    return BackTracker::ChannelToTrackIDEs(clockData, channel, hit_start_time, hit_end_time);
  }

  //---------------------------------------------------------------------
  std::vector<sim::TrackIDE> BackTrackerService::HitToTrackIDEs(
    detinfo::DetectorClocksData const& clockData,
    recob::Hit const& hit) const
  {
    return BackTracker::HitToTrackIDEs(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<sim::TrackIDE> BackTrackerService::HitToTrackIDEs(
    detinfo::DetectorClocksData const& clockData,
    art::Ptr<recob::Hit> const& hit) const
  {
    return BackTracker::HitToTrackIDEs(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<int> BackTrackerService::HitToTrackIds(detinfo::DetectorClocksData const& clockData,
                                                     recob::Hit const& hit) const
  {
    return BackTracker::HitToTrackIds(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<sim::TrackIDE> BackTrackerService::HitToEveTrackIDEs(
    detinfo::DetectorClocksData const& clockData,
    recob::Hit const& hit) const
  {
    return BackTracker::HitToEveTrackIDEs(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<sim::TrackIDE> BackTrackerService::HitToEveTrackIDEs(
    detinfo::DetectorClocksData const& clockData,
    art::Ptr<recob::Hit> const& hit) const
  {
    return BackTracker::HitToEveTrackIDEs(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<art::Ptr<recob::Hit>> BackTrackerService::TrackIdToHits_Ps(
    detinfo::DetectorClocksData const& clockData,
    const int tkId,
    std::vector<art::Ptr<recob::Hit>> const& hitsIn) const
  {
    return BackTracker::TrackIdToHits_Ps(clockData, tkId, hitsIn);
  }

  //---------------------------------------------------------------------
  std::vector<std::vector<art::Ptr<recob::Hit>>> BackTrackerService::TrackIdsToHits_Ps(
    detinfo::DetectorClocksData const& clockData,
    std::vector<int> const& tkIds,
    std::vector<art::Ptr<recob::Hit>> const& hitsIn) const
  {
    return BackTracker::TrackIdsToHits_Ps(clockData, tkIds, hitsIn);
  }

  //---------------------------------------------------------------------
  std::vector<sim::IDE> BackTrackerService::HitToAvgSimIDEs(
    detinfo::DetectorClocksData const& clockData,
    recob::Hit const& hit) const
  {
    return BackTracker::HitToAvgSimIDEs(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<sim::IDE> BackTrackerService::HitToAvgSimIDEs(
    detinfo::DetectorClocksData const& clockData,
    art::Ptr<recob::Hit> hit) const
  {
    return BackTracker::HitToAvgSimIDEs(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<const sim::IDE*> BackTrackerService::HitToSimIDEs_Ps(
    detinfo::DetectorClocksData const& clockData,
    recob::Hit const& hit) const
  {
    return BackTracker::HitToSimIDEs_Ps(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<const sim::IDE*> BackTrackerService::HitToSimIDEs_Ps(
    detinfo::DetectorClocksData const& clockData,
    art::Ptr<recob::Hit> const& hit) const
  {
    return BackTracker::HitToSimIDEs_Ps(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::SimIDEsToXYZ(std::vector<sim::IDE> const& ides) const
  {
    return BackTracker::SimIDEsToXYZ(ides);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::SimIDEsToXYZ(
    std::vector<const sim::IDE*> const& ide_Ps) const
  {
    return BackTracker::SimIDEsToXYZ(ide_Ps);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                                   const recob::Hit& hit) const
  {
    return BackTracker::HitToXYZ(clockData, hit);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                                   art::Ptr<recob::Hit> const& hit) const
  {
    return HitToXYZ(clockData, *hit);
  }

  //---------------------------------------------------------------------
  double BackTrackerService::HitCollectionPurity(
    detinfo::DetectorClocksData const& clockData,
    std::set<int> const& trackIds,
    std::vector<art::Ptr<recob::Hit>> const& hits) const
  {
    return BackTracker::HitCollectionPurity(clockData, trackIds, hits);
  }

  //---------------------------------------------------------------------
  double BackTrackerService::HitChargeCollectionPurity(
    detinfo::DetectorClocksData const& clockData,
    std::set<int> const& trackIds,
    std::vector<art::Ptr<recob::Hit>> const& hits) const
  {
    return BackTracker::HitChargeCollectionPurity(clockData, trackIds, hits);
  }

  //---------------------------------------------------------------------
  double BackTrackerService::HitCollectionEfficiency(
    detinfo::DetectorClocksData const& clockData,
    std::set<int> const& trackIds,
    std::vector<art::Ptr<recob::Hit>> const& hits,
    std::vector<art::Ptr<recob::Hit>> const& allhits,
    geo::View_t const& view) const
  {
    return BackTracker::HitCollectionEfficiency(clockData, trackIds, hits, allhits, view);
  }

  //---------------------------------------------------------------------
  double BackTrackerService::HitChargeCollectionEfficiency(
    detinfo::DetectorClocksData const& clockData,
    std::set<int> const& trackIds,
    std::vector<art::Ptr<recob::Hit>> const& hits,
    std::vector<art::Ptr<recob::Hit>> const& allhits,
    geo::View_t const& view) const
  {
    return BackTracker::HitChargeCollectionEfficiency(clockData, trackIds, hits, allhits, view);
  }

  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfTrackIds() const
  {
    return BackTracker::GetSetOfTrackIds();
  }

  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfEveIds() const
  {
    return BackTracker::GetSetOfEveIds();
  }

  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfTrackIds(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::Hit>> const& hits) const
  {
    return BackTracker::GetSetOfTrackIds(clockData, hits);
  }

  //---------------------------------------------------------------------
  std::set<int> BackTrackerService::GetSetOfEveIds(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::Hit>> const& hits) const
  {
    return BackTracker::GetSetOfEveIds(clockData, hits);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::SpacePointHitsToWeightedXYZ(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::Hit>> const& hits) const
  {
    return BackTracker::SpacePointHitsToWeightedXYZ(clockData, hits);
  }

  //---------------------------------------------------------------------
  std::vector<art::Ptr<recob::Hit>> BackTrackerService::SpacePointToHits_Ps(
    art::Ptr<recob::SpacePoint> const& spt) const
  {
    // Note, this function is goofy in that it doesn't use the hits prepared by
    // the module. That should likely be changed to make things more uniform.
    if (!fEvt) {
      throw cet::exception("BackTrackerService")
        << "This function is not yet implimented pending the implimentation of "
           "backtracker lazy loading.";
    }
    return BackTracker::SpacePointToHits_Ps(spt, *fEvt);
  }

  //---------------------------------------------------------------------
  std::vector<double> BackTrackerService::SpacePointToXYZ(
    detinfo::DetectorClocksData const& clockData,
    art::Ptr<recob::SpacePoint> const& spt) const
  {
    // if( !this->priv_CanRun(*fEvt)) { this->priv_PrepFailed(); }
    if (!fEvt) {
      throw cet::exception("BackTrackerService")
        << "This function is not yet implimented pending the implimentation of "
           "backtracker lazy loading.";
    }
    return BackTracker::SpacePointToXYZ(clockData, spt, *fEvt);
  }

} // end namespace cheat
