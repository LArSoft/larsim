/////////////////////////////////////////////////////////////////////////////////
// \file BackTrackerService.h
// \A service for backtracking reconstruction information to its truth information
//
// \author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov
//
/////////////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_BACKTRACKERSERVICESERVICE_H
#define CHEAT_BACKTRACKERSERVICESERVICE_H

#include "larsim/MCCheater/BackTracker.h"

namespace detinfo {
  class DetectorClocksData;
}

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Persistency/Provenance/ScheduleContext.h"

namespace art {
  class ActivityRegistry;
  class Event;
}

#include "canvas/Persistency/Common/Ptr.h"

#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"

namespace fhicl {
  class ParameterSet;
}

#include <set>
#include <vector>

namespace cheat {
  class BackTrackerService : private BackTracker {
  public:
    struct fhiclConfig {
      fhicl::Table<BackTracker::fhiclConfig> BackTrackerTable{
        fhicl::Name("BackTracker"),
        fhicl::Comment("This is the fhicl configuration of the BackTracker service provider.")};
    };

    using provider_type = BackTracker;
    const provider_type* provider() const { return static_cast<provider_type const*>(this); }

    BackTrackerService(const fhicl::ParameterSet& pSet, art::ActivityRegistry& reg);
    BackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg);

    //Temporrary rebuild function (until I can make this lazy again.

    void Rebuild(const art::Event& evt);

    const std::vector<art::Ptr<sim::SimChannel>>& SimChannels() const;

    std::vector<const sim::IDE*> TrackIdToSimIDEs_Ps(int const& id) const;
    std::vector<const sim::IDE*> TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const;

    art::Ptr<sim::SimChannel> FindSimChannel(raw::ChannelID_t channel) const;

    std::vector<sim::TrackIDE> ChannelToTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                                  raw::ChannelID_t channel,
                                                  const double hit_start_time,
                                                  const double hit_end_time) const;

    std::vector<sim::TrackIDE> HitToTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                              recob::Hit const& hit) const;
    std::vector<sim::TrackIDE> HitToTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                              art::Ptr<recob::Hit> const& hit) const;

    std::vector<int> HitToTrackIds(detinfo::DetectorClocksData const& clockData,
                                   recob::Hit const& hit) const;

    std::vector<sim::TrackIDE> HitToEveTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                                 recob::Hit const& hit) const;
    std::vector<sim::TrackIDE> HitToEveTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                                 art::Ptr<recob::Hit> const& hit) const;

    std::vector<art::Ptr<recob::Hit>> TrackIdToHits_Ps(
      detinfo::DetectorClocksData const& clockData,
      int tkId,
      std::vector<art::Ptr<recob::Hit>> const& hitsIn) const;
    //    std::vector<art::Ptr<recob::Hit>> TrackIdToHits_Ps(const int& tkId)
    //    const;

    std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIdsToHits_Ps(
      detinfo::DetectorClocksData const& clockData,
      std::vector<int> const& tkIds,
      std::vector<art::Ptr<recob::Hit>> const& hitsIn) const;

    std::vector<sim::IDE> HitToAvgSimIDEs(detinfo::DetectorClocksData const& clockData,
                                          recob::Hit const& hit) const;
    std::vector<sim::IDE> HitToAvgSimIDEs(detinfo::DetectorClocksData const& clockData,
                                          art::Ptr<recob::Hit> hit) const;

    std::vector<const sim::IDE*> HitToSimIDEs_Ps(detinfo::DetectorClocksData const& clockData,
                                                 recob::Hit const& hit) const;
    std::vector<const sim::IDE*> HitToSimIDEs_Ps(detinfo::DetectorClocksData const& clockData,
                                                 art::Ptr<recob::Hit> const& hit) const;

    std::vector<double> SimIDEsToXYZ(std::vector<sim::IDE> const& ides) const;
    std::vector<double> SimIDEsToXYZ(std::vector<const sim::IDE*> const& ide_Ps) const;

    std::vector<double> HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                 const recob::Hit& hit) const;
    std::vector<double> HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                 art::Ptr<recob::Hit> const& hit) const;

    double HitCollectionPurity(detinfo::DetectorClocksData const& clockData,
                               std::set<int> const& trackIds,
                               std::vector<art::Ptr<recob::Hit>> const& hits) const;
    double HitChargeCollectionPurity(detinfo::DetectorClocksData const& clockData,
                                     std::set<int> const& trackIds,
                                     std::vector<art::Ptr<recob::Hit>> const& hits) const;

    double HitCollectionEfficiency(detinfo::DetectorClocksData const& clockData,
                                   std::set<int> const& trackIds,
                                   std::vector<art::Ptr<recob::Hit>> const& hits,
                                   std::vector<art::Ptr<recob::Hit>> const& allhits,
                                   geo::View_t const& view) const;
    double HitChargeCollectionEfficiency(detinfo::DetectorClocksData const& clockData,
                                         std::set<int> const& trackIds,
                                         std::vector<art::Ptr<recob::Hit>> const& hits,
                                         std::vector<art::Ptr<recob::Hit>> const& allhits,
                                         geo::View_t const& view) const;

    std::set<int> GetSetOfTrackIds() const;
    std::set<int> GetSetOfEveIds() const;

    std::set<int> GetSetOfTrackIds(detinfo::DetectorClocksData const& clockData,
                                   std::vector<art::Ptr<recob::Hit>> const& hits) const;
    std::set<int> GetSetOfEveIds(detinfo::DetectorClocksData const& clockData,
                                 std::vector<art::Ptr<recob::Hit>> const& hits) const;

    std::vector<double> SpacePointHitsToWeightedXYZ(
      detinfo::DetectorClocksData const& clockData,
      std::vector<art::Ptr<recob::Hit>> const& hits) const;
    std::vector<art::Ptr<recob::Hit>> SpacePointToHits_Ps(
      art::Ptr<recob::SpacePoint> const& spt) const;
    std::vector<double> SpacePointToXYZ(detinfo::DetectorClocksData const& clockData,
                                        art::Ptr<recob::SpacePoint> const& spt) const;

  private:
    const art::Event* fEvt = nullptr;

    // Prep functions go here.
    void priv_PrepEvent(const art::Event& evt, art::ScheduleContext);
    void priv_PrepSimChannels(const art::Event& evt);
    //      void priv_PrepAllHitList ();
    void priv_PrepFailed();

    bool priv_CanRun(const art::Event& evt);

    bool priv_SimChannelsReady() { return BackTracker::SimChannelsReady(); }
    //      bool priv_AllHitListReady() { return
    //      BackTracker::AllHitListReady();}

  }; // class BackTrackerService

} // end namespace cheat
DECLARE_ART_SERVICE(cheat::BackTrackerService, LEGACY)

#endif // CHEAT_BACKTRACKERSERVICESERVICE_H
