////////////////////////////////////////////////////////////////////////////
//
// \file BackTracker.h
// \brief Functions needed by the BackTracker service in order to connect truth information with reconstruction.
//
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_BACKTRACKER_H
#define CHEAT_BACKTRACKER_H

#include <vector>

#include "fhiclcpp/fwd.h"
#include "fhiclcpp/types/Atom.h"

#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventory.h"

namespace detinfo {
  class DetectorClocksData;
}
namespace recob {
  class SpacePoint;
}

namespace cheat {

  class BackTracker {
  public:
    // Structure for configuration parameters. (fhicl validation)
    struct fhiclConfig {
      fhicl::Atom<art::InputTag> G4ModuleLabel{
        fhicl::Name("G4ModuleLabel"),
        fhicl::Comment("The label of the LArG4   module used to produce the "
                       "art file we will be using."),
        "largeant"};
      fhicl::Atom<art::InputTag> SimChannelModuleLabel{
        fhicl::Name("SimChannelModuleLabel"),
        fhicl::Comment("The label of the module containing the sim::SimChannel product."),
        G4ModuleLabel()}; // -- D.R. label not required, if not provided
                          // defaults to the value of G4ModuleLabel
      fhicl::Atom<art::InputTag> DefaultHitModuleLabel{
        fhicl::Name("DefaultHitModuleLabel"),
        fhicl::Comment("The label  of the module used to produce the hits in the art file "
                       "we will default to when no hitlist is provided."),
        "hitfd"};
      fhicl::Atom<double> MinHitEnergyFraction{
        fhicl::Name("MinHitEnergyFraction"),
        fhicl::Comment("The minimum     contribution an energy deposit must "
                       "make to a Hit to be considered part of that hit."),
        0.010};
      fhicl::Atom<bool> OverrideRealData{
        fhicl::Name("OverrideRealData"),
        fhicl::Comment("Option when overlaying simulation on real data, to tell the "
                       "backtracker to continue even if event looks like data."),
        false};
      fhicl::Atom<double> HitTimeRMS{fhicl::Name("HitTimeRMS"),
                                     fhicl::Comment("The number of RMS units to move away"
                                                    "From a Hit peak time for searching IDE."),
                                     1.0};
    };

    BackTracker(const fhiclConfig& config,
                const cheat::ParticleInventory* partInv,
                geo::GeometryCore const* geom,
                geo::WireReadoutGeom const* wireReadoutGeom);
    BackTracker(const fhicl::ParameterSet& pSet,
                const cheat::ParticleInventory* partInv,
                geo::GeometryCore const* geom,
                geo::WireReadoutGeom const* wireReadoutGeom);

    // I may need to include this to delete copy of service providers.
    BackTracker(BackTracker const&) = delete;

    template <typename Evt>
    void PrepEvent(const Evt& evt);

    template <typename Evt>
    void PrepSimChannels(const Evt& evt);

    //-----------------------------------------------------
    template <typename Evt>
    bool CanRun(const Evt& evt)
    {
      return !evt.isRealData() || fOverrideRealData;
    }

    //-----------------------------------------------------
    template <typename Evt>
    std::vector<art::Ptr<recob::Hit>> SpacePointToHits_Ps(art::Ptr<recob::SpacePoint> const& spt,
                                                          const Evt& evt) const;

    //-----------------------------------------------------
    template <typename Evt>
    std::vector<double> SpacePointToXYZ(detinfo::DetectorClocksData const& clockData,
                                        art::Ptr<recob::SpacePoint> const& spt,
                                        const Evt& evt) const;

    //-----------------------------------------------------
    void ClearEvent();

    bool SimChannelsReady() const { return !fSimChannels.empty(); }

    const std::vector<art::Ptr<sim::SimChannel>>& SimChannels() const { return fSimChannels; }
    // All Hit List would go here. We explicitly choose not to include it, as the user
    // should not be using backtracker to access Hits. This could change in a concievable
    // future use case where we also allow the user to define what the "AllHitList" should
    // be, though this would have ramifications on other functions.

    std::vector<const sim::IDE*> TrackIdToSimIDEs_Ps(int const& id) const;
    std::vector<const sim::IDE*> TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const;

    /**
     * @brief Returns the cached `sim::SimChannel` on the specified `channel`.
     * @param channel ID of the TPC channel to find
     * @return _art_ pointer to `sim::SimChannel`, or an null pointer if none
     * @see FindSimChannel()
     */
    art::Ptr<sim::SimChannel> FindSimChannelPtr(raw::ChannelID_t channel) const;

    /**
     * @brief Returns the cached `sim::SimChannel` on the specified `channel`.
     * @param channel ID of the TPC channel to find
     * @return _art_ pointer to `sim::SimChannel`
     * @throw cet::exception (category: `"BackTracker"`) if no `sim::SimChannel`
     *                       for the requested  `channel` found
     * @see FindSimChannelPtr()
     */
    art::Ptr<sim::SimChannel> FindSimChannel(raw::ChannelID_t channel) const;

    std::vector<sim::TrackIDE> ChannelToTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                                  raw::ChannelID_t channel,
                                                  const double hit_start_time,
                                                  const double hit_end_time) const;

    // Track IDEs cannot be returned as pointers, as they dont exist in the data product,
    // and we will not be storing them.
    std::vector<sim::TrackIDE> HitToTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                              recob::Hit const& hit) const;
    std::vector<sim::TrackIDE> HitToTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                              art::Ptr<recob::Hit> const& hit) const
    {
      return HitToTrackIDEs(clockData, *hit);
    }

    std::vector<int> HitToTrackIds(detinfo::DetectorClocksData const& clockData,
                                   recob::Hit const& hit) const;

    std::vector<sim::TrackIDE> HitToEveTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                                 recob::Hit const& hit) const;
    std::vector<sim::TrackIDE> HitToEveTrackIDEs(detinfo::DetectorClocksData const& clockData,
                                                 art::Ptr<recob::Hit> const& hit) const
    {
      return HitToEveTrackIDEs(clockData, *hit);
    }

    // I will not return these by copy, as that could get very large very quickly.
    std::vector<art::Ptr<recob::Hit>> TrackIdToHits_Ps(
      detinfo::DetectorClocksData const& clockData,
      int tkId,
      std::vector<art::Ptr<recob::Hit>> const& hitsIn) const;

    std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIdsToHits_Ps(
      detinfo::DetectorClocksData const& clockData,
      std::vector<int> const& tkIds,
      std::vector<art::Ptr<recob::Hit>> const& hitsIn) const;

    std::vector<sim::IDE> HitToAvgSimIDEs(detinfo::DetectorClocksData const& clockData,
                                          recob::Hit const& hit) const;
    std::vector<sim::IDE> HitToAvgSimIDEs(detinfo::DetectorClocksData const& clockData,
                                          art::Ptr<recob::Hit> const& hit) const
    {
      return HitToAvgSimIDEs(clockData, *hit);
    }

    std::vector<const sim::IDE*> HitToSimIDEs_Ps(detinfo::DetectorClocksData const& clockData,
                                                 recob::Hit const& hit) const;
    std::vector<const sim::IDE*> HitToSimIDEs_Ps(detinfo::DetectorClocksData const& clockData,
                                                 art::Ptr<recob::Hit> const& hit) const
    {
      return HitToSimIDEs_Ps(clockData, *hit);
    }

    std::vector<double> SimIDEsToXYZ(std::vector<sim::IDE> const& ides) const;
    std::vector<double> SimIDEsToXYZ(std::vector<const sim::IDE*> const& ide_Ps) const;

    std::vector<double> HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                 const recob::Hit& hit) const;
    std::vector<double> HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                 art::Ptr<recob::Hit> const& hit) const
    {
      return HitToXYZ(clockData, *hit);
    }

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

    std::set<int> GetSetOfTrackIds() const { return fPartInv->GetSetOfTrackIds(); }
    std::set<int> GetSetOfEveIds() const { return fPartInv->GetSetOfEveIds(); }

    std::set<int> GetSetOfTrackIds(detinfo::DetectorClocksData const& clockData,
                                   std::vector<art::Ptr<recob::Hit>> const& hits) const;
    std::set<int> GetSetOfEveIds(detinfo::DetectorClocksData const& clockData,
                                 std::vector<art::Ptr<recob::Hit>> const& hits) const;

    std::vector<double> SpacePointHitsToWeightedXYZ(
      detinfo::DetectorClocksData const& clockData,
      std::vector<art::Ptr<recob::Hit>> const& hits) const;

  private:
    const cheat::ParticleInventory* fPartInv; // The constructor needs to put something in here
    const geo::GeometryCore* fGeom;
    const geo::WireReadoutGeom* fWireReadoutGeom;
    const art::InputTag fG4ModuleLabel;
    const art::InputTag fSimChannelModuleLabel;
    const art::InputTag fHitLabel;
    const double fMinHitEnergyFraction;
    const bool fOverrideRealData;
    const double fHitTimeRMS;

    mutable std::vector<art::Ptr<sim::SimChannel>> fSimChannels;

  }; // end class BackTracker

} // end namespace cheat

#include "BackTracker.tcc"

#endif // CHEAT_BACKTRACKER_H
