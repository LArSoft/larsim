////////////////////////////////////////////////////////////////////////
//
// \file: PhotonBackTrackerServiceService_service.cc
//
//jason.stock@mines.sdsmt.edu
//Based on the BackTracker_service by Brian Rebel
//
////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_PHOTONBACKTRACKERSERVICESERVICE_H
#define CHEAT_PHOTONBACKTRACKERSERVICESERVICE_H

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/MCCheater/PhotonBackTracker.h"

namespace art {
  class ActivityRegistry;
  class Event;
}

#include "canvas/Persistency/Common/Ptr.h"

#include "fhiclcpp/fwd.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"

#include <set>
#include <vector>

namespace cheat {
  class PhotonBackTrackerService : private PhotonBackTracker {
  public:
    struct fhiclConfig {
      fhicl::Table<PhotonBackTracker::fhiclConfig> PhotonBackTrackerTable{
        fhicl::Name("PhotonBackTracker"),
        fhicl::Comment(
          "This if the fhicl configuration of the PhotonBackTracker service provider.")};
    };

    using provider_type = PhotonBackTracker;
    const provider_type* provider() const { return static_cast<provider_type const*>(this); }

    PhotonBackTrackerService(fhicl::ParameterSet const& pSet, art::ActivityRegistry& reg);
    PhotonBackTrackerService(fhiclConfig const& config, art::ActivityRegistry& reg);

    // Temporarily include a rebuild function until the lazy rebuild works.

    void Rebuild(art::Event const& evt);

    //----------------------------------------------------------------------
    // BackTracking implementation functions
    std::vector<art::Ptr<sim::OpDetBacktrackerRecord>> const& OpDetBTRs();
    double GetDelay();
    std::vector<const sim::SDP*> TrackIdToSimSDPs_Ps(int id);
    art::Ptr<sim::OpDetBacktrackerRecord> FindOpDetBTR(int opDetNum);
    std::vector<sim::TrackSDP> OpDetToTrackSDPs(int OpDetNum,
                                                double opHit_start_time,
                                                double opHit_end_time);
    std::vector<sim::TrackSDP> OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P);
    std::vector<sim::TrackSDP> OpHitToTrackSDPs(recob::OpHit const& opHit);
    std::vector<int> OpHitToTrackIds(recob::OpHit const& opHit);
    std::vector<int> OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit_P);
    std::vector<int> OpHitToEveTrackIds(recob::OpHit const& opHit);
    std::vector<int> OpHitToEveTrackIds(art::Ptr<recob::OpHit> const& opHit_P);
    std::vector<sim::TrackSDP> OpHitToEveTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P);
    std::vector<sim::TrackSDP> OpHitToEveTrackSDPs(recob::OpHit const& opHit);
    std::vector<art::Ptr<recob::OpHit>> TrackIdToOpHits_Ps(
      int tkId,
      std::vector<art::Ptr<recob::OpHit>> const& hitsIn);
    std::vector<std::vector<art::Ptr<recob::OpHit>>> TrackIdsToOpHits_Ps(
      std::vector<int> const& tkIds,
      std::vector<art::Ptr<recob::OpHit>> const& hitsIn);
    std::vector<const sim::SDP*> OpHitToSimSDPs_Ps(recob::OpHit const& opHit);
    std::vector<const sim::SDP*> OpHitToSimSDPs_Ps(art::Ptr<recob::OpHit> const& opHit_P);
    std::unordered_set<const sim::SDP*> OpHitToEveSimSDPs_Ps(recob::OpHit const& opHit);
    std::unordered_set<const sim::SDP*> OpHitToEveSimSDPs_Ps(art::Ptr<recob::OpHit>& opHit_P);
    std::vector<double> SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&;
    std::vector<double> SimSDPsToXYZ(std::vector<const sim::SDP*> const& sdps_Ps);
    std::vector<double> OpHitToXYZ(recob::OpHit const& opHit);
    std::vector<double> OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit_P);
    std::set<int> GetSetOfEveIds();
    std::set<int> GetSetOfTrackIds();
    std::set<int> GetSetOfEveIds(std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps);
    std::set<int> GetSetOfEveIds(const std::vector<recob::OpHit>& opHits);
    std::set<int> GetSetOfTrackIds(std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps);
    std::set<int> GetSetOfTrackIds(std::vector<recob::OpHit> const& opHits);
    double OpHitCollectionPurity(std::set<int> const& tkIds,
                                 std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps);
    double OpHitLightCollectionPurity(std::set<int> const& tkIds,
                                      std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps);
    double OpHitCollectionEfficiency(std::set<int> const& tkIds,
                                     std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps,
                                     std::vector<art::Ptr<recob::OpHit>> const& opHitsIn_Ps);
    double OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
                                          std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps,
                                          std::vector<art::Ptr<recob::OpHit>> const& opHitsIn_Ps);
    double OpHitChargeCollectionEfficiency(std::set<int> const& tkIds,
                                           std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps,
                                           std::vector<art::Ptr<recob::OpHit>> const& opHitsIn_Ps)
    {
      return PhotonBackTracker::OpHitLightCollectionEfficiency(tkIds, opHits_Ps, opHitsIn_Ps);
    } //Exists only temporarily. Is deprecated.
    std::set<int> OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P) const;
    std::vector<art::Ptr<recob::OpHit>> OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P);
    std::vector<double> OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P);

  private:
    //The PhotonBackTrackerService has no parameters.

    void priv_PrepEvent(art::Event const& evt, art::ScheduleContext);
    void priv_PrepFailed();
    void priv_PrepOpDetBTRs(art::Event const& evt);
    void priv_PrepOpFlashToOpHits(art::Event const& evt);

    bool priv_CanRun(art::Event const& evt);
    bool priv_OpDetBTRsReady() { return PhotonBackTracker::BTRsReady(); }
    bool priv_OpFlashToOpHitsReady() { return PhotonBackTracker::OpFlashToOpHitsReady(); }

  }; //Class PhotonBackTrackerService

  //----------------------------------------------------------------------
} // namespace

DECLARE_ART_SERVICE(cheat::PhotonBackTrackerService, LEGACY)

#endif //CHEAT_PHOTONBACKTRACKERSERVICESERVICE_H
