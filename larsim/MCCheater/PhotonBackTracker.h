////////////////////////////////////////////////////////////////////
//
// \file PhotonBackTracker.h
// \brief The functions needed for the PhotonBackTracker class needed by the PhotonBackTrackerService in order to connect truth information with reconstruction.
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
///////////////////////////////////////////////////////////////////
//
//TODO: Impliment alternate backtracking scheme developed by T. Usher
//
///////////////////////////////////////////////////////////////////
#ifndef CHEAT_PHOTONBACKTRACKER_H
#define CHEAT_PHOTONBACKTRACKER_H

//CPP
#include <map>
#include <unordered_set>
#include <vector>

//Framework
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/fwd.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

//LArSoft
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"

namespace cheat {
  class ParticleInventory;
}
namespace sim {
  struct SDP;
  struct TrackSDP;
}

namespace cheat {
  class PhotonBackTracker {
  public:
    struct fhiclConfig {
      fhicl::Atom<double> Delay{
        fhicl::Name("Delay"),
        fhicl::Comment("The delay time needed to correctly account for the optical simulation and "
                       "optical systems simulation. (The time between when a g4partcile was made, "
                       "and when the simulation write out says a signal was recorded)."),
        0};
      fhicl::Atom<art::InputTag> G4ModuleLabel{
        fhicl::Name("G4ModuleLabel"),
        fhicl::Comment(
          "The label of the LArG4 module used to produce the art file we will be using."),
        "largeant"};
      fhicl::Sequence<art::InputTag> G4ModuleLabels{
        fhicl::Name("G4ModuleLabels"),
        fhicl::Comment(
          "The labels of the LArG4 modules used to produce the art file we will be using."),
        fhicl::Sequence<art::InputTag>::default_type{}};
      // The following should be removed and replaced with some way to access the
      // OpHitLabel given by the user in their own analysis module to avoid differing
      // definitions.
      fhicl::Atom<art::InputTag> OpHitLabel{
        fhicl::Name("OpHitLabel"),
        fhicl::Comment("The default label for the module to use when grabbing OpHits"),
        "ophit"};
      // The following should be removed and replaced with some way to access the
      // OpFlashLabel given by the user in their own analysis module to avoid differing
      // definitions.
      fhicl::Atom<art::InputTag> OpFlashLabel{
        fhicl::Name("OpFlashLabel"),
        fhicl::Comment("The default label for the module to use when grabbing OpFlash"),
        "opflash"};
      fhicl::Atom<double> MinOpHitEnergyFraction{
        fhicl::Name("MinOpHitEnergyFraction"),
        fhicl::Comment("The minimum contribution an energy deposit must make to a Hit to be "
                       "considered part of that hit."),
        0.010};
    };

    //----------------------------------------------------------------
    PhotonBackTracker(fhiclConfig const& config,
                      const cheat::ParticleInventory* partInv,
                      geo::WireReadoutGeom const* wireReadoutGeom);

    //----------------------------------------------------------------
    PhotonBackTracker(fhicl::ParameterSet const& pSet,
                      const cheat::ParticleInventory* partInv,
                      geo::WireReadoutGeom const* wireReadoutGeom);

    //----------------------------------------------------------------
    PhotonBackTracker(PhotonBackTracker const&) = delete;

    //-----------------------------------------------------
    template <typename Evt>
    void PrepEvent(Evt const& evt);

    //-----------------------------------------------------
    template <typename Evt>
    bool CanRun(Evt const& evt);

    //-----------------------------------------------------
    template <typename Evt>
    void PrepOpDetBTRs(Evt const& evt);

    //-----------------------------------------------------
    template <typename Evt>
    void PrepOpFlashToOpHits(Evt const& evt);

    //-----------------------------------------------------
    std::vector<art::Ptr<recob::OpHit>> OpFlashToOpHits_Ps(art::Ptr<recob::OpFlash>& flash_P) const;

    //-----------------------------------------------------
    std::vector<double> OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P) const;

    //-----------------------------------------------------
    std::set<int> OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P) const;

    double GetDelay() const;
    void ClearEvent();
    bool BTRsReady() const;
    bool OpFlashToOpHitsReady() const;

    std::vector<art::Ptr<sim::OpDetBacktrackerRecord>> const& OpDetBTRs() const;
    std::vector<const sim::SDP*> TrackIdToSimSDPs_Ps(int id) const;
    art::Ptr<sim::OpDetBacktrackerRecord> FindOpDetBTR(int opDetNum) const;

    //-----------------------------------------------------

    std::vector<sim::TrackSDP> OpDetToTrackSDPs(int OpDetNum,
                                                double opHit_start_time,
                                                double opHit_end_time) const;
    std::vector<sim::TrackSDP> OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P) const;
    std::vector<sim::TrackSDP> OpHitToTrackSDPs(recob::OpHit const& opHit) const;
    std::vector<int> OpHitToTrackIds(recob::OpHit const& opHit) const;
    std::vector<int> OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit) const;
    std::vector<int> OpHitToEveTrackIds(recob::OpHit const& opHit);
    std::vector<int> OpHitToEveTrackIds(art::Ptr<recob::OpHit> const& opHit_P);
    std::vector<sim::TrackSDP> OpHitToEveTrackSDPs(art::Ptr<recob::OpHit> const& opHit_P) const;
    std::vector<sim::TrackSDP> OpHitToEveTrackSDPs(recob::OpHit const& opHit) const;

    //-----------------------------------------------------
    std::vector<art::Ptr<recob::OpHit>> TrackIdToOpHits_Ps(
      int tkId,
      std::vector<art::Ptr<recob::OpHit>> const& hitsIn) const;
    std::vector<std::vector<art::Ptr<recob::OpHit>>> TrackIdsToOpHits_Ps(
      std::vector<int> const& tkIds,
      std::vector<art::Ptr<recob::OpHit>> const& hitsIn) const;

    //-----------------------------------------------------
    std::vector<const sim::SDP*> OpHitToSimSDPs_Ps(recob::OpHit const& opHit) const;
    std::vector<const sim::SDP*> OpHitToSimSDPs_Ps(art::Ptr<recob::OpHit> const& opHit_P) const;
    std::vector<const sim::SDP*> OpHitsToSimSDPs_Ps(
      std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps) const;

    //-----------------------------------------------------
    std::vector<double> SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&;
    std::vector<double> SimSDPsToXYZ(std::vector<const sim::SDP*> const& sdps_Ps) const&;

    //-----------------------------------------------------
    std::vector<double> OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit);
    std::vector<double> OpHitToXYZ(recob::OpHit const& opHit);
    std::vector<double> OpHitsToXYZ(std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps) const;

    //-----------------------------------------------------
    std::unordered_set<const sim::SDP*> OpHitToEveSimSDPs_Ps(recob::OpHit const& opHit);
    std::unordered_set<const sim::SDP*> OpHitToEveSimSDPs_Ps(art::Ptr<recob::OpHit>& opHit_P);

    //-----------------------------------------------------
    std::set<int> GetSetOfEveIds() const;
    std::set<int> GetSetOfTrackIds() const;
    std::set<int> GetSetOfEveIds(std::vector<art::Ptr<recob::OpHit>> const& opHits) const;
    std::set<int> GetSetOfEveIds(std::vector<recob::OpHit> const& opHits) const;
    std::set<int> GetSetOfTrackIds(std::vector<art::Ptr<recob::OpHit>> const& opHits) const;
    std::set<int> GetSetOfTrackIds(std::vector<recob::OpHit> const& opHits) const;

    //-----------------------------------------------------
    double OpHitCollectionPurity(std::set<int> const& tkIds,
                                 std::vector<art::Ptr<recob::OpHit>> const& opHits);
    double OpHitLightCollectionPurity(std::set<int> const& tkIds,
                                      std::vector<art::Ptr<recob::OpHit>> const& opHits);
    double OpHitCollectionEfficiency(std::set<int> const& tkIds,
                                     std::vector<art::Ptr<recob::OpHit>> const& opHits,
                                     std::vector<art::Ptr<recob::OpHit>> const& opHitsIn);
    double OpHitLightCollectionEfficiency(std::set<int> const& tkIds,
                                          std::vector<art::Ptr<recob::OpHit>> const& opHits,
                                          std::vector<art::Ptr<recob::OpHit>> const& opHitsIn);

  private:
    const cheat::ParticleInventory* fPartInv; //The constructor needs to put something in here
    const geo::WireReadoutGeom* fWireReadoutGeom;
    const double fDelay;
    const art::InputTag fG4ModuleLabel;
    const std::vector<art::InputTag> fG4ModuleLabels;
    const art::InputTag fOpHitLabel;
    const art::InputTag fOpFlashLabel;
    const double fMinOpHitEnergyFraction;
    mutable std::vector<art::Ptr<sim::OpDetBacktrackerRecord>> priv_OpDetBTRs;
    std::map<art::Ptr<recob::OpFlash>, std::vector<art::Ptr<recob::OpHit>>> priv_OpFlashToOpHits;
  };
}

#include "PhotonBackTracker.tcc"

#endif
