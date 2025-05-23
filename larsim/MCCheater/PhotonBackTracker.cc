////////////////////////////////////////////////////////////////////////////////
//
// \file PhotonBackTracker.cc
// \brief The functions needed for the PhotonBackTracker class needed by the
//   PhotonBackTrackerService in order to connect truth information with
//   reconstruction.
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////////
//
//TODO: Impliment alternate backtracking scheme developed by T. Usher
//TODO: OpChanToOpDetSDPs (Expanded Clone of OpDetNumToOpDetSDPs
//
////////////////////////////////////////////////////////////////////////////////

// LArSoft
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/MCCheater/ParticleInventory.h"

// Framework
#include "messagefacility/MessageLogger/MessageLogger.h"

// STL
#include <map>

namespace cheat {

  //----------------------------------------------------------------
  PhotonBackTracker::PhotonBackTracker(fhiclConfig const& config,
                                       cheat::ParticleInventory const* partInv,
                                       geo::WireReadoutGeom const* wireReadoutGeom)
    : fPartInv(partInv)
    , fWireReadoutGeom(wireReadoutGeom)
    , fDelay(config.Delay())
    , fG4ModuleLabel(config.G4ModuleLabel())
    , fG4ModuleLabels(config.G4ModuleLabels())
    , fOpHitLabel(config.OpHitLabel())
    , fOpFlashLabel(config.OpFlashLabel())
    , fMinOpHitEnergyFraction(config.MinOpHitEnergyFraction())
  {}

  //----------------------------------------------------------------
  PhotonBackTracker::PhotonBackTracker(fhicl::ParameterSet const& pSet,
                                       cheat::ParticleInventory const* partInv,
                                       geo::WireReadoutGeom const* wireReadoutGeom)
    : fPartInv(partInv)
    , fWireReadoutGeom(wireReadoutGeom)
    , fDelay(pSet.get<double>("Delay"))
    , fG4ModuleLabel(pSet.get<art::InputTag>("G4ModuleLabel", "largeant"))
    , fG4ModuleLabels(pSet.get<std::vector<art::InputTag>>("G4ModuleLabels", {}))
    , fOpHitLabel(pSet.get<art::InputTag>("OpHitLabel", "ophit"))
    , fOpFlashLabel(pSet.get<art::InputTag>("OpFlashLabel", "opflash"))
    , fMinOpHitEnergyFraction(pSet.get<double>("MinimumOpHitEnergyFraction", 0.1))
  {}

  //----------------------------------------------------------------
  double PhotonBackTracker::GetDelay() const
  {
    return fDelay;
  }

  //----------------------------------------------------------------
  void PhotonBackTracker::ClearEvent()
  {
    priv_OpDetBTRs.clear();
    priv_OpFlashToOpHits.clear();
  }

  //----------------------------------------------------------------
  bool PhotonBackTracker::BTRsReady() const
  {
    return !priv_OpDetBTRs.empty();
  }

  //----------------------------------------------------------------
  bool PhotonBackTracker::OpFlashToOpHitsReady() const
  {
    return !priv_OpFlashToOpHits.empty();
  }

  //----------------------------------------------------------------
  const std::vector<art::Ptr<sim::OpDetBacktrackerRecord>>& PhotonBackTracker::OpDetBTRs() const
  {
    return priv_OpDetBTRs;
  }

  //----------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTracker::TrackIdToSimSDPs_Ps(int const id) const
  {
    std::vector<const sim::SDP*> sdp_Ps;
    for (size_t odet = 0; odet < priv_OpDetBTRs.size(); ++odet) {
      const auto& pdTimeSDPmap = priv_OpDetBTRs[odet]->timePDclockSDPsMap();
      for (auto mapitr = pdTimeSDPmap.begin(); mapitr != pdTimeSDPmap.end(); mapitr++) {
        std::vector<sim::SDP> const& sdpvec = (*mapitr).second;
        for (size_t iv = 0; iv < sdpvec.size(); ++iv) {
          if (abs(sdpvec[iv].trackID) == id) sdp_Ps.push_back(&(sdpvec[iv]));
        }
      } // end loop over map from sim::OpDetBacktrackerRecord
    }   // end loop over sim::OpDetBacktrackerRecords
    return sdp_Ps;
  }

  //----------------------------------------------------------------
  art::Ptr<sim::OpDetBacktrackerRecord> PhotonBackTracker::FindOpDetBTR(int const opDetNum) const
  {
    art::Ptr<sim::OpDetBacktrackerRecord> opDet;
    for (size_t sc = 0; sc < priv_OpDetBTRs.size(); ++sc) {
      // This could become a bug. What if it occurs twice (shouldn't happen in correct
      // records, but still, no error handeling included for the situation
      if (priv_OpDetBTRs.at(sc)->OpDetNum() == opDetNum) opDet = priv_OpDetBTRs.at(sc);
    }
    if (!opDet) {
      throw cet::exception("PhotonBackTracker2") << "No sim:: OpDetBacktrackerRecord corresponding "
                                                 << "to opDetNum: " << opDetNum << "\n";
    }
    return opDet;
  }

  //----------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTracker::OpDetToTrackSDPs(int const OpDetNum,
                                                                 double const opHit_start_time,
                                                                 double const opHit_end_time) const
  {
    std::vector<sim::TrackSDP> tSDPs;
    double totalE = 0;
    try {
      const art::Ptr<sim::OpDetBacktrackerRecord> opDetBTR = FindOpDetBTR(OpDetNum);
      std::vector<sim::SDP> simSDPs =
        opDetBTR->TrackIDsAndEnergies(opHit_start_time, opHit_end_time);
      for (size_t e = 0; e < simSDPs.size(); ++e)
        totalE += simSDPs[e].energy;
      if (totalE < 1.e-5) totalE = 1.;
      for (size_t e = 0; e < simSDPs.size(); ++e) {
        if (simSDPs[e].trackID == sim::NoParticleId) continue;
        sim::TrackSDP info;
        info.trackID = std::abs(simSDPs[e].trackID);
        info.energyFrac = simSDPs[e].energy / totalE;
        info.energy = simSDPs[e].energy;
        tSDPs.push_back(info);
      }
    }
    catch (cet::exception const& e) {
      // This needs to go. Make it specific if there is a really an exception we would
      // like to catch.
      mf::LogWarning("PhotonBackTracker") << "Exception caught\n" << e;
    }
    return tSDPs;
  }

  //----------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToTrackSDPs(
    art::Ptr<recob::OpHit> const& opHit_P) const
  {
    auto OpDetNum = fWireReadoutGeom->OpDetFromOpChannel(opHit_P->OpChannel());
    const double pTime = opHit_P->PeakTime();
    const double pWidth = opHit_P->Width();
    const double start = (pTime - pWidth) * 1000 - fDelay;
    const double end = (pTime + pWidth) * 1000 - fDelay;
    return OpDetToTrackSDPs(OpDetNum, start, end);
  }

  //----------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToTrackSDPs(recob::OpHit const& opHit) const
  {
    auto OpDetNum = fWireReadoutGeom->OpDetFromOpChannel(opHit.OpChannel());
    const double pTime = opHit.PeakTime();
    const double pWidth = opHit.Width();
    const double start = (pTime - pWidth) * 1000 - fDelay;
    const double end = (pTime + pWidth) * 1000 - fDelay;
    return OpDetToTrackSDPs(OpDetNum, start, end);
  }

  //----------------------------------------------------------------
  std::vector<int> PhotonBackTracker::OpHitToTrackIds(recob::OpHit const& opHit) const
  {
    std::vector<int> retVec;
    for (auto const trackSDP : OpHitToTrackSDPs(opHit)) {
      retVec.push_back(trackSDP.trackID);
    }
    return retVec;
  }

  //----------------------------------------------------------------
  std::vector<int> PhotonBackTracker::OpHitToTrackIds(art::Ptr<recob::OpHit> const& opHit) const
  {
    return OpHitToTrackIds(*opHit);
  }

  //----------------------------------------------------------------
  std::vector<int> PhotonBackTracker::OpHitToEveTrackIds(recob::OpHit const& opHit)
  {
    std::vector<int> retVec;
    for (auto const trackSDP : OpHitToEveTrackSDPs(opHit)) {
      retVec.push_back(trackSDP.trackID);
    }
    return retVec;
  }

  //----------------------------------------------------------------
  std::vector<int> PhotonBackTracker::OpHitToEveTrackIds(art::Ptr<recob::OpHit> const& opHit_P)
  {
    return OpHitToEveTrackIds(*opHit_P);
  }

  //----------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToEveTrackSDPs(
    art::Ptr<recob::OpHit> const& opHit_P) const
  {
    return OpHitToEveTrackSDPs(*opHit_P);
  }

  //----------------------------------------------------------------
  std::vector<sim::TrackSDP> PhotonBackTracker::OpHitToEveTrackSDPs(recob::OpHit const& opHit) const
  {
    std::vector<sim::TrackSDP> trackSDPs = OpHitToTrackSDPs(opHit);

    // make a map of evd ID values and fraction of energy represented by that eve id in
    // this opHit
    std::map<int, float> eveIDtoEfrac;

    double totalE = 0.;
    for (size_t t = 0; t < trackSDPs.size(); ++t) {
      eveIDtoEfrac[(fPartInv->ParticleList()).EveId(trackSDPs[t].trackID)] += trackSDPs[t].energy;
      totalE += trackSDPs[t].energy;
    }

    // now fill the eveSDPs vector from the map
    std::vector<sim::TrackSDP> eveSDPs;
    eveSDPs.reserve(eveIDtoEfrac.size());
    for (auto itr = eveIDtoEfrac.begin(); itr != eveIDtoEfrac.end(); itr++) {
      sim::TrackSDP temp;
      temp.trackID = (*itr).first;
      temp.energyFrac = (*itr).second / totalE;
      temp.energy = (*itr).second;
      eveSDPs.push_back(std::move(temp));
    }
    return eveSDPs;
  }

  //----------------------------------------------------------------
  //TODO: Make a copy of this function that uses an allOpHits list.
  std::vector<art::Ptr<recob::OpHit>> PhotonBackTracker::TrackIdToOpHits_Ps(
    int const tkId,
    std::vector<art::Ptr<recob::OpHit>> const& hitsIn) const
  {
    // One would think we would want to have this function defined, and call this function
    // in the std::vector<tkids> to opHits, but that would require more loops (and a
    // higher overhead.) Instead, to provide this, we will just call the existing
    // std::vector<tkids>ToOpHits with an input of 1.
    std::vector<int> tkidFake(1, tkId);
    return TrackIdsToOpHits_Ps(tkidFake, hitsIn).at(0);
  }

  //----------------------------------------------------------------
  std::vector<std::vector<art::Ptr<recob::OpHit>>> PhotonBackTracker::TrackIdsToOpHits_Ps(
    std::vector<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& hitsIn) const
  {
    std::vector<std::pair<int, art::Ptr<recob::OpHit>>> opHitList;
    for (auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) {
      art::Ptr<recob::OpHit> const& opHit = *itr;
      auto OpDetNum = fWireReadoutGeom->OpDetFromOpChannel(opHit->OpChannel());
      const double pTime = opHit->PeakTime(), pWidth = opHit->Width();
      const double start = (pTime - pWidth) * 1000.0 - fDelay,
                   end = (pTime + pWidth) * 1000.0 - fDelay;
      std::vector<sim::TrackSDP> tids = OpDetToTrackSDPs(OpDetNum, start, end);
      for (auto itid = tids.begin(); itid != tids.end(); ++itid) {
        for (auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
          if (itid->trackID == *itkid and itid->energyFrac > fMinOpHitEnergyFraction)
            opHitList.push_back(std::make_pair(*itkid, opHit));
        } // itkid
      }   // itid
    }     // itr
    // now build the truOpHits vector that will be returned to the caller
    std::vector<std::vector<art::Ptr<recob::OpHit>>> truOpHits;
    // temporary vector containing opHits assigned to one MC particle
    std::vector<art::Ptr<recob::OpHit>> tmpOpHits;
    for (auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
      tmpOpHits.clear();
      for (auto itr = opHitList.begin(); itr != opHitList.end(); ++itr) {
        if (*itkid == (*itr).first) tmpOpHits.push_back((*itr).second);
      }
      truOpHits.push_back(tmpOpHits);
    }
    return truOpHits;
  }

  //----------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTracker::OpHitToSimSDPs_Ps(recob::OpHit const& opHit) const
  {
    std::vector<const sim::SDP*> retVec;
    double fPeakTime = opHit.PeakTime();
    double fWidth = opHit.Width();
    sim::OpDetBacktrackerRecord::timePDclock_t start_time =
      ((fPeakTime - fWidth) * 1000.0) - fDelay;
    sim::OpDetBacktrackerRecord::timePDclock_t end_time = ((fPeakTime + fWidth) * 1000.0) - fDelay;
    if (start_time > end_time) { throw; }

    auto const& timeSDPMap = FindOpDetBTR(fWireReadoutGeom->OpDetFromOpChannel(opHit.OpChannel()))
                               ->timePDclockSDPsMap(); //Not guranteed to be sorted.

    std::vector<const std::pair<double, std::vector<sim::SDP>>*> timePDclockSDPMap_SortedPointers;
    for (auto& pair : timeSDPMap) {
      timePDclockSDPMap_SortedPointers.push_back(&pair);
    }
    auto pairSort = [](auto& a, auto& b) { return a->first < b->first; };
    if (!std::is_sorted(timePDclockSDPMap_SortedPointers.begin(),
                        timePDclockSDPMap_SortedPointers.end(),
                        pairSort)) {
      std::sort(
        timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), pairSort);
    }
    // auto const& timeSDPMap =
    //   FindOpDetBTR(fWireReadoutGeom->OpDetFromOpChannel(opHit.OpChannel()))->timePDclockSDPsMap();

    //This section is a hack to make comparisons work right.
    std::vector<sim::SDP> dummyVec;
    std::pair<double, std::vector<sim::SDP>> start_timePair = std::make_pair(start_time, dummyVec);
    std::pair<double, std::vector<sim::SDP>> end_timePair = std::make_pair(end_time, dummyVec);
    auto start_timePair_P = &start_timePair;
    auto end_timePair_P = &end_timePair;
    //First interesting iterator.
    // auto mapFirst = timeSDPMap.lower_bound(start_time);
    auto mapFirst = std::lower_bound(timePDclockSDPMap_SortedPointers.begin(),
                                     timePDclockSDPMap_SortedPointers.end(),
                                     start_timePair_P,
                                     pairSort);
    //Last interesting iterator.
    // auto mapLast = timeSDPMap.upper_bound(end_time);
    auto mapLast =
      std::upper_bound(mapFirst, timePDclockSDPMap_SortedPointers.end(), end_timePair_P, pairSort);

    for (auto& mapitr = mapFirst; mapitr != mapLast; ++mapitr)
      // for (auto& sdp : mapitr->second)
      for (auto& sdp : (*mapitr)->second)
        retVec.push_back(&sdp);

    return retVec;
  }

  //----------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTracker::OpHitToSimSDPs_Ps(
    art::Ptr<recob::OpHit> const& opHit_P) const
  {
    return OpHitToSimSDPs_Ps(*opHit_P);
  }

  //----------------------------------------------------------------
  std::vector<double> PhotonBackTracker::SimSDPsToXYZ(std::vector<sim::SDP> const& sdps) const&
  {
    std::vector<double> xyz(3, -999.);
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double w = 0.;
    // loop over photons.
    for (auto const& sdp : sdps) {
      double weight = sdp.numPhotons;
      w += weight;
      x += weight * sdp.x;
      y += weight * sdp.y;
      z += weight * sdp.z;
    } // end loop over sim::SDPs
    // If the sum of the weights is still zero, then fail to return a value.  A hit with
    // no contributing photons does't make sense.
    if (w < 1.e-5)
      throw cet::exception("PhotonBackTracker")
        << "No sim::SDPs providing non-zero number of photons"
        << " can't determine originating location from truth\n";
    xyz[0] = x / w;
    xyz[1] = y / w;
    xyz[2] = z / w;
    return xyz;
  }

  //----------------------------------------------------------------
  std::vector<double> PhotonBackTracker::SimSDPsToXYZ(
    std::vector<const sim::SDP*> const& sdps_Ps) const&
  {
    std::vector<double> xyz(3, -999.);
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double w = 0.;
    // loop over photons.
    for (const sim::SDP* sdp_P : sdps_Ps) {
      auto& sdp = *sdp_P;
      double weight = sdp.numPhotons;
      w += weight;
      x += weight * sdp.x;
      y += weight * sdp.y;
      z += weight * sdp.z;
    } // end loop over sim::SDPs
    // If the sum of the weights is still zero, then fail to return a value.  A hit with
    // no contributing photons does't make sense.
    if (w < 1.e-5)
      throw cet::exception("PhotonBackTracker")
        << "No sim::SDPs providing non-zero number of photons"
        << " can't determine originating location from truth\n";
    xyz[0] = x / w;
    xyz[1] = y / w;
    xyz[2] = z / w;
    return xyz;
  }

  //----------------------------------------------------------------
  std::vector<double> PhotonBackTracker::OpHitToXYZ(recob::OpHit const& opHit)
  {
    return SimSDPsToXYZ(OpHitToSimSDPs_Ps(opHit));
  }

  //----------------------------------------------------------------
  std::vector<double> PhotonBackTracker::OpHitToXYZ(art::Ptr<recob::OpHit> const& opHit)
  {
    return SimSDPsToXYZ(OpHitToSimSDPs_Ps(*opHit));
  }

  //----------------------------------------------------------------
  std::vector<const sim::SDP*> PhotonBackTracker::OpHitsToSimSDPs_Ps(
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps) const
  {
    std::vector<const sim::SDP*> sdps_Ps;
    for (auto opHit_P : opHits_Ps) {
      std::vector<const sim::SDP*> to_add_sdps_Ps = OpHitToSimSDPs_Ps(opHit_P);
      sdps_Ps.insert(sdps_Ps.end(), to_add_sdps_Ps.begin(), to_add_sdps_Ps.end());
    }
    return sdps_Ps;
  }

  //----------------------------------------------------------------
  std::vector<double> PhotonBackTracker::OpHitsToXYZ(
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps) const
  {
    const std::vector<const sim::SDP*> SDPs_Ps = OpHitsToSimSDPs_Ps(opHits_Ps);
    return SimSDPsToXYZ(SDPs_Ps);
  }

  //----------------------------------------------------------------
  std::unordered_set<const sim::SDP*> PhotonBackTracker::OpHitToEveSimSDPs_Ps(
    recob::OpHit const& opHit_P)
  {
    const std::vector<int> ids = OpHitToEveTrackIds(opHit_P);
    std::unordered_set<const sim::SDP*> sdps;
    for (auto const& id : ids) {
      std::vector<const sim::SDP*> tmp_sdps = TrackIdToSimSDPs_Ps(id);
      for (const sim::SDP* tmp_sdp : tmp_sdps) {
        sdps.insert(tmp_sdp); //emplace not needed here.
      }
    }
    return sdps;
  }

  //----------------------------------------------------------------
  std::unordered_set<const sim::SDP*> PhotonBackTracker::OpHitToEveSimSDPs_Ps(
    art::Ptr<recob::OpHit>& opHit)
  {
    const std::vector<int> ids = OpHitToEveTrackIds(opHit);
    std::unordered_set<const sim::SDP*> sdps;
    for (auto const& id : ids) {
      std::vector<const sim::SDP*> tmp_sdps = TrackIdToSimSDPs_Ps(id);
      for (const sim::SDP* tmp_sdp : tmp_sdps) {
        sdps.insert(tmp_sdp); //emplace not needed here.
      }
    }
    return sdps;
  }

  //----------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfEveIds() const
  {
    return fPartInv->GetSetOfEveIds();
  }

  //----------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfTrackIds() const
  {
    return fPartInv->GetSetOfTrackIds();
  }

  //----------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfEveIds(
    std::vector<art::Ptr<recob::OpHit>> const& opHits_Ps) const
  {
    std::set<int> eveIds;
    for (auto const& opHit_P : opHits_Ps) {
      const std::vector<sim::TrackSDP> sdps = OpHitToEveTrackSDPs(opHit_P);
      for (auto const& sdp : sdps) {
        eveIds.insert(sdp.trackID);
      } //end sdps
    }   //End for hits
    return eveIds;
  }

  //----------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfEveIds(std::vector<recob::OpHit> const& opHits) const
  {
    std::set<int> eveIds;
    for (auto const& opHit : opHits) {
      const std::vector<sim::TrackSDP> sdps = OpHitToEveTrackSDPs(opHit);
      for (auto const& sdp : sdps) {
        eveIds.insert(sdp.trackID);
      } //end sdps
    }   //End for hits
    return eveIds;
  }

  //----------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfTrackIds(
    std::vector<art::Ptr<recob::OpHit>> const& opHits) const
  {
    std::set<int> tids;
    for (auto const& opHit : opHits) {
      for (auto const& sdp : OpHitToTrackSDPs(opHit)) {
        tids.insert(sdp.trackID);
      } //End for TrackSDPs
    }   //End for hits
    return tids;
  }

  //----------------------------------------------------------------
  std::set<int> PhotonBackTracker::GetSetOfTrackIds(std::vector<recob::OpHit> const& opHits) const
  {
    std::set<int> tids;
    for (auto const& opHit : opHits) {
      for (auto const& sdp : OpHitToTrackSDPs(opHit)) {
        tids.insert(sdp.trackID);
      } //End for TrackSDPs
    }   //End for hits
    return tids;
  }

  //----------------------------------------------------------------
  double PhotonBackTracker::OpHitCollectionPurity(std::set<int> const& tkIds,
                                                  std::vector<art::Ptr<recob::OpHit>> const& opHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection if the
    // EveID shows up in the input list of tkIds, then it counts
    float total = 1. * opHits.size();
    float desired = 0.;
    for (size_t h = 0; h < opHits.size(); ++h) {
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackSDPs = OpHitToTrackSDPs(opHit);
      // don't double count if this opHit has more than one of the desired track IDs
      // associated with it
      for (size_t e = 0; e < opHitTrackSDPs.size(); ++e) {
        if (tkIds.find(opHitTrackSDPs[e].trackID) != tkIds.end()) {
          desired += 1.;
          break;
        }
      }
    } // end loop over opHits
    double purity = 0.;
    if (total > 0) purity = desired / total;
    return purity;
  }

  //----------------------------------------------------------------
  double PhotonBackTracker::OpHitLightCollectionPurity(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits)
  {
    // get the list of EveIDs that correspond to the opHits in this collection if the
    // EveID shows up in the input list of tkIds, then it counts
    float total = 0;
    float desired = 0.;
    // don't have to check the view in the opHits collection because those are assumed to
    // be from the object we are testing and will the correct view by definition then.
    for (size_t h = 0; h < opHits.size(); ++h) {
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = OpHitToTrackSDPs(opHit);
      total += opHit->Area(); // sum up the charge in the cluster
      // don't double count if this opHit has more than one of the desired track IDs
      // associated with it
      for (size_t e = 0; e < opHitTrackIDs.size(); ++e) {
        if (tkIds.find(opHitTrackIDs[e].trackID) != tkIds.end()) {
          desired += opHit->Area();
          break;
        }
      }
    } // end loop over opHits
    double purity = 0.;
    if (total > 0) purity = desired / total;
    return purity;
  }

  //----------------------------------------------------------------
  double PhotonBackTracker::OpHitCollectionEfficiency(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits,
    std::vector<art::Ptr<recob::OpHit>> const& opHitsIn)
  {
    float desired = 0.;
    float total = 0.;
    for (size_t h = 0; h < opHits.size(); ++h) {
      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackSDPs = OpHitToTrackSDPs(opHit);
      // also don't double count if this opHit has more than one of the desired track IDs
      // associated with it
      for (size_t e = 0; e < opHitTrackSDPs.size(); ++e) {
        if (tkIds.find(opHitTrackSDPs[e].trackID) != tkIds.end() &&
            opHitTrackSDPs[e].energyFrac >= fMinOpHitEnergyFraction) {
          desired += 1.;
          break;
        }
      }
    } // end loop over opHits
    // now figure out how many opHits in the whole collection are associated with this id
    for (size_t h = 0; h < opHitsIn.size(); ++h) {
      art::Ptr<recob::OpHit> opHit = opHitsIn[h];
      std::vector<sim::TrackSDP> opHitTrackSDPs = OpHitToTrackSDPs(opHit);
      for (size_t e = 0; e < opHitTrackSDPs.size(); ++e) {
        // don't worry about opHits where the energy fraction for the chosen trackID is <
        // 0.1; also don't double count if this opHit has more than one of the desired
        // track IDs associated with it
        if (tkIds.find(opHitTrackSDPs[e].trackID) != tkIds.end() &&
            opHitTrackSDPs[e].energyFrac >= fMinOpHitEnergyFraction) {
          total += 1.;
          break;
        }
      }
    } // end loop over all opHits
    double efficiency = 0.;
    if (total > 0.) efficiency = desired / total;
    return efficiency;
  }

  //----------------------------------------------------------------
  double PhotonBackTracker::OpHitLightCollectionEfficiency(
    std::set<int> const& tkIds,
    std::vector<art::Ptr<recob::OpHit>> const& opHits,
    std::vector<art::Ptr<recob::OpHit>> const& opHitsIn)
  {
    float desired = 0.;
    float total = 0.;

    // don't have to check the view in the opHits collection because those are assumed to
    // be from the object we are testing and will the correct view by definition then.
    for (size_t h = 0; h < opHits.size(); ++h) {

      art::Ptr<recob::OpHit> opHit = opHits[h];
      std::vector<sim::TrackSDP> opHitTrackIDs = OpHitToTrackSDPs(opHit);

      // don't worry about opHits where the energy fraction for the chosen trackID is <
      // 0.1; also don't double count if this opHit has more than one of the desired track
      // IDs associated with it
      for (size_t e = 0; e < opHitTrackIDs.size(); ++e) {
        if (tkIds.find(opHitTrackIDs[e].trackID) != tkIds.end() &&
            opHitTrackIDs[e].energyFrac >= fMinOpHitEnergyFraction) {
          desired += opHit->Area();
          break;
        }
      }
    } // end loop over opHits
    for (size_t h = 0; h < opHitsIn.size(); ++h) {
      art::Ptr<recob::OpHit> opHit = opHitsIn[h];
      // check that we are looking at the appropriate view here in the case of 3D objects
      // we take all opHits
      std::vector<sim::TrackSDP> opHitTrackIDs = OpHitToTrackSDPs(opHit);
      for (size_t e = 0; e < opHitTrackIDs.size(); ++e) {
        // don't worry about opHits where the energy fraction for the chosen trackID is <
        // 0.1; also don't double count if this opHit has more than one of the desired
        // track IDs associated with it
        if (tkIds.find(opHitTrackIDs[e].trackID) != tkIds.end() &&
            opHitTrackIDs[e].energyFrac >= fMinOpHitEnergyFraction) {
          total += opHit->Area();
          break;
        }
      }
    } // end loop over all opHits
    double efficiency = 0.;
    if (total > 0.) efficiency = desired / total;
    return efficiency;
  }
  //--------------------------------------------------
  std::vector<art::Ptr<recob::OpHit>> PhotonBackTracker::OpFlashToOpHits_Ps(
    art::Ptr<recob::OpFlash>& flash_P) const
  {
    // There is not "non-pointer" version of this because the art::Ptr is needed to look
    // up the assn. One could loop the Ptrs and dereference them, but I will not encourage
    // the behavior by building the tool to do it.
    return priv_OpFlashToOpHits.at(flash_P);
  }

  //--------------------------------------------------
  std::vector<double> PhotonBackTracker::OpFlashToXYZ(art::Ptr<recob::OpFlash>& flash_P) const
  {
    return OpHitsToXYZ(OpFlashToOpHits_Ps(flash_P));
  }

  //--------------------------------------------------
  std::set<int> PhotonBackTracker::OpFlashToTrackIds(art::Ptr<recob::OpFlash>& flash_P) const
  {
    std::vector<art::Ptr<recob::OpHit>> opHits_Ps = OpFlashToOpHits_Ps(flash_P);
    std::set<int> ids;
    for (auto& opHit_P : opHits_Ps) {
      for (const int& id : OpHitToTrackIds(opHit_P)) {
        ids.insert(id);
      } // end for ids
    }   // end for opHits
    return ids;
  }
} // namespace
