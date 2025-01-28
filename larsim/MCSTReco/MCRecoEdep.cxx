////////////////////////////////////////////////////////////////////////
//
//  MCRecoEdep source
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "MCRecoEdep.h"

#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace sim {

  namespace details {
    std::map<geo::PlaneID, size_t> createPlaneIndexMap()
    {
      auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
      std::map<geo::PlaneID, size_t> m;
      size_t i = 0;
      for (auto const& pid : wireReadoutGeom.Iterate<geo::PlaneID>()) {
        m[pid] = i;
        i++;
      }
      return m;
    }
  }

  //##################################################################
  MCRecoEdep::MCRecoEdep(fhicl::ParameterSet const& pset)
  //##################################################################
  {
    _debug_mode = pset.get<bool>("DebugMode");
    _save_mchit = pset.get<bool>("SaveMCHit");
  }

  const std::vector<sim::MCEdep>& MCRecoEdep::GetEdepArrayAt(size_t edep_index) const
  {
    if (edep_index >= _mc_edeps.size())
      throw cet::exception(__FUNCTION__) << Form("Track ID %zu not found!", edep_index);
    return _mc_edeps.at(edep_index);
  }

  std::vector<sim::MCEdep>& MCRecoEdep::__GetEdepArray__(unsigned int track_id)
  {
    if (ExistTrack(track_id)) return _mc_edeps.at((*_track_index.find(track_id)).second);
    _track_index.insert(std::pair<unsigned int, size_t>(track_id, _mc_edeps.size()));
    _mc_edeps.push_back(std::vector<sim::MCEdep>());
    return (*(_mc_edeps.rbegin()));
  }

  void MCRecoEdep::MakeMCEdep(const std::vector<sim::SimChannel>& schArray)
  {
    _mc_edeps.clear();
    _track_index.clear();

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

    // Key map to identify a unique particle energy deposition point
    std::map<std::pair<UniquePosition, unsigned int>, int> hit_index_m;

    auto pindex = details::createPlaneIndexMap();

    if (_debug_mode) std::cout << "Processing " << schArray.size() << " channels..." << std::endl;
    // Loop over channels
    for (size_t i = 0; i < schArray.size(); ++i) {

      // Get data to loop over
      auto const& sch = schArray[i];
      const auto& sch_map(sch.TDCIDEMap());
      // Channel
      UInt_t ch = sch.Channel();
      // Loop over ticks
      for (auto tdc_iter = sch_map.begin(); tdc_iter != sch_map.end(); ++tdc_iter) {
        // Loop over IDEs
        for (auto const& ide : (*tdc_iter).second) {

          int track_id = ide.trackID;
          if (track_id < 0) track_id = track_id * (-1);
          unsigned int real_track_id = track_id;

          UniquePosition pos(ide.x, ide.y, ide.z);

          int hit_index = -1;
          auto key = std::make_pair(pos, real_track_id);
          auto hit_index_track_iter = hit_index_m.find(key);
          if (hit_index_track_iter == hit_index_m.end()) {
            int new_hit_index = __GetEdepArray__(real_track_id).size();
            hit_index_m[key] = new_hit_index;
          }
          else {
            hit_index = (*hit_index_track_iter).second;
          }
          auto const pid = wireReadoutGeom.ChannelToWire(ch)[0].planeID();
          auto const channel_id = pindex[pid];
          double charge = ide.numElectrons;
          if (hit_index < 0) {
            // This particle energy deposition is never recorded so far. Create a new Edep
            __GetEdepArray__(real_track_id)
              .emplace_back(pos, pid, pindex.size(), ide.energy, charge, channel_id);
          }
          else {
            // Append charge to the relevant edep (@ hit_index)
            MCEdep& edep = __GetEdepArray__(real_track_id).at(hit_index);
            edep.deps[channel_id].charge += charge;
            edep.deps[channel_id].energy += ide.energy;
          }
        } // end looping over ides in this tick
      }   // end looping over ticks in this channel
    }     // end looping over channels

    if (_debug_mode) {
      std::cout << Form("  Collected %zu particles' energy depositions...", _mc_edeps.size())
                << std::endl;
    }
  }

  void MCRecoEdep::MakeMCEdep(const std::vector<sim::SimEnergyDeposit>& sedArray)
  {
    _mc_edeps.clear();
    _track_index.clear();

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    art::ServiceHandle<geo::Geometry const> geom;

    // Key map to identify a unique particle energy deposition point
    std::map<std::pair<UniquePosition, unsigned int>, int> hit_index_m;

    auto pindex = details::createPlaneIndexMap();

    if (_debug_mode)
      std::cout << "Processing " << sedArray.size() << " energy deposits..." << std::endl;
    // Loop over channels
    for (size_t i = 0; i < sedArray.size(); ++i) {

      // Get data to loop over
      auto const& sed = sedArray[i];

      // David Caratelli: much of the code below is taken from the module:
      // https://cdcvs.fnal.gov/redmine/projects/larsim/repository/revisions/develop/entry/larsim/ElectronDrift/SimDriftElectrons_module.cc

      // given this SimEnergyDeposit, find the TPC channel information "xyz" is the
      // position of the energy deposit in world coordinates. Note that the units of
      // distance in sim::SimEnergyDeposit are supposed to be cm.
      auto const mp = sed.MidPoint();
      // From the position in world coordinates, determine the cryostat and tpc. If
      // somehow the step is outside a tpc (e.g., cosmic rays in rock) just move on to the
      // next one.
      unsigned int cryostat = 0;
      try {
        geom->PositionToCryostatID(mp);
      }
      catch (cet::exception& e) {
        mf::LogWarning("SimDriftElectrons") << "step " // << energyDeposit << "\n"
                                            << "cannot be found in a cryostat\n"
                                            << e;
        continue;
      }
      unsigned int tpc = 0;
      try {
        geom->PositionToTPCID(mp);
      }
      catch (cet::exception& e) {
        mf::LogWarning("SimDriftElectrons") << "step " // << energyDeposit << "\n"
                                            << "cannot be found in a TPC\n"
                                            << e;
        continue;
      }
      geo::TPCID const tpcid{cryostat, tpc};

      // Define charge drift direction: driftcoordinate (x, y or z) and driftsign
      // (positive or negative). Also define coordinates perpendicular to drift direction.

      // make a collection of electrons for each plane
      for (auto const& planeid : wireReadoutGeom.Iterate<geo::PlaneID>()) {

        // grab the nearest channel to the fDriftClusterPos position
        // David Caratelli, comment begin:
        // NOTE: the below code works only when the drift coordinate is indeed in x (i.e. 0th coordinate)
        // see code linked above for a much more careful treatment of the coordinate system
        // David Caratelli, comment end.
        raw::ChannelID_t ch;
        try {
          ch = wireReadoutGeom.NearestChannel(mp, planeid);
        }
        catch (cet::exception& e) {
          mf::LogWarning("SimDriftElectrons") << "step " // << energyDeposit << "\n"
                                              << "nearest wire not in TPC\n"
                                              << e;
          continue;
        }

        int track_id = sed.TrackID();

        if (track_id < 0) track_id = track_id * (-1);
        unsigned int real_track_id = track_id;

        UniquePosition pos(mp.X(), mp.Y(), mp.Z());

        int hit_index = -1;
        auto key = std::make_pair(pos, real_track_id);
        auto hit_index_track_iter = hit_index_m.find(key);
        if (hit_index_track_iter == hit_index_m.end()) {
          int new_hit_index = __GetEdepArray__(real_track_id).size();
          hit_index_m[key] = new_hit_index;
        }
        else {
          hit_index = (*hit_index_track_iter).second;
        }
        auto const pid = wireReadoutGeom.ChannelToWire(ch)[0].planeID();
        auto const channel_id = pindex[pid];
        double charge = sed.NumElectrons();
        if (hit_index < 0) {
          // This particle energy deposition is never recorded so far. Create a new Edep
          __GetEdepArray__(real_track_id)
            .emplace_back(pos, pid, pindex.size(), sed.Energy(), charge, channel_id);
        }
        else {
          // Append charge to the relevant edep (@ hit_index)
          MCEdep& edep = __GetEdepArray__(real_track_id).at(hit_index);
          edep.deps[channel_id].charge += charge;
          edep.deps[channel_id].energy += sed.Energy();
        }
      } // end looping over planes
    }   // end looping over SimEnergyDeposits

    if (_debug_mode) {
      std::cout << Form("  Collected %zu particles' energy depositions...", _mc_edeps.size())
                << std::endl;
    }
  }

  void MCRecoEdep::MakeMCEdep(const std::vector<sim::SimEnergyDepositLite>& sedArray)
  {
    // Create a substitute array of sim::SimEnergyDeposit to avoid duplicating code...
    // Note that this makes use of the explicit conversion operator defined in
    // SimEnergyDepositLite class. Information will be partial.  Most notably for
    // MakeMCEdep, charge (numElectrons) will be 0.
    std::vector<sim::SimEnergyDeposit> new_sedArray(sedArray.begin(), sedArray.end());
    MakeMCEdep(new_sedArray);
  }

}
