/*!
 * Title:   MergeSimSources Utility Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class that merges different simulation sources together to created a combined sim list.
 * Typically just merges vectors/maps/etc together. But, if anything as a G4 trackID, applies
 * a user-defined offset to those IDs.
 *
*/

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "MergeSimSources.h"

sim::MergeSimSourcesUtility::MergeSimSourcesUtility(const std::vector<int>& offsets)
  : fG4TrackIDOffsets(offsets)
{
  Reset();
}

void sim::MergeSimSourcesUtility::Reset()
{
  fG4TrackIDRanges.clear();
  fG4TrackIDRanges.resize(
    fG4TrackIDOffsets.size(),
    std::make_pair(std::numeric_limits<int>::max(), std::numeric_limits<int>::min()));
  fMCParticleListMap.clear();
  fMCParticleListMap.resize(fG4TrackIDOffsets.size(), std::vector<size_t>());
}

void sim::MergeSimSourcesUtility::MergeMCParticles(
  std::vector<simb::MCParticle>& merged_vector,
  const std::vector<simb::MCParticle>& input_vector,
  size_t source_index)
{

  if (source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  fMCParticleListMap[source_index].resize(input_vector.size());
  merged_vector.reserve(merged_vector.size() + input_vector.size());

  std::pair<int, int> range_trackID(std::numeric_limits<int>::max(),
                                    std::numeric_limits<int>::min());

  for (size_t i_p = 0; i_p < input_vector.size(); i_p++) {
    merged_vector.emplace_back(input_vector[i_p], fG4TrackIDOffsets[source_index]);

    fMCParticleListMap[source_index][i_p] = merged_vector.size() - 1;

    if (std::abs(merged_vector.back().TrackId()) < range_trackID.first)
      range_trackID.first = std::abs(merged_vector.back().TrackId());
    if (std::abs(merged_vector.back().TrackId()) > range_trackID.second)
      range_trackID.second = std::abs(merged_vector.back().TrackId());
  }

  UpdateG4TrackIDRange(range_trackID, source_index);
}

void sim::MergeSimSourcesUtility::MergeSimChannels(std::vector<sim::SimChannel>& merged_vector,
                                                   const std::vector<sim::SimChannel>& input_vector,
                                                   size_t source_index)
{
  if (source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  merged_vector.reserve(merged_vector.size() + input_vector.size());

  std::pair<int, int> range_trackID(std::numeric_limits<int>::max(),
                                    std::numeric_limits<int>::min());

  for (auto const& simchannel : input_vector) {
    std::vector<sim::SimChannel>::iterator it =
      std::find(merged_vector.begin(), merged_vector.end(), simchannel);

    if (it == merged_vector.end()) {
      merged_vector.emplace_back(simchannel.Channel());
      it = merged_vector.end() - 1;
    }

    std::pair<int, int> thisrange =
      it->MergeSimChannel(simchannel, fG4TrackIDOffsets[source_index]);
    if (std::abs(thisrange.first) < std::abs(range_trackID.first))
      range_trackID.first = std::abs(thisrange.first);
    if (std::abs(thisrange.second) > std::abs(range_trackID.second))
      range_trackID.second = std::abs(thisrange.second);
  }

  UpdateG4TrackIDRange(range_trackID, source_index);
}

void sim::MergeSimSourcesUtility::MergeAuxDetSimChannels(
  std::vector<sim::AuxDetSimChannel>& merged_vector,
  const std::vector<sim::AuxDetSimChannel>& input_vector,
  size_t source_index)
{
  if (source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  merged_vector.reserve(merged_vector.size() + input_vector.size());

  std::pair<int, int> range_trackID(std::numeric_limits<int>::max(),
                                    std::numeric_limits<int>::min());

  for (auto const& simchannel : input_vector) {
    std::vector<sim::AuxDetSimChannel>::iterator it =
      std::find(merged_vector.begin(), merged_vector.end(), simchannel);

    if (it == merged_vector.end()) {
      merged_vector.emplace_back(simchannel.AuxDetID(), simchannel.AuxDetSensitiveID());
      it = merged_vector.end() - 1;
    }

    // re-make the AuxDetSimChannel with both pairs of AuxDetIDEs
    int offset = fG4TrackIDOffsets[source_index];
    std::vector<sim::AuxDetIDE> all_ides = it->AuxDetIDEs();
    for (const sim::AuxDetIDE& ide : simchannel.AuxDetIDEs()) {
      all_ides.emplace_back(ide, offset);

      auto tid = std::abs(ide.trackID) + offset;

      if (tid < range_trackID.first) range_trackID.first = tid;
      if (tid > range_trackID.second) range_trackID.second = tid;
    }

    *it = sim::AuxDetSimChannel(
      simchannel.AuxDetID(), std::move(all_ides), simchannel.AuxDetSensitiveID());
  }

  UpdateG4TrackIDRange(range_trackID, source_index);
}

void sim::MergeSimSourcesUtility::MergeSimPhotons(std::vector<sim::SimPhotons>& merged_vector,
                                                  const std::vector<sim::SimPhotons>& input_vector)
{

  merged_vector.reserve(merged_vector.size() + input_vector.size());

  for (auto const& simphotons : input_vector) {
    std::vector<sim::SimPhotons>::iterator it =
      std::find(merged_vector.begin(), merged_vector.end(), simphotons);

    if (it == merged_vector.end()) {
      merged_vector.emplace_back(simphotons.OpChannel());
      it = merged_vector.end() - 1;
    }

    *it += simphotons;
  }
}

void sim::MergeSimSourcesUtility::MergeSimPhotonsLite(
  std::vector<sim::SimPhotonsLite>& merged_vector,
  const std::vector<sim::SimPhotonsLite>& input_vector)
{

  merged_vector.reserve(merged_vector.size() + input_vector.size());

  for (auto const& simphotons : input_vector) {
    std::vector<sim::SimPhotonsLite>::iterator it =
      std::find(merged_vector.begin(), merged_vector.end(), simphotons);

    if (it == merged_vector.end()) {
      merged_vector.emplace_back(simphotons.OpChannel);
      it = merged_vector.end() - 1;
    }

    *it += simphotons;
  }
}

void sim::MergeSimSourcesUtility::MergeSimEnergyDeposits(
  std::vector<sim::SimEnergyDeposit>& dest,
  const std::vector<sim::SimEnergyDeposit>& src,
  std::size_t source_index) const
{

  int const offset = fG4TrackIDOffsets.at(source_index);
  auto const offsetEDepID = [offset](sim::SimEnergyDeposit const& edep) {
    return sim::MergeSimSourcesUtility::offsetSimEnergyDepositTrackID(edep, offset);
  };

  dest.reserve(dest.size() + src.size());
  std::transform(begin(src), end(src), back_inserter(dest), offsetEDepID);
}

void sim::MergeSimSourcesUtility::MergeAuxDetHits(std::vector<sim::AuxDetHit>& dest,
                                                  const std::vector<sim::AuxDetHit>& src,
                                                  std::size_t source_index) const
{

  int const offset = fG4TrackIDOffsets.at(source_index);
  auto const offsetAuxDetHitID = [offset](sim::AuxDetHit const& adh) {
    return sim::MergeSimSourcesUtility::offsetAuxDetHitTrackID(adh, offset);
  };

  dest.reserve(dest.size() + src.size());
  std::transform(begin(src), end(src), back_inserter(dest), offsetAuxDetHitID);
}

void sim::MergeSimSourcesUtility::UpdateG4TrackIDRange(std::pair<int, int> newrange,
                                                       size_t source_index)
{
  if (source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesUtility: Source index out of range!");

  if (newrange.first >= fG4TrackIDRanges[source_index].first &&
      newrange.second <= fG4TrackIDRanges[source_index].second)
    return;

  for (size_t i = 0; i < fG4TrackIDRanges.size(); i++) {
    if (i == source_index) continue;

    if ((newrange.first >= fG4TrackIDRanges[i].first &&
         newrange.first <= fG4TrackIDRanges[i].second) ||
        (newrange.second >= fG4TrackIDRanges[i].first &&
         newrange.second <= fG4TrackIDRanges[i].second)) {
      std::stringstream ss;
      ss << "ERROR in MergeSimSourcesUtility: Source trackIDs overlap!"
         << "\n\t" << i << "\t" << fG4TrackIDRanges[i].first << " " << fG4TrackIDRanges[i].second
         << "\n\t"
         << "n\t" << newrange.first << " " << newrange.second;
      throw std::runtime_error(ss.str());
    }
  }

  if (newrange.first < fG4TrackIDRanges[source_index].first)
    fG4TrackIDRanges[source_index].first = newrange.first;
  if (newrange.second > fG4TrackIDRanges[source_index].second)
    fG4TrackIDRanges[source_index].second = newrange.second;
}

sim::SimEnergyDeposit sim::MergeSimSourcesUtility::offsetSimEnergyDepositTrackID(
  sim::SimEnergyDeposit const& edep,
  int offset)
{

  auto tid = (edep.TrackID() >= 0) ? (edep.TrackID() + offset) : (edep.TrackID() - offset);

  return sim::SimEnergyDeposit{
    edep.NumPhotons(),      // np
    edep.NumElectrons(),    // ne
    edep.ScintYieldRatio(), // sy
    edep.Energy(),          // e
    edep.Start(),           // start
    edep.End(),             // end
    edep.T0(),              // t0
    edep.T1(),              // t1
    tid,                    // id
    edep.PdgCode()          // pdg
  };
} // sim::MergeSimSourcesUtility::offsetTrackID()

sim::AuxDetHit sim::MergeSimSourcesUtility::offsetAuxDetHitTrackID(sim::AuxDetHit const& adh,
                                                                   int offset)
{

  auto tid = (adh.GetTrackID() >= 0) ? (adh.GetTrackID() + offset) : (adh.GetTrackID() - offset);

  return sim::AuxDetHit{
    adh.GetID(), // copy number
    tid,         // g4 track id
    adh.GetEnergyDeposited(),
    adh.GetEntryX(),
    adh.GetEntryY(),
    adh.GetEntryZ(),
    adh.GetEntryT(),
    adh.GetExitX(),
    adh.GetExitY(),
    adh.GetExitZ(),
    adh.GetExitT(),
    adh.GetExitMomentumX(),
    adh.GetExitMomentumY(),
    adh.GetExitMomentumZ(),
  };
} // sim::MergeSimSourcesUtility::offsetAuxDetHitTrackID()
