#ifndef MERGESIMSOURCES_H
#define MERGESIMSOURCES_H

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

#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/ParticleAncestryMap.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include <utility> // std::pair<>
#include <vector>

namespace sim {

  class MergeSimSourcesUtility {

  public:
    MergeSimSourcesUtility(std::vector<int> const&);

    void Reset();

    void MergeMCParticles(std::vector<simb::MCParticle>&,
                          const std::vector<simb::MCParticle>&,
                          size_t);

    void MergeSimChannels(std::vector<sim::SimChannel>&,
                          const std::vector<sim::SimChannel>&,
                          size_t,
			  bool);

    void MergeAuxDetSimChannels(std::vector<sim::AuxDetSimChannel>&,
                                const std::vector<sim::AuxDetSimChannel>&,
                                size_t);

    void MergeSimPhotons(std::vector<sim::SimPhotons>&, const std::vector<sim::SimPhotons>&);

    void MergeSimPhotonsLite(std::vector<sim::SimPhotonsLite>&,
                             const std::vector<sim::SimPhotonsLite>&);

    void MergeSimEnergyDeposits(std::vector<sim::SimEnergyDeposit>&,
                                const std::vector<sim::SimEnergyDeposit>&,
                                size_t) const;

    void MergeAuxDetHits(std::vector<sim::AuxDetHit>&,
                         const std::vector<sim::AuxDetHit>&,
                         size_t) const;

    void MergeParticleAncestryMaps(std::vector<sim::ParticleAncestryMap>&,
                                   const sim::ParticleAncestryMap&,
                                   size_t) const;

    const std::vector<std::vector<size_t>>& GetMCParticleListMap() { return fMCParticleListMap; }

  private:
    std::vector<int> fG4TrackIDOffsets;
    std::vector<std::pair<int, int>> fG4TrackIDRanges;

    std::vector<std::vector<size_t>> fMCParticleListMap;

    void UpdateG4TrackIDRange(std::pair<int, int>, size_t);

    static sim::SimEnergyDeposit offsetSimEnergyDepositTrackID(sim::SimEnergyDeposit const&, int);

    static sim::AuxDetHit offsetAuxDetHitTrackID(sim::AuxDetHit const&, int);

    static sim::ParticleAncestryMap offsetParticleAncestryMapTrackID(
      sim::ParticleAncestryMap const&,
      int);

  }; //end MergeSimSourcesUtility class

} //end namespace sim

#endif
