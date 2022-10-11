/**
 * @file larsim/Utils/TruthMatchUtils.h
 *
 * @brief Utilities for matching a recob::Hit or vector of recob::Hit to the ID of the most significantly contributing G4 particle.  Various functions using different definitions of 'most significantly contributing' are provided.
 *
 * @author Dom Brailsford (d.brailsford@lancaster.ac.uk)
 *
 * $log: $
 */
#ifndef TRUTHMATCHUTILS_H_SEEN
#define TRUTHMATCHUTILS_H_SEEN

// LArSoft
namespace recob {
  class Hit;
}

// art
#include "canvas/Persistency/Common/Ptr.h"

// c++
#include <map>
#include <vector>

namespace detinfo {
  class DetectorClocksData;
}

namespace TruthMatchUtils {

  typedef int G4ID;
  typedef float EDeposit;
  typedef std::map<G4ID, EDeposit> IDToEDepositMap;

  /**
 *  @brief  Test whether a G4ID returned by the TruthMatchUtils functions is valid
 *
 *  @param  g4ID the G4ID to be tested
 *
 *  @return validity bool
 */
  bool Valid(const G4ID g4ID) noexcept;

  /**
 *  @brief  The G4 ID of the true particle which deposits the most energy in the recob::Hit
 *
 *  @param  pHit the recob::Hit to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
  G4ID TrueParticleID(detinfo::DetectorClocksData const& clockData,
                      const art::Ptr<recob::Hit>& pHit,
                      const bool rollupUnsavedIDs);

  /**
 *  @brief  The G4 ID of the true particle which deposits the most energy in a vector of recob::Hit
 *
 *  @param  pHits the recob::Hit vector to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
  G4ID TrueParticleIDFromTotalTrueEnergy(detinfo::DetectorClocksData const& clockData,
                                         const std::vector<art::Ptr<recob::Hit>>& pHits,
                                         const bool rollupUnsavedIDs);

  /**
 *  @brief  The G4 ID of the true particle whose matched hits have produced the largest amount of reconstructed charge in a recob::Hit vector
 *
 *  @param  pHits the recob::Hit vector to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
  G4ID TrueParticleIDFromTotalRecoCharge(detinfo::DetectorClocksData const& clockData,
                                         const std::vector<art::Ptr<recob::Hit>>& pHits,
                                         const bool rollupUnsavedIDs);

  /**
 *  @brief  The G4 ID of the true particle who has been truth-matched to the most hits in a recob::Hit vector
 *
 *  @param  pHits the recob::Hit vector to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
  G4ID TrueParticleIDFromTotalRecoHits(detinfo::DetectorClocksData const& clockData,
                                       const std::vector<art::Ptr<recob::Hit>>& pHits,
                                       const bool rollupUnsavedIDs);

  /**
 *  @brief  Fill an energy deposition map (maps G4 ID to true energy deposition) for a recob::Hit
 *
 *  @param  idToEDepMap the energy deposition map to be filled
 *  @param  pHit the recob::Hit to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 */
  void FillG4IDToEnergyDepositMap(IDToEDepositMap& idToEDepMap,
                                  detinfo::DetectorClocksData const& clockData,
                                  const art::Ptr<recob::Hit>& pHit,
                                  const bool rollupUnsavedIDs);

} // namespace TruthMatchUtils

#endif // #ifndef TRUTHMATCHUTILS_H_SEEN
