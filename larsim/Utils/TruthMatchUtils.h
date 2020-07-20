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
#include "lardataobj/RecoBase/Hit.h"

// c++
#include <map>
#include <limits>
#include <vector>

namespace TruthMatchUtils {

typedef int G4ID;
typedef float EDeposit;
typedef std::map<G4ID, EDeposit> IDToEDepositMap;

const G4ID kNoG4ID = std::numeric_limits<G4ID>::lowest();           ///< The value used when no G4 ID has been found

/**
 *  @brief  The G4 ID of the true particle which deposits the most energy in the recob::Hit
 *
 *  @param  pHit the recob::Hit to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
G4ID TrueParticleID(const art::Ptr<recob::Hit> &pHit, const bool rollupUnsavedIDs);

/**
 *  @brief  The G4 ID of the true particle which deposits the most energy in a vector of recob::Hit
 *
 *  @param  pHits the recob::Hit vector to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
G4ID TrueParticleIDFromTotalTrueEnergy(const std::vector<art::Ptr<recob::Hit> >& pHits, const bool rollupUnsavedIDs);

/**
 *  @brief  The G4 ID of the true particle whose matched hits have produced the largest amount of reconstructed charge in a recob::Hit vector
 *
 *  @param  pHits the recob::Hit vector to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
G4ID TrueParticleIDFromTotalRecoCharge(const std::vector<art::Ptr<recob::Hit> >& pHits, const bool rollupUnsavedIDs);

/**
 *  @brief  The G4 ID of the true particle who has been truth-matched to the most hits in a recob::Hit vector
 *
 *  @param  pHits the recob::Hit vector to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 *
 *  @return the ID of the G4 particle
 */
G4ID TrueParticleIDFromTotalRecoHits(const std::vector<art::Ptr<recob::Hit> >& pHits, const bool rollupUnsavedIDs);

/**
 *  @brief  Fill an energy deposition map (maps G4 ID to true energy deposition) for a recob::Hit
 *
 *  @param  idToEDepMap the energy deposition map to be filled
 *  @param  pHit the recob::Hit to be assessed
 *  @param  rollupUnsavedIDs whether to squash energy deposits for non-saved G4 particles (e.g. shower secondaries) its saved ancestor particle
 */
void FillG4IDToEnergyDepositMap(IDToEDepositMap &idToEDepMap, const art::Ptr<recob::Hit> &pHit, const bool rollupUnsavedIDs);

} // namespace TruthMatchUtils

#endif // #ifndef TRUTHMATCHUTILS_H_SEEN
