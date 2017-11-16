/**
 * @file   larsim/MCDumpers/MCDumperUtils.h
 * @brief  Utility functions to print MC truth information.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 2, 2017
 * 
 * Utilities for dumping on screen of Monte Carlo data.
 * 
 */

#ifndef LARSIM_MCDUMPERS_MCDUMPERUTILS_H
#define LARSIM_MCDUMPERS_MCDUMPERUTILS_H

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h" // simb::Origin_t

// C/C++ standard libraries
#include <string>


namespace sim {
  
  /// @{
  /// @name Function to name `simb::MCParticle` enumerators and codes.
  
  /// Returns a string with the name of particle the specified with PDG ID.
  std::string ParticleName(int pigid);
  
  /// Describes the status of a particle (`simb::MCParticle::StatusCode()`).
  std::string ParticleStatusName(int code);
  
  /// @}
  
  
  /// @{
  /// @name Functions to name `simb::MCTruth` enumerators and codes.
  
  /// Returns a string representing the specified process origin.
  std::string TruthOriginName(simb::Origin_t origin);
  
  /// Returns a string representing the specified process from `simb::MCTruth`
  /// (CC or NC, nothing fancy).
  std::string TruthCCNCname(int ccnc);
  
  /// Returns the "mode" of the reaction (a lesser version of interaction type).
  std::string TruthReactionMode(int mode);
  
  /// Returns a string representing the specified interaction type as in 
  /// `simb::MCTruth` convention.
  std::string TruthInteractionTypeName(int type);
  
  /// @}
  
  
} // namespace sim


#endif // LARSIM_MCDUMPERS_MCDUMPERUTILS_H
