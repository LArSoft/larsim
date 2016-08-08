

#ifndef BEAMTYPES_H
#define BEAMTYPES_H

namespace sim {

  /// Defines category of beams to be stored in sim::BeamGateInfo
  enum BeamType_t {
    kUnknown=0,  ///< Unknown beam type
    kBNB,        ///< BNB
    kNuMI,       ///< NuMI 
    kBeamTypeMax ///< Max value of enum for iteration
  };
}

#endif
