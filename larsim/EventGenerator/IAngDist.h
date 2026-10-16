/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : IAngDist.h
 * @created     : Friday Oct 16, 2026 15:30:33 CEST
 * @description : Interface for algorithms that generate a particle direction.
 */

#ifndef LARSIM_IANGDIST_H

#define LARSIM_IANGDIST_H

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Vector_t
#include "CLHEP/Random/RandomEngine.h"

namespace evgen {
  /// Interface for algorithms that generate a particle direction.
  class IAngDist {
  public:
    virtual ~IAngDist() noexcept = default;

    /// Returns a unit direction vector in the detector frame (dimensionless).
    virtual geo::Vector_t generate(CLHEP::HepRandomEngine& engine) = 0;
  };
} 

#endif /* end of include guard LARSIM_IANGDIST_H */

