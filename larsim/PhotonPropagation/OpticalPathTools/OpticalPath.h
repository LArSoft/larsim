// Defines whether a photon detector is visible from a given scintillation emission point
// used for the semi-analytical model fast optical simulation (PDFastSimPAR)

#ifndef OpticalPath_H
#define OpticalPath_H

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

namespace phot {
  class OpticalPath {
  public:
    virtual ~OpticalPath() noexcept = default;
    virtual const bool isOpDetVisible(geo::Point_t const& ScintPoint,
                                      geo::Point_t const& OpDetPoint) = 0;
  };
}

#endif