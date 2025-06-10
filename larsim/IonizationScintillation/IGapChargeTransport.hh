#ifndef LARSIM_IGAPCHARGETRANSPORTTOOL_HH
#define LARSIM_IGAPCHARGETRANSPORTTOOL_HH

#include "larcorealg/Geometry/geo_vectors_utils.h"
#include <utility>

namespace gap {
  class IGapChargeTransport { // used to recover charge from the gap between CRPs
  public:
    virtual ~IGapChargeTransport() = default;

    /// Whether to apply gap charge transport function
    virtual bool EnableGapChargeTransport() const = 0;

    /// Return shifted position and number of moved electrons
    virtual std::pair<geo::Point_t, int> GetOffset(
      double x, double y, double z, int n) const = 0;
  };
}

#endif
