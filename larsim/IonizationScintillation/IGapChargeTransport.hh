#ifndef LARSIM_IGAPCHARGETRANSPORTTOOL_HH
#define LARSIM_IGAPCHARGETRANSPORTTOOL_HH

#include "larcorealg/Geometry/geo_vectors_utils.h"
#include <utility>

namespace gap {
  class IGapChargeTransport { // used to recover charge from the gap between CRPs
  public:
    virtual ~IGapChargeTransport() = default;

    /// Return shifted position and number of moved electrons
    virtual std::pair<geo::Point_t, float> GetOffset(double x,
                                                     double y,
                                                     double z,
                                                     float n) const = 0;

    //Apply the tool in this volume
    virtual std::string Volume() const = 0;
    virtual double MaxGap() const = 0;
  };
}

#endif
