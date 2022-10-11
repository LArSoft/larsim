////////////////////////////////////////////////////////////////////////
//Class: ISTPC
//File: ISTPC.h and ISTPC.cxx
//Description: Class including common functions needed for using the
//Hybrid Model for Photon Propagation
//UGR, 2021
////////////////////////////////////////////////////////////////////////

#ifndef IS_ISTPC_H
#define IS_ISTPC_H

#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

namespace geo {
  class GeometryCore;
}

#include <vector>

namespace larg4 {
  class ISTPC {
  public:
    explicit ISTPC(geo::GeometryCore const& geom);
    bool isScintInActiveVolume(geo::Point_t const& ScintPoint);
    static std::vector<geo::BoxBoundedGeo> extractActiveLArVolume(geo::GeometryCore const& geom);

  private:
    std::vector<geo::BoxBoundedGeo> fActiveVolumes;
  };
}
#endif
