#ifndef IS_ISTPC_H
#define IS_ISTPC_H


#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {
  class ISTPC {
  public:

    ISTPC();
    std::vector<geo::BoxBoundedGeo> fActiveVolumes;

     bool isScintInActiveVolume(geo::Point_t const& ScintPoint);
 
    static std::vector<geo::BoxBoundedGeo> extractActiveVolumes(geo::GeometryCore const& geom);
 
      
  };
}
#endif

