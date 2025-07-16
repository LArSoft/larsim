// default optical path tool
// all photon detectors are visible from all scintillation emission points

#ifndef StandardOpticalPath_H
#define StandardOpticalPath_H

#include "art/Utilities/ToolMacros.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"

#include <iostream>

namespace phot {
  class StandardOpticalPath : public phot::OpticalPath {
  public:
    explicit StandardOpticalPath(fhicl::ParameterSet const& ps){};
    ~StandardOpticalPath() noexcept override = default;

    const bool isOpDetVisible(geo::Point_t const& ScintPoint,
                              geo::Point_t const& OpDetPoint) override
    {
      return true;
    }
  };
}

#endif