////////////////////////////////////////////////////////////////////////
// Description:
// Utility functions
////////////////////////////////////////////////////////////////////////
#ifndef PhotonPropagationUtils_H
#define PhotonPropagationUtils_H

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace geo {
  class OpDetGeo;
}

namespace phot {
  // Light-collecting window of a box-shaped optical detector, in world coordinates.
  struct OpDetWindow {
    geo::Vector_t normal;   // unit normal (sign not meaningful)
    geo::Vector_t u_h, u_w; // unit vectors along the two in-plane edges
    double h, w;            // edge lengths along u_h and u_w
    int orientation;        // world axis the normal lies along: 0 = X, 1 = Y, 2 = Z
  };

  // OpDetGeo sizes are those of the box in its local frame and do not change when
  // the volume is rotated. The normal is taken to be the thinnest box axis, which
  // holds under any rotation, and is then transformed to world coordinates.
  // Only valid for box-shaped optical detectors (geo::OpDetGeo::isBar()).
  OpDetWindow opDetWindow(geo::OpDetGeo const& opDet);

  double fast_acos(double x);
  double interpolate(const std::vector<double>& xData,
                     const std::vector<double>& yData,
                     const double x,
                     const bool extrapolate,
                     size_t i = 0);
  double interpolate2(const std::vector<double>& xDistances,
                      const std::vector<double>& rDistances,
                      const std::vector<std::vector<std::vector<double>>>& parameters,
                      const double x,
                      const double r,
                      const size_t k);
  void interpolate3(std::array<double, 3>& inter,
                    const std::vector<double>& xData,
                    const std::vector<double>& yData1,
                    const std::vector<double>& yData2,
                    const std::vector<double>& yData3,
                    const double x,
                    const bool extrapolate);

  // implements relative method - do not use for comparing with zero
  // use this most of the time, tolerance needs to be meaningful in your context
  template <typename TReal>
  inline constexpr static bool
  isApproximatelyEqual(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    TReal diff = std::fabs(a - b);
    if (diff <= tolerance) return true;
    if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
    return false;
  }

  // supply tolerance that is meaningful in your context
  // for example, default tolerance may not work if you are comparing double with
  // float
  template <typename TReal>
  inline constexpr static bool isApproximatelyZero(
    TReal a,
    TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    if (std::fabs(a) <= tolerance) return true;
    return false;
  }

  // use this when you want to be on safe side
  // for example, don't start rover unless signal is above 1
  template <typename TReal>
  inline constexpr static bool
  isDefinitelyLessThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    TReal diff = a - b;
    if (diff < tolerance) return true;
    if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
    return false;
  }

  template <typename TReal>
  inline constexpr static bool
  isDefinitelyGreaterThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
  {
    TReal diff = a - b;
    if (diff > tolerance) return true;
    if (diff > std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
    return false;
  }

} // namespace phot
#endif
