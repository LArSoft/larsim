////////////////////////////////////////////////////////////////////////
// Description:
// Utility functions
////////////////////////////////////////////////////////////////////////
#ifndef PhotonPropagationUtils_H
#define PhotonPropagationUtils_H

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace phot {
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
