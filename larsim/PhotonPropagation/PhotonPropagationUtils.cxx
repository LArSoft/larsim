////////////////////////////////////////////////////////////////////////
// Description:
// Utility functions
////////////////////////////////////////////////////////////////////////
#include "larsim/PhotonPropagation/PhotonPropagationUtils.h"
#include <cmath>

namespace phot {

  //......................................................................
  // fast approximation to acos, This implementation foillows the functional
  // form from C. Hastings, Jr., "Approximations for Digital Computers",
  // Princeton University Press, 1955.
  // This differs from that version by using a new fit to the same functional
  // form, but done in double precision. The resulting parameters differ, and
  // result in the same calculation speed but yield a better approximation.
  //
  // The functional form for the fit is:
  //    acos(x) = sqrt(1-x) * a0 * (x + a1 * (x + a2 * (x + a3)))
  // The fit done is *not* a least squares fit.
  // It is a fit that minimizes the value of the maximum absolute deviation
  // between the fitted function and the C math library's implementation of
  // acos(double).
  //
  double fast_acos(double xin)
  {
    double const a3 = -2.08730442907856008e-02;
    double const a2 =  7.68769404161671888e-02;
    double const a1 = -2.12871094165645952e-01;
    double const a0 =  1.57075835365209659e+00;
    double const x = std::abs(xin);
    double ret = a3;
    ret *= x;
    ret += a2;
    ret *= x;
    ret += a1;
    ret *= x;
    ret += a0;
    ret *= std::sqrt(1.0-x);
    if (xin >= 0) return ret;
    return M_PI - ret;
  }

  //......................................................................
  // Returns interpolated value at x from parallel arrays ( xData, yData )
  // Assumes that xData has at least two elements, is sorted and is strictly
  // monotonic increasing boolean argument extrapolate determines behaviour
  // beyond ends of array (if needed)
  double interpolate(const std::vector<double>& xData,
                     const std::vector<double>& yData,
                     const double x,
                     const bool extrapolate,
                     size_t i)
  {
    if (i == 0) {
      size_t size = xData.size();
      if (x >= xData[size - 2]) { // special case: beyond right end
        i = size - 2;
      }
      else {
        while (x > xData[i + 1])
          i++;
      }
    }
    double xL = xData[i];
    double xR = xData[i + 1];
    double yL = yData[i];
    double yR = yData[i + 1]; // points on either side (unless beyond ends)
    if (!extrapolate) {       // if beyond ends of array and not extrapolating
      if (x < xL) return yL;
      if (x > xR) return yL;
    }
    const double dydx = (yR - yL) / (xR - xL); // gradient
    return yL + dydx * (x - xL);               // linear interpolation
  }

  double interpolate2(const std::vector<double>& xDistances,
                      const std::vector<double>& rDistances,
                      const std::vector<std::vector<std::vector<double>>>& parameters,
                      const double x,
                      const double r,
                      const size_t k)
  {
    // interpolate in x for each r bin, for angle bin k
    const size_t nbins_r = parameters[k].size();
    std::vector<double> interp_vals(nbins_r, 0.);

    size_t idx = 0;
    size_t size = xDistances.size();
    if (x >= xDistances[size - 2])
      idx = size - 2;
    else {
      while (x > xDistances[idx + 1])
        idx++;
    }
    for (size_t i = 0; i < nbins_r; ++i) {
      interp_vals[i] = interpolate(xDistances, parameters[k][i], x, false, idx);
    }

    // interpolate in r
    double border_correction = interpolate(rDistances, interp_vals, r, false);
    return border_correction;
  }

  //......................................................................
  void interpolate3(std::array<double, 3>& inter,
                    const std::vector<double>& xData,
                    const std::vector<double>& yData1,
                    const std::vector<double>& yData2,
                    const std::vector<double>& yData3,
                    const double x,
                    const bool extrapolate)
  {
    size_t size = xData.size();
    size_t i = 0;               // find left end of interval for interpolation
    if (x >= xData[size - 2]) { // special case: beyond right end
      i = size - 2;
    }
    else {
      while (x > xData[i + 1])
        i++;
    }
    double xL = xData[i];
    double xR = xData[i + 1]; // points on either side (unless beyond ends)
    double yL1 = yData1[i];
    double yR1 = yData1[i + 1];
    double yL2 = yData2[i];
    double yR2 = yData2[i + 1];
    double yL3 = yData3[i];
    double yR3 = yData3[i + 1];

    if (!extrapolate) { // if beyond ends of array and not extrapolating
      if (x < xL) {
        inter[0] = yL1;
        inter[1] = yL2;
        inter[2] = yL3;
        return;
      }
      if (x > xR) {
        inter[0] = yL1;
        inter[1] = yL2;
        inter[2] = yL3;
        return;
      }
    }
    const double m = (x - xL) / (xR - xL);
    inter[0] = m * (yR1 - yL1) + yL1;
    inter[1] = m * (yR2 - yL2) + yL2;
    inter[2] = m * (yR3 - yL3) + yL3;
  }

} // namespace phot
