#include "PropagationTimeModel.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h" // geo::vect::toTVector3()

// support libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

#include "CLHEP/Random/RandFlat.h"

#include <iostream>

// constructor
PropagationTimeModel::PropagationTimeModel(fhicl::ParameterSet VUVTimingParams, fhicl::ParameterSet VISTimingParams, CLHEP::HepRandomEngine& ScintTimeEngine, bool doReflectedLight, bool GeoPropTimeOnly)
  :
  fVUVTimingParams(VUVTimingParams)
  , fVISTimingParams(VISTimingParams)
  , fdoReflectedLight(doReflectedLight)
  , fGeoPropTimeOnly(GeoPropTimeOnly)
  , fISTPC{*(lar::providerFrom<geo::Geometry>())}
  , fScintTimeEngine(ScintTimeEngine)
  , fUniformGen(fScintTimeEngine)
{
  // initialise parameters and geometry
  mf::LogInfo("PropagationTimeModel") << "Photon propagation time model initalized." << std::endl;
  Initialization();
}

// initialization
void
PropagationTimeModel::Initialization()
{
  // load photon propagation time model parameters
  if (!fGeoPropTimeOnly) {

    // Direct / VUV
    mf::LogInfo("PropagationTimeModel") << "Using VUV timing parameterization";

    fparameters[0] = std::vector(1, fVUVTimingParams.get<std::vector<double>>("Distances_landau"));
    fparameters[1] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Norm_over_entries");
    fparameters[2] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Mpv");
    fparameters[3] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Width");
    fparameters[4] = std::vector(1, fVUVTimingParams.get<std::vector<double>>("Distances_exp"));
    fparameters[5] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Slope");
    fparameters[6] = fVUVTimingParams.get<std::vector<std::vector<double>>>("Expo_over_Landau_norm");

    fstep_size                = fVUVTimingParams.get<double>("step_size");
    fmax_d                    = fVUVTimingParams.get<double>("max_d");
    fmin_d                    = fVUVTimingParams.get<double>("min_d");
    fvuv_vgroup_mean          = fVUVTimingParams.get<double>("vuv_vgroup_mean");
    fvuv_vgroup_max           = fVUVTimingParams.get<double>("vuv_vgroup_max");
    finflexion_point_distance = fVUVTimingParams.get<double>("inflexion_point_distance");
    fangle_bin_timing_vuv     = fVUVTimingParams.get<double>("angle_bin_timing_vuv");

    // create vector of empty TF1s that will be replaced with the sampled
    // parameterisations that are generated as they are required
    const size_t num_params = (fmax_d - fmin_d) / fstep_size; // for d < fmin_d, no parameterisaton, a delta function is used instead
    const size_t num_angles = std::round(90/fangle_bin_timing_vuv);
    fVUV_timing = std::vector(num_angles, std::vector(num_params, TF1()));

    // initialise vectors to contain range parameterisations sampled to in each case
    // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling
    // is regenerated, this is the slow part!
    fVUV_max = std::vector(num_angles, std::vector(num_params, 0.0));
    fVUV_min = std::vector(num_angles, std::vector(num_params, 0.0));

    // generate VUV parameters
    for (size_t angle_bin=0; angle_bin < num_angles; ++angle_bin) {
      for (size_t index=0; index < num_params; ++index) {
        generateParam(index, angle_bin);
      }
    }

    // Reflected / Visible
    if (fdoReflectedLight) {
      mf::LogInfo("PropagationTimeModel") << "Using VIS (reflected) timing parameterization";

      fdistances_refl        = fVISTimingParams.get<std::vector<double>>("Distances_refl");
      fradial_distances_refl = fVISTimingParams.get<std::vector<double>>("Distances_radial_refl");
      fcut_off_pars          = fVISTimingParams.get<std::vector<std::vector<std::vector<double>>>>("Cut_off");
      ftau_pars              = fVISTimingParams.get<std::vector<std::vector<std::vector<double>>>>("Tau");
      fvis_vmean             = fVISTimingParams.get<double>("vis_vmean");
      fangle_bin_timing_vis  = fVISTimingParams.get<double>("angle_bin_timing_vis");
    }

  }

  if (fGeoPropTimeOnly) {
    mf::LogInfo("PropagationTimeModel") << "Using geometric VUV time propagation";
      fvuv_vgroup_mean          = fVUVTimingParams.get<double>("vuv_vgroup_mean");
  }

  // access information from the geometry service
  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

  // get TPC information
  fplane_depth = std::abs(geom.TPC(0, 0).GetCathodeCenter().X());
  fActiveVolumes = fISTPC.extractActiveLArVolume(geom);
  fcathode_centre = {geom.TPC(0, 0).GetCathodeCenter().X(),
                     fActiveVolumes[0].CenterY(),
                     fActiveVolumes[0].CenterZ()};

  // get PDS information
  nOpDets = geom.NOpDets();

  fOpDetOrientation.reserve(nOpDets); fOpDetCenter.reserve(nOpDets);
  for (size_t const i : util::counter(nOpDets)) {

    geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpDet(i);
    fOpDetCenter.push_back(opDet.GetCenter());

    if (opDet.isSphere()) {  // dome PMTs
      fOpDetOrientation.push_back(0);
    }
    else if (opDet.isBar()) {
      // determine orientation to get correction OpDet dimensions
      if (opDet.Width() > opDet.Height()) { // laterals, Y dimension smallest
        fOpDetOrientation.push_back(1);
      }
      else {  // anode/cathode (default), X dimension smallest
        fOpDetOrientation.push_back(0);
      }
    }
    else {
      fOpDetOrientation.push_back(0);
    }
  }
}

//......................................................................
// Propagation time calculation function
void
PropagationTimeModel:: propagationTime(std::vector<double>& arrival_time_dist,
                                       geo::Point_t const& x0,
                                       const size_t OpChannel,
                                       bool Reflected)
{
  if (!fGeoPropTimeOnly) {
    // Get VUV photons transport time distribution from the parametrization
    geo::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
    if (!Reflected) {
      double distance = std::hypot(x0.X() - opDetCenter.X(), x0.Y() - opDetCenter.Y(), x0.Z() - opDetCenter.Z());
      double cosine;
      if (fOpDetOrientation[OpChannel] == 1) cosine = std::abs(x0.Y() - opDetCenter.Y()) / distance;
      else cosine = std::abs(x0.X() - opDetCenter.X()) / distance;

      double theta = fast_acos(cosine)*180./CLHEP::pi;
      int angle_bin = theta/fangle_bin_timing_vuv;
      getVUVTimes(arrival_time_dist, distance, angle_bin); // in ns
    }
    else {
      getVISTimes(arrival_time_dist, geo::vect::toTVector3(x0),
                  geo::vect::toTVector3(opDetCenter)); // in ns
    }
  }
  else if (fGeoPropTimeOnly && !Reflected) {
    // Get VUV photons arrival time geometrically
    geo::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
    double distance = std::hypot(x0.X() - opDetCenter.X(), x0.Y() - opDetCenter.Y(), x0.Z() - opDetCenter.Z());
    getVUVTimesGeo(arrival_time_dist, distance); // in ns
  }
  else {
    throw cet::exception("PropagationTimeModel")
      << "Propagation time model not found.";
  }
}

//......................................................................
// VUV propagation times calculation function
void
PropagationTimeModel::getVUVTimes(std::vector<double>& arrivalTimes, const double distance, const size_t angle_bin)
{
  if (distance < fmin_d) {
    // times are fixed shift i.e. direct path only
    double t_prop_correction = distance / fvuv_vgroup_mean;
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] = t_prop_correction;
    }
  }
  else {
    // determine nearest parameterisation in discretisation
    int index = std::round((distance - fmin_d) / fstep_size);
    // randomly sample parameterisation for each photon
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] = fVUV_timing[angle_bin][index].GetRandom(fVUV_min[angle_bin][index], fVUV_max[angle_bin][index]);
    }
  }
}

//......................................................................
// VUV arrival times calculation function - pure geometric approximation for use in Xenon doped scenarios
void
PropagationTimeModel::getVUVTimesGeo(std::vector<double>& arrivalTimes, const double distance)
{
  // times are fixed shift i.e. direct path only
  double t_prop_correction = distance / fvuv_vgroup_mean;
  for (size_t i = 0; i < arrivalTimes.size(); ++i) {
    arrivalTimes[i] = t_prop_correction;
  }
}

//......................................................................
// VUV propagation times parameterization generation function
void
PropagationTimeModel::generateParam(const size_t index, const size_t angle_bin)
{
  // get distance
  double distance_in_cm = (index * fstep_size) + fmin_d;

  // time range
  const double signal_t_range = 5000.;

  // parameterisation TF1
  TF1 VUVTiming;

  // direct path transport time
  double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
  double t_direct_min = distance_in_cm / fvuv_vgroup_max;

  // Defining the model function(s) describing the photon transportation timing vs distance
  // Getting the landau parameters from the time parametrization
  std::array<double, 3> pars_landau;
  interpolate3(pars_landau,
               fparameters[0][0],
               fparameters[2][angle_bin],
               fparameters[3][angle_bin],
               fparameters[1][angle_bin],
               distance_in_cm,
               true);
  // Deciding which time model to use (depends on the distance)
  // defining useful times for the VUV arrival time shapes
  if (distance_in_cm >= finflexion_point_distance) {
    double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
    // Set model: Landau
    VUVTiming = TF1("VUVTiming", model_far, 0, signal_t_range, 4);
    VUVTiming.SetParameters(pars_far);
  }
  else {
    // Set model: Landau + Exponential
    VUVTiming = TF1("VUVTiming", model_close, 0, signal_t_range, 7);
    // Exponential parameters
    double pars_expo[2];
    // Getting the exponential parameters from the time parametrization
    pars_expo[1] = interpolate(fparameters[4][0], fparameters[5][angle_bin], distance_in_cm, true);
    pars_expo[0] = interpolate(fparameters[4][0], fparameters[6][angle_bin], distance_in_cm, true);
    pars_expo[0] *= pars_landau[2];
    pars_expo[0] = std::log(pars_expo[0]);
    // this is to find the intersection point between the two functions:
    TF1 fint = TF1("fint", finter_d, pars_landau[0], 4 * t_direct_mean, 5);
    double parsInt[5] = {
      pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
    fint.SetParameters(parsInt);
    double t_int = fint.GetMinimumX();
    double minVal = fint.Eval(t_int);
    // the functions must intersect - output warning if they don't
    if (minVal > 0.015) {
      mf::LogWarning("PropagationTimeModel")
        << "WARNING: Parametrization of VUV light discontinuous for distance = "
        << distance_in_cm << std::endl;
      mf::LogWarning("PropagationTimeModel")
        << "WARNING: This shouldn't be happening " << std::endl;
    }
    double parsfinal[7] = {t_int,
                           pars_landau[0],
                           pars_landau[1],
                           pars_landau[2],
                           pars_expo[0],
                           pars_expo[1],
                           t_direct_min};
    VUVTiming.SetParameters(parsfinal);
  }

  // set the number of points used to sample parameterisation
  // for shorter distances, peak is sharper so more sensitive sampling required
  int fsampling;
  if (distance_in_cm < 50) fsampling = 10000;
  else if (distance_in_cm < 100) fsampling = 5000;
  else fsampling = 1000;
  VUVTiming.SetNpx(fsampling);

  // calculate max and min distance relevant to sample parameterisation
  // max
  const size_t nq_max = 1;
  double xq_max[nq_max];
  double yq_max[nq_max];
  xq_max[0] = 0.975; // include 97.5% of tail
  VUVTiming.GetQuantiles(nq_max, yq_max, xq_max);
  double max = yq_max[0];
  // min
  double min = t_direct_min;

  // store TF1 and min/max, this allows identical TF1 to be used every time sampling
  // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
  // all subsequent calls check if it has been generated previously and are ~100+ times quicker
  fVUV_timing[angle_bin][index] = VUVTiming;
  fVUV_max[angle_bin][index] = max;
  fVUV_min[angle_bin][index] = min;
}

//......................................................................
// VIS arrival times calculation functions
void
PropagationTimeModel::getVISTimes(std::vector<double>& arrivalTimes,
                                const TVector3 &ScintPoint,
                                const TVector3 &OpDetPoint)
{
  // *************************************************************************************************
  //     Calculation of earliest arrival times and corresponding unsmeared
  //     distribution
  // *************************************************************************************************

  // set plane_depth for correct TPC:
  double plane_depth;
  if (ScintPoint[0] < 0) plane_depth = -fplane_depth;
  else plane_depth = fplane_depth;

  // calculate point of reflection for shortest path
  TVector3 bounce_point(plane_depth,ScintPoint[1],ScintPoint[2]);

  // calculate distance travelled by VUV light and by vis light
  double VUVdist = (bounce_point - ScintPoint).Mag();
  double Visdist = (OpDetPoint - bounce_point).Mag();

  // calculate times taken by VUV part of path
  int angle_bin_vuv = 0; // on-axis by definition
  getVUVTimes(arrivalTimes, VUVdist, angle_bin_vuv);

  // sum parts to get total transport times times
  for (size_t i = 0; i < arrivalTimes.size(); ++i) {
    arrivalTimes[i] += Visdist / fvis_vmean;
  }

  // *************************************************************************************************
  //      Smearing of arrival time distribution
  // *************************************************************************************************
  // calculate fastest time possible
  // vis part
  double vis_time = Visdist / fvis_vmean;
  // vuv part
  double vuv_time;
  if (VUVdist < fmin_d) {
    vuv_time = VUVdist / fvuv_vgroup_max;
  }
  else {
    // find index of required parameterisation
    const size_t index = std::round((VUVdist - fmin_d) / fstep_size);
    // find shortest time
    vuv_time = fVUV_min[angle_bin_vuv][index];
  }
  // sum
  double fastest_time = vis_time + vuv_time;

  // calculate angle theta between bound_point and optical detector
  double cosine_theta = std::abs(OpDetPoint[0] - bounce_point[0]) / Visdist;
  double theta = fast_acos(cosine_theta) * 180. / CLHEP::pi;

  // determine smearing parameters using interpolation of generated points:
  // 1). tau = exponential smearing factor, varies with distance and angle
  // 2). cutoff = largest smeared time allowed, preventing excessively large
  //     times caused by exponential distance to cathode
  double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
  // angular bin
  size_t theta_bin = theta / fangle_bin_timing_vis;
  // radial distance from centre of TPC (y,z plane)
  double r = std::hypot(ScintPoint[1] - fcathode_centre[1], ScintPoint[2] - fcathode_centre[2]);

  // cut-off and tau
  // cut-off
  // interpolate in d_c for each r bin
  std::vector<double> interp_vals(fcut_off_pars[theta_bin].size(), 0.0);
  for (size_t i = 0; i < fcut_off_pars[theta_bin].size(); i++){
    interp_vals[i] = interpolate(fdistances_refl, fcut_off_pars[theta_bin][i], distance_cathode_plane, true);
  }
  // interpolate in r
  double cutoff = interpolate(fradial_distances_refl, interp_vals, r, true);

  // tau
  // interpolate in x for each r bin
  std::vector<double> interp_vals_tau(ftau_pars[theta_bin].size(), 0.0);
  for (size_t i = 0; i < ftau_pars[theta_bin].size(); i++){
    interp_vals_tau[i] = interpolate(fdistances_refl, ftau_pars[theta_bin][i], distance_cathode_plane, true);
  }
  // interpolate in r
  double tau = interpolate(fradial_distances_refl, interp_vals_tau, r, true);

  // apply smearing:
  for (size_t i = 0; i < arrivalTimes.size(); ++i) {
    double arrival_time_smeared;
    // if time is already greater than cutoff, do not apply smearing
    if (arrivalTimes[i] >= cutoff) { continue; }
    // otherwise smear
    else {
      unsigned int counter = 0;
      // loop until time generated is within cutoff limit
      // most are within single attempt, very few take more than two
      do {
        // don't attempt smearings too many times
        if (counter >= 10) {
          arrival_time_smeared = arrivalTimes[i]; // don't smear
          break;
        }
        else {
          // generate random number in appropriate range
          double x = fUniformGen.fire(0.5, 1.0);
          // apply the exponential smearing
          arrival_time_smeared =
            arrivalTimes[i] + (arrivalTimes[i] - fastest_time) * (std::pow(x, -tau) - 1);
        }
        counter++;
      } while (arrival_time_smeared > cutoff);
    }
    arrivalTimes[i] = arrival_time_smeared;
  }
}

//......................................................................
double
PropagationTimeModel::fast_acos(double x) const
{
  double negate = double(x < 0);
  x = std::abs(x);
  x -= double(x > 1.0) * (x - 1.0); // <- equivalent to min(1.0,x), but faster
  double ret = -0.0187293;
  ret = ret * x;
  ret = ret + 0.0742610;
  ret = ret * x;
  ret = ret - 0.2121144;
  ret = ret * x;
  ret = ret + 1.5707288;
  ret = ret * std::sqrt(1.0 - x);
  ret = ret - 2. * negate * ret;
  return negate * 3.14159265358979 + ret;
}

//......................................................................
// Returns interpolated value at x from parallel arrays ( xData, yData )
// Assumes that xData has at least two elements, is sorted and is strictly
// monotonic increasing boolean argument extrapolate determines behaviour
// beyond ends of array (if needed)
double
PropagationTimeModel::interpolate(const std::vector<double>& xData,
                                  const std::vector<double>& yData,
                                  double x,
                                  bool extrapolate,
                                  size_t i) const
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

//......................................................................
void
PropagationTimeModel::interpolate3(std::array<double, 3>& inter,
                         const std::vector<double>& xData,
                         const std::vector<double>& yData1,
                         const std::vector<double>& yData2,
                         const std::vector<double>& yData3,
                         double x,
                         bool extrapolate)
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

//......................................................................
double
PropagationTimeModel::finter_d(const double* x, const double* par)
{
  double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
  double y2 = TMath::Exp(par[3] + x[0] * par[4]);

  return TMath::Abs(y1 - y2);
}

//......................................................................
double
PropagationTimeModel::model_close(const double* x, const double* par)
{
  // par0 = joining point
  // par1 = Landau MPV
  // par2 = Landau width
  // par3 = normalization
  // par4 = Expo cte
  // par5 = Expo tau
  // par6 = t_min

  double y1 = par[3] * TMath::Landau(x[0], par[1], par[2]);
  double y2 = TMath::Exp(par[4] + x[0] * par[5]);
  if (x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
  if (x[0] < par[0]) y2 = 0.;

  return (y1 + y2);
}

//......................................................................
double
PropagationTimeModel::model_far(const double* x, const double* par)
{
  // par1 = Landau MPV
  // par2 = Landau width
  // par3 = normalization
  // par0 = t_min

  double y = par[3] * TMath::Landau(x[0], par[1], par[2]);
  if (x[0] <= par[0]) y = 0.;

  return y;
}
