#include "PropagationTimeModel.h"
#include "larsim/PhotonPropagation/PhotonPropagationUtils.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"

// support libraries
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

#include <iostream>

namespace {
  //......................................................................
  double finter_d(const double* x, const double* par)
  {
    double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
    double y2 = TMath::Exp(par[3] + x[0] * par[4]);
    return TMath::Abs(y1 - y2);
  }

  //......................................................................
  double model_close(const double* x, const double* par)
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
  double model_far(const double* x, const double* par)
  {
    // par1 = Landau MPV
    // par2 = Landau width
    // par3 = normalization
    // par0 = t_min
    if (x[0] <= par[0]) return 0.;
    return par[3] * TMath::Landau(x[0], par[1], par[2]);
  }
}

namespace phot {

  PropagationTimeModel::PropagationTimeModel(const fhicl::ParameterSet& VUVTimingParams,
                                             const fhicl::ParameterSet& VISTimingParams,
                                             CLHEP::HepRandomEngine& scintTimeEngine,
                                             const bool doReflectedLight,
                                             const bool GeoPropTimeOnly)
    : fGeoPropTimeOnly(GeoPropTimeOnly)
    , fUniformGen(scintTimeEngine)
    , fGeom(*(lar::providerFrom<geo::Geometry>()))
    , fplane_depth(std::abs(fGeom.TPC().GetCathodeCenter().X()))
    , fOpDetCenter(opDetCenters())
    , fOpDetOrientation(opDetOrientations())
  {
    mf::LogInfo("PropagationTimeModel")
      << "Initializing Photon propagation time model." << std::endl;
    if (!fGeoPropTimeOnly) {
      // Direct / VUV
      mf::LogInfo("PropagationTimeModel") << "Using VUV timing parameterization";
      fstep_size = VUVTimingParams.get<double>("step_size");
      fmin_d = VUVTimingParams.get<double>("min_d");
      fvuv_vgroup_mean = VUVTimingParams.get<double>("vuv_vgroup_mean");
      fvuv_vgroup_max = VUVTimingParams.get<double>("vuv_vgroup_max");
      finflexion_point_distance = VUVTimingParams.get<double>("inflexion_point_distance");
      fangle_bin_timing_vuv = VUVTimingParams.get<double>("angle_bin_timing_vuv");
      generateVUVParams(VUVTimingParams, scintTimeEngine);

      // Reflected / Visible
      if (doReflectedLight) {
        mf::LogInfo("PropagationTimeModel") << "Using VIS (reflected) timing parameterization";
        fdistances_refl = VISTimingParams.get<std::vector<double>>("Distances_refl");
        fradial_distances_refl = VISTimingParams.get<std::vector<double>>("Distances_radial_refl");
        fcut_off_pars =
          VISTimingParams.get<std::vector<std::vector<std::vector<double>>>>("Cut_off");
        ftau_pars = VISTimingParams.get<std::vector<std::vector<std::vector<double>>>>("Tau");
        fvis_vmean = VISTimingParams.get<double>("vis_vmean");
        fangle_bin_timing_vis = VISTimingParams.get<double>("angle_bin_timing_vis");
      }
    }

    if (fGeoPropTimeOnly) {
      mf::LogInfo("PropagationTimeModel") << "Using geometric VUV time propagation";
      fvuv_vgroup_mean = VUVTimingParams.get<double>("vuv_vgroup_mean");
    }

    mf::LogInfo("PropagationTimeModel") << "Photon propagation time model initalized." << std::endl;
  }

  //......................................................................
  // Propagation time calculation function
  void PropagationTimeModel::propagationTime(std::vector<double>& arrivalTimes,
                                             const geo::Point_t& x0,
                                             const size_t OpChannel,
                                             const bool Reflected)
  {
    if (!fGeoPropTimeOnly) {
      // Get VUV photons transport time distribution from the parametrization
      geo::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
      if (!Reflected) {
        double distance =
          std::hypot(x0.X() - opDetCenter.X(), x0.Y() - opDetCenter.Y(), x0.Z() - opDetCenter.Z());
        double cosine;
        if (fOpDetOrientation[OpChannel] == 1)
          cosine = std::abs(x0.Y() - opDetCenter.Y()) / distance;
        else
          cosine = std::abs(x0.X() - opDetCenter.X()) / distance;

        double theta = fast_acos(cosine) * 180. / CLHEP::pi;
        int angle_bin = theta / fangle_bin_timing_vuv;
        getVUVTimes(arrivalTimes, distance, angle_bin); // in ns
      }
      else {
        getVISTimes(arrivalTimes, x0, opDetCenter); // in ns
      }
    }
    else if (fGeoPropTimeOnly && !Reflected) {
      // Get VUV photons arrival time geometrically
      geo::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
      double distance =
        std::hypot(x0.X() - opDetCenter.X(), x0.Y() - opDetCenter.Y(), x0.Z() - opDetCenter.Z());
      getVUVTimesGeo(arrivalTimes, distance); // in ns
    }
    else {
      throw cet::exception("PropagationTimeModel") << "Propagation time model not found.";
    }
  }

  //......................................................................
  // VUV propagation times calculation function
  void PropagationTimeModel::getVUVTimes(std::vector<double>& arrivalTimes,
                                         const double distance,
                                         const size_t angle_bin)
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
        arrivalTimes[i] = fVUVTimingGen[angle_bin][index].fire() *
                            (fVUV_max[angle_bin][index] - fVUV_min[angle_bin][index]) +
                          fVUV_min[angle_bin][index];
      }
    }
  }

  //......................................................................
  // VUV arrival times calculation function
  // - pure geometric approximation for use in Xenon doped scenarios
  void PropagationTimeModel::getVUVTimesGeo(std::vector<double>& arrivalTimes,
                                            const double distance) const
  {
    // times are fixed shift i.e. direct path only
    double t_prop_correction = distance / fvuv_vgroup_mean;
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] = t_prop_correction;
    }
  }

  //......................................................................
  // VUV propagation times parameterization generation function
  void PropagationTimeModel::generateVUVParams(const fhicl::ParameterSet& VUVTimingParams,
                                               CLHEP::HepRandomEngine& scintTimeEngine)
  {
    mf::LogInfo("PropagationTimeModel") << "Generating VUV parameters";
    double max_d = VUVTimingParams.get<double>("max_d");
    const size_t num_params =
      (max_d - fmin_d) / fstep_size; // for d < fmin_d, no parameterisaton, a
                                     // delta function is used instead
    const size_t num_angles = std::round(90. / fangle_bin_timing_vuv);

    // initialise vectors to contain range parameterisations
    double dummy[1] = {1.};
    fVUVTimingGen = std::vector(num_angles, std::vector(num_params, CLHEP::RandGeneral(dummy, 1)));
    fVUV_max = std::vector(num_angles, std::vector(num_params, 0.0));
    fVUV_min = std::vector(num_angles, std::vector(num_params, 0.0));

    std::vector<std::vector<double>> parameters[7];
    parameters[0] = std::vector(1, VUVTimingParams.get<std::vector<double>>("Distances_landau"));
    parameters[1] = VUVTimingParams.get<std::vector<std::vector<double>>>("Norm_over_entries");
    parameters[2] = VUVTimingParams.get<std::vector<std::vector<double>>>("Mpv");
    parameters[3] = VUVTimingParams.get<std::vector<std::vector<double>>>("Width");
    parameters[4] = std::vector(1, VUVTimingParams.get<std::vector<double>>("Distances_exp"));
    parameters[5] = VUVTimingParams.get<std::vector<std::vector<double>>>("Slope");
    parameters[6] = VUVTimingParams.get<std::vector<std::vector<double>>>("Expo_over_Landau_norm");

    // time range
    const double signal_t_range = 5000.;

    for (size_t index = 0; index < num_params; ++index) {
      double distance_in_cm = (index * fstep_size) + fmin_d;

      // direct path transport time
      double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
      double t_direct_min = distance_in_cm / fvuv_vgroup_max;

      // number of sampling points, for shorter distances, peak is
      // sharper so more sensitive sampling required
      int sampling;
      if (distance_in_cm < 2. * fmin_d)
        sampling = 10000;
      else if (distance_in_cm < 4. * fmin_d)
        sampling = 5000;
      else
        sampling = 1000;

      for (size_t angle_bin = 0; angle_bin < num_angles; ++angle_bin) {
        // Defining the model function(s) describing the photon
        // transportation timing vs distance.
        // Getting the landau parameters from the time parametrization
        std::array<double, 3> pars_landau;
        interpolate3(pars_landau,
                     parameters[0][0],
                     parameters[2][angle_bin],
                     parameters[3][angle_bin],
                     parameters[1][angle_bin],
                     distance_in_cm,
                     true);
        TF1 VUVTiming;
        // Deciding which time model to use (depends on the distance)
        // defining useful times for the VUV arrival time shapes
        if (distance_in_cm >= finflexion_point_distance) {
          // Set model: Landau
          double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
          VUVTiming = TF1("VUVTiming", model_far, 0., signal_t_range, 4);
          VUVTiming.SetParameters(pars_far);
        }
        else {
          // Set model: Landau + Exponential
          double pars_expo[2];
          // Getting the exponential parameters from the time parametrization
          pars_expo[1] =
            interpolate(parameters[4][0], parameters[5][angle_bin], distance_in_cm, true);
          pars_expo[0] =
            interpolate(parameters[4][0], parameters[6][angle_bin], distance_in_cm, true);
          pars_expo[0] *= pars_landau[2];
          pars_expo[0] = std::log(pars_expo[0]);
          // this is to find the intersection point between the two functions:
          TF1 fint = TF1("fint", finter_d, pars_landau[0], 4. * t_direct_mean, 5);
          double parsInt[5] = {
            pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
          fint.SetParameters(parsInt);
          double t_int = fint.GetMinimumX();
          double minVal = fint.Eval(t_int);
          // the functions must intersect - output warning if they don't
          if (minVal > 0.015) {
            mf::LogWarning("PropagationTimeModel")
              << "WARNING: Parametrization of VUV light discontinuous for "
                 "distance = "
              << distance_in_cm << "\n"
              << "This shouldn't be happening " << std::endl;
          }
          double parsfinal[7] = {t_int,
                                 pars_landau[0],
                                 pars_landau[1],
                                 pars_landau[2],
                                 pars_expo[0],
                                 pars_expo[1],
                                 t_direct_min};
          VUVTiming = TF1("VUVTiming", model_close, 0., signal_t_range, 7);
          VUVTiming.SetParameters(parsfinal);
        }

        // store min/max, necessary to transform to the domain of the
        // original distribution
        // const size_t nq_max = 1;
        double xq_max[1] = {0.975}; // include 97.5% of tail
        double yq_max[1];
        VUVTiming.SetNpx(sampling);
        VUVTiming.GetQuantiles(1, yq_max, xq_max);
        double max = yq_max[0];
        double min = t_direct_min;
        VUVTiming.SetRange(min, max);
        fVUV_max[angle_bin][index] = max;
        fVUV_min[angle_bin][index] = min;

        // create the distributions that represent the parametrised timing,
        // and use those to create a RNG that samples said distributions
        auto hh = (TH1D*)VUVTiming.GetHistogram();
        std::vector<double> vuv_timings(sampling, 0.);
        for (int i = 0; i < sampling; ++i)
          vuv_timings[i] = hh->GetBinContent(i + 1);
        fVUVTimingGen[angle_bin][index] =
          CLHEP::RandGeneral(scintTimeEngine, vuv_timings.data(), vuv_timings.size());
      } // index < num_params
    }   // angle_bin < num_angles
  }

  //......................................................................
  // VIS arrival times calculation functions
  void PropagationTimeModel::getVISTimes(std::vector<double>& arrivalTimes,
                                         const geo::Point_t& scintPoint,
                                         const geo::Point_t& opDetPoint)
  {
    // ***************************************************************************
    //     Calculation of earliest arrival times and corresponding unsmeared
    //     distribution
    // ***************************************************************************

    // set plane_depth for correct TPC:
    double plane_depth;
    if (scintPoint.X() < 0.)
      plane_depth = -fplane_depth;
    else
      plane_depth = fplane_depth;

    // calculate point of reflection for shortest path
    geo::Point_t bounce_point(plane_depth, scintPoint.Y(), scintPoint.Z());

    // calculate distance travelled by VUV light and by vis light
    double VUVdist = std::sqrt((bounce_point - scintPoint).Mag2());
    double Visdist = std::sqrt((opDetPoint - bounce_point).Mag2());

    // calculate times taken by VUV part of path
    int angle_bin_vuv = 0; // on-axis by definition
    getVUVTimes(arrivalTimes, VUVdist, angle_bin_vuv);

    // sum parts to get total transport times times
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] += Visdist / fvis_vmean;
    }

    // ***************************************************************************
    //      Smearing of arrival time distribution
    // ***************************************************************************
    // calculate fastest time possible
    // vis part
    double vis_time = Visdist / fvis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < fmin_d) { vuv_time = VUVdist / fvuv_vgroup_max; }
    else {
      // find index of required parameterisation
      const size_t index = std::round((VUVdist - fmin_d) / fstep_size);
      // find shortest time
      vuv_time = fVUV_min[angle_bin_vuv][index];
    }
    double fastest_time = vis_time + vuv_time;

    // calculate angle theta between bound_point and optical detector
    double cosine_theta = std::abs(opDetPoint.X() - bounce_point.X()) / Visdist;
    double theta = fast_acos(cosine_theta) * 180. / CLHEP::pi;

    // determine smearing parameters using interpolation of generated points:
    // 1). tau = exponential smearing factor, varies with distance and angle
    // 2). cutoff = largest smeared time allowed, preventing excessively large
    //     times caused by exponential distance to cathode
    double distance_cathode_plane = std::abs(plane_depth - scintPoint.X());
    // angular bin
    size_t theta_bin = theta / fangle_bin_timing_vis;
    // radial distance from centre of TPC (y,z plane)
    double r =
      std::hypot(scintPoint.Y() - fcathode_centre.Y(), scintPoint.Z() - fcathode_centre.Z());

    // cut-off and tau
    // cut-off
    // interpolate in d_c for each r bin
    std::vector<double> interp_vals(fcut_off_pars[theta_bin].size(), 0.0);
    for (size_t i = 0; i < fcut_off_pars[theta_bin].size(); i++) {
      interp_vals[i] =
        interpolate(fdistances_refl, fcut_off_pars[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double cutoff = interpolate(fradial_distances_refl, interp_vals, r, true);

    // tau
    // interpolate in x for each r bin
    std::vector<double> interp_vals_tau(ftau_pars[theta_bin].size(), 0.0);
    for (size_t i = 0; i < ftau_pars[theta_bin].size(); i++) {
      interp_vals_tau[i] =
        interpolate(fdistances_refl, ftau_pars[theta_bin][i], distance_cathode_plane, true);
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

  geo::Point_t PropagationTimeModel::cathodeCentre() const
  {
    larg4::ISTPC is_tpc = larg4::ISTPC{fGeom};
    std::vector<geo::BoxBoundedGeo> activeVolumes = is_tpc.extractActiveLArVolume(fGeom);
    geo::Point_t cathode_centre = {
      fGeom.TPC().GetCathodeCenter().X(), activeVolumes[0].CenterY(), activeVolumes[0].CenterZ()};
    return cathode_centre;
  }

  std::vector<geo::Point_t> PropagationTimeModel::opDetCenters() const
  {
    std::vector<geo::Point_t> opDetCenter;
    for (size_t const i : util::counter(fGeom.NOpDets())) {
      geo::OpDetGeo const& opDet = fGeom.OpDetGeoFromOpDet(i);
      opDetCenter.push_back(opDet.GetCenter());
    }
    return opDetCenter;
  }

  std::vector<int> PropagationTimeModel::opDetOrientations() const
  {
    std::vector<int> opDetOrientation;
    for (size_t const i : util::counter(fGeom.NOpDets())) {
      geo::OpDetGeo const& opDet = fGeom.OpDetGeoFromOpDet(i);
      if (opDet.isSphere()) { // dome PMTs
        opDetOrientation.push_back(0);
      }
      else if (opDet.isBar()) {
        // determine orientation to get correction OpDet dimensions
        if (opDet.Width() > opDet.Height()) { // laterals, Y dimension smallest
          opDetOrientation.push_back(1);
        }
        else { // anode/cathode (default), X dimension smallest
          opDetOrientation.push_back(0);
        }
      }
      else {
        opDetOrientation.push_back(0);
      }
    }
    return opDetOrientation;
  }

} // namespace phot
