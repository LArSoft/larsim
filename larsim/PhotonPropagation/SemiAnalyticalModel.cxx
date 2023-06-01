#include "SemiAnalyticalModel.h"
#include "larsim/PhotonPropagation/PhotonPropagationUtils.h"

// LArSoft Libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// support libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

#include <iostream>
#include <vector>

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

namespace phot {

  // constructor
  SemiAnalyticalModel::SemiAnalyticalModel(const fhicl::ParameterSet& VUVHitsParams,
                                           const fhicl::ParameterSet& VISHitsParams,
                                           const bool doReflectedLight,
                                           const bool includeAnodeReflections,
                                           const bool useXeAbsorption)
    : fGeom{*(lar::providerFrom<geo::Geometry>())}
    , fISTPC{fGeom}
    , fNTPC(fGeom.NTPC())
    , fActiveVolumes(fISTPC.extractActiveLArVolume(fGeom))
    , fcathode_centre{fGeom.TPC().GetCathodeCenter().X(),
                      fActiveVolumes[0].CenterY(),
                      fActiveVolumes[0].CenterZ()}
    , fanode_centre{fGeom.TPC().FirstPlane().GetCenter().X(),
                    fActiveVolumes[0].CenterY(),
                    fActiveVolumes[0].CenterZ()}
    , fNOpDets(fGeom.NOpDets())
    , fOpDetector(opticalDetectors())
    , fDoReflectedLight(doReflectedLight)
    , fIncludeAnodeReflections(includeAnodeReflections)
    , fUseXeAbsorption(useXeAbsorption)
  {
    mf::LogInfo("SemiAnalyticalModel") << "Initializing Semi-analytical model." << std::endl;

    // Load Gaisser-Hillas corrections for VUV semi-analytic hits
    mf::LogInfo("SemiAnalyticalModel") << "Using VUV visibility parameterization";

    fIsFlatPDCorr = VUVHitsParams.get<bool>("FlatPDCorr", false);
    fIsFlatPDCorrLat = VUVHitsParams.get<bool>("FlatPDCorrLat", false);
    fIsDomePDCorr = VUVHitsParams.get<bool>("DomePDCorr", false);
    fdelta_angulo_vuv = VUVHitsParams.get<double>("delta_angulo_vuv", 10);
    fradius = VUVHitsParams.get<double>("PMT_radius", 10.16);
    fApplyFieldCageTransparency = VUVHitsParams.get<bool>("ApplyFieldCageTransparency", false);
    fFieldCageTransparencyLateral = VUVHitsParams.get<double>("FieldCageTransparencyLateral", 1.0);
    fFieldCageTransparencyCathode = VUVHitsParams.get<double>("FieldCageTransparencyCathode", 1.0);

    if (!fIsFlatPDCorr && !fIsDomePDCorr && !fIsFlatPDCorrLat) {
      throw cet::exception("SemiAnalyticalModel")
        << "Both isFlatPDCorr/isFlatPDCorrLat and isDomePDCorr parameters are "
           "false, at least one type of parameterisation is required for the "
           "semi-analytic light simulation."
        << "\n";
    }
    if (fIsFlatPDCorr) {
      fGHvuvpars_flat = VUVHitsParams.get<std::vector<std::vector<double>>>("GH_PARS_flat");
      fborder_corr_angulo_flat = VUVHitsParams.get<std::vector<double>>("GH_border_angulo_flat");
      fborder_corr_flat = VUVHitsParams.get<std::vector<std::vector<double>>>("GH_border_flat");
    }
    if (fIsFlatPDCorrLat) {
      fGHvuvpars_flat_lateral =
        VUVHitsParams.get<std::vector<std::vector<double>>>("GH_PARS_flat_lateral");
      fborder_corr_angulo_flat_lateral =
        VUVHitsParams.get<std::vector<double>>("GH_border_angulo_flat_lateral");
      fborder_corr_flat_lateral =
        VUVHitsParams.get<std::vector<std::vector<double>>>("GH_border_flat_lateral");
    }
    if (fIsDomePDCorr) {
      fGHvuvpars_dome = VUVHitsParams.get<std::vector<std::vector<double>>>("GH_PARS_dome");
      fborder_corr_angulo_dome = VUVHitsParams.get<std::vector<double>>("GH_border_angulo_dome");
      fborder_corr_dome = VUVHitsParams.get<std::vector<std::vector<double>>>("GH_border_dome");
    }

    // Load corrections for VIS semi-analytic hits
    if (fDoReflectedLight) {
      mf::LogInfo("SemiAnalyticalModel") << "Using VIS (reflected) visibility parameterization";
      fdelta_angulo_vis = VISHitsParams.get<double>("delta_angulo_vis");

      if (fIsFlatPDCorr) {
        fvis_distances_x_flat = VISHitsParams.get<std::vector<double>>("VIS_distances_x_flat");
        fvis_distances_r_flat = VISHitsParams.get<std::vector<double>>("VIS_distances_r_flat");
        fvispars_flat =
          VISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_flat");
      }
      if (fIsDomePDCorr) {
        fvis_distances_x_dome = VISHitsParams.get<std::vector<double>>("VIS_distances_x_dome");
        fvis_distances_r_dome = VISHitsParams.get<std::vector<double>>("VIS_distances_r_dome");
        fvispars_dome =
          VISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_dome");
      }

      // set cathode plane struct for solid angle function
      fcathode_plane.h = fActiveVolumes[0].SizeY();
      fcathode_plane.w = fActiveVolumes[0].SizeZ();
      fplane_depth = std::abs(fcathode_centre[0]);
    }

    // Load corrections for Anode reflections configuration
    if (fIncludeAnodeReflections) {
      mf::LogInfo("SemiAnalyticalModel") << "Using anode reflections parameterization";
      fdelta_angulo_vis = VISHitsParams.get<double>("delta_angulo_vis");
      fAnodeReflectivity = VISHitsParams.get<double>("AnodeReflectivity");

      if (fIsFlatPDCorr) {
        fvis_distances_x_flat = VISHitsParams.get<std::vector<double>>("VIS_distances_x_flat");
        fvis_distances_r_flat = VISHitsParams.get<std::vector<double>>("VIS_distances_r_flat");
        fvispars_flat =
          VISHitsParams.get<std::vector<std::vector<std::vector<double>>>>("VIS_correction_flat");
      }

      if (fIsFlatPDCorrLat) {
        fvis_distances_x_flat_lateral =
          VISHitsParams.get<std::vector<double>>("VIS_distances_x_flat_lateral");
        fvis_distances_r_flat_lateral =
          VISHitsParams.get<std::vector<double>>("VIS_distances_r_flat_lateral");
        fvispars_flat_lateral = VISHitsParams.get<std::vector<std::vector<std::vector<double>>>>(
          "VIS_correction_flat_lateral");
      }

      // set anode plane struct for solid angle function
      fanode_plane.h = fActiveVolumes[0].SizeY();
      fanode_plane.w = fActiveVolumes[0].SizeZ();
      fanode_plane_depth = fanode_centre[0];
    }

    // set absorption length
    fvuv_absorption_length = VUVAbsorptionLength();
    mf::LogInfo("SemiAnalyticalModel")
      << "Setting absorption length to: " << fvuv_absorption_length << std::endl;

    mf::LogInfo("SemiAnalyticalModel") << "Semi-analytical model initialized." << std::endl;
  }

  double SemiAnalyticalModel::VUVAbsorptionLength() const
  {
    // determine LAr absorption length in cm
    std::map<double, double> abs_length_spectrum =
      lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
    std::vector<double> x_v, y_v;
    for (auto elem : abs_length_spectrum) {
      x_v.push_back(elem.first);
      y_v.push_back(elem.second);
    }
    double vuv_absorption_length;
    if (fUseXeAbsorption)
      vuv_absorption_length =
        interpolate(x_v, y_v, 7.1, false); // 7.1 eV: peak of Xe VUV emission spectrum
    else
      vuv_absorption_length =
        interpolate(x_v, y_v, 9.7, false); // 9.7 eV: peak of Ar VUV emission spectrum
    if (vuv_absorption_length <= 0) {
      throw cet::exception("SemiAnalyticalModel")
        << "Error: VUV Absorption Length is 0 or negative.\n";
    }
    return vuv_absorption_length;
  }

  //......................................................................
  // VUV semi-analytical model visibility calculation
  void SemiAnalyticalModel::detectedDirectVisibilities(std::vector<double>& DetectedVisibilities,
                                                       geo::Point_t const& ScintPoint) const
  {
    DetectedVisibilities.resize(fNOpDets);
    for (size_t const OpDet : util::counter(fNOpDets)) {
      if (!isOpDetInSameTPC(ScintPoint, fOpDetector[OpDet].center)) {
        DetectedVisibilities[OpDet] = 0.;
        continue;
      }

      DetectedVisibilities[OpDet] = VUVVisibility(ScintPoint, fOpDetector[OpDet]);
      ;
    }
  }

  double SemiAnalyticalModel::VUVVisibility(geo::Point_t const& ScintPoint,
                                            OpticalDetector const& opDet) const
  {
    // distance and angle between ScintPoint and OpDet center
    geo::Vector_t const relative = ScintPoint - opDet.center;
    const double distance = relative.R();
    double cosine;
    if (opDet.orientation == 2) 
      cosine = std::abs(relative.Z()) / distance;
    else if (opDet.orientation == 1)
      cosine = std::abs(relative.Y()) / distance;
    else
      cosine = std::abs(relative.X()) / distance;
    const double theta = fast_acos(cosine) * 180. / CLHEP::pi;

    double solid_angle = 0.;
    // ARAPUCAS/Bars (rectangle)
    if (opDet.type == 0) {
      // get scintillation point coordinates relative to arapuca window centre
      geo::Vector_t const abs_relative{
        std::abs(relative.X()), std::abs(relative.Y()), std::abs(relative.Z())};
      solid_angle = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_relative, opDet.orientation);
    }
    // PMTs (dome)
    else if (opDet.type == 1) {
      solid_angle = Omega_Dome_Model(distance, theta);
    }
    // PMTs (disk)
    else if (opDet.type == 2) {
      const double zy_offset = std::sqrt(relative.Y() * relative.Y() + relative.Z() * relative.Z());
      const double x_distance = std::abs(relative.X());
      solid_angle = Disk_SolidAngle(zy_offset, x_distance, fradius);
    }
    else {
      throw cet::exception("SemiAnalyticalModel")
        << "Error: Invalid optical detector shape requested - configuration "
           "error in semi-analytical model, only rectangular, dome or disk "
           "optical detectors are supported."
        << "\n";
    }

    // calculate visibility by geometric acceptance
    // accounting for solid angle and LAr absorbtion length
    double visibility_geo =
      std::exp(-1. * distance / fvuv_absorption_length) * (solid_angle / (4 * CLHEP::pi));

    // apply Gaisser-Hillas correction for Rayleigh scattering distance
    // and angular dependence offset angle bin
    const size_t j = (theta / fdelta_angulo_vuv);

    // determine GH parameters, accounting for border effects
    // radial distance from centre of detector (Y-Z standard / X-Z laterals)
    double r = 0;
    if (opDet.orientation == 2)
      r = std::hypot(ScintPoint.X() - fcathode_centre[0], ScintPoint.Y() - fcathode_centre[1]);
    else if (opDet.orientation == 1)
      r = std::hypot(ScintPoint.X() - fcathode_centre[0], ScintPoint.Z() - fcathode_centre[2]);
    else
      r = std::hypot(ScintPoint.Y() - fcathode_centre[1], ScintPoint.Z() - fcathode_centre[2]);

    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    // flat PDs
    if ((opDet.type == 0 || opDet.type == 2) && (fIsFlatPDCorr || fIsFlatPDCorrLat)) {
      if ((opDet.orientation == 1 && fIsFlatPDCorrLat) ||
          (opDet.orientation == 2 && fIsFlatPDCorrLat)) { // laterals
        pars_ini[0] = fGHvuvpars_flat_lateral[0][j];
        pars_ini[1] = fGHvuvpars_flat_lateral[1][j];
        pars_ini[2] = fGHvuvpars_flat_lateral[2][j];
        pars_ini[3] = fGHvuvpars_flat_lateral[3][j];
        s1 =
          interpolate(fborder_corr_angulo_flat_lateral, fborder_corr_flat_lateral[0], theta, true);
        s2 =
          interpolate(fborder_corr_angulo_flat_lateral, fborder_corr_flat_lateral[1], theta, true);
        s3 =
          interpolate(fborder_corr_angulo_flat_lateral, fborder_corr_flat_lateral[2], theta, true);
      }
      else if (opDet.orientation == 0 && fIsFlatPDCorr) { // cathode/anode
        pars_ini[0] = fGHvuvpars_flat[0][j];
        pars_ini[1] = fGHvuvpars_flat[1][j];
        pars_ini[2] = fGHvuvpars_flat[2][j];
        pars_ini[3] = fGHvuvpars_flat[3][j];
        s1 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[0], theta, true);
        s2 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[1], theta, true);
        s3 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[2], theta, true);
      }
      else {
        throw cet::exception("SemiAnalyticalModel")
          << "Error: flat optical detectors are found, but parameters are "
             "missing - configuration error in semi-analytical model."
          << "\n";
      }
    }
    // dome PDs
    else if (opDet.type == 1 && fIsDomePDCorr) {
      pars_ini[0] = fGHvuvpars_dome[0][j];
      pars_ini[1] = fGHvuvpars_dome[1][j];
      pars_ini[2] = fGHvuvpars_dome[2][j];
      pars_ini[3] = fGHvuvpars_dome[3][j];
      s1 = interpolate(fborder_corr_angulo_dome, fborder_corr_dome[0], theta, true);
      s2 = interpolate(fborder_corr_angulo_dome, fborder_corr_dome[1], theta, true);
      s3 = interpolate(fborder_corr_angulo_dome, fborder_corr_dome[2], theta, true);
    }
    else {
      throw cet::exception("SemiAnalyticalModel")
        << "Error: Invalid optical detector shape requested or corrections are "
           "missing - configuration error in semi-analytical model."
        << "\n";
    }

    // add border correction to parameters
    pars_ini[0] = pars_ini[0] + s1 * r;
    pars_ini[1] = pars_ini[1] + s2 * r;
    pars_ini[2] = pars_ini[2] + s3 * r;
    pars_ini[3] = pars_ini[3];

    // calculate correction
    double GH_correction = Gaisser_Hillas(distance, pars_ini);

    // determine corrected visibility of photo-detector
    return GH_correction * visibility_geo / cosine;
  }

  //......................................................................
  // VIS semi-analytical model visibility calculation
  void SemiAnalyticalModel::detectedReflectedVisibilities(
    std::vector<double>& ReflDetectedVisibilities,
    geo::Point_t const& ScintPoint,
    bool AnodeMode) const
  {
    // 1). calculate visibility of VUV photons on
    // reflective foils via solid angle + Gaisser-Hillas
    // corrections:

    // get scintpoint coords relative to centre of cathode plane and set plane
    // dimensions
    geo::Vector_t ScintPoint_relative;
    Dims plane_dimensions;
    double plane_depth;
    if (AnodeMode) {
      plane_dimensions = fanode_plane;
      plane_depth = fanode_plane_depth;
      ScintPoint_relative.SetCoordinates(std::abs(ScintPoint.X() - fanode_plane_depth),
                                         std::abs(ScintPoint.Y() - fanode_centre[1]),
                                         std::abs(ScintPoint.Z() - fanode_centre[2]));
    }
    else {
      plane_dimensions = fcathode_plane;
      plane_depth = ScintPoint.X() < 0. ? -fplane_depth : fplane_depth;
      ScintPoint_relative.SetCoordinates(std::abs(ScintPoint.X() - plane_depth),
                                         std::abs(ScintPoint.Y() - fcathode_centre[1]),
                                         std::abs(ScintPoint.Z() - fcathode_centre[2]));
    }

    // calculate solid angle of anode/cathode from the scintillation point,
    // orientation always = 0 (anode/cathode)
    double solid_angle_cathode = Rectangle_SolidAngle(plane_dimensions, ScintPoint_relative, 0);

    // calculate distance and angle between ScintPoint and hotspot
    // vast majority of hits in hotspot region directly infront of scintpoint,
    // therefore consider attenuation for this distance and on axis GH instead of
    // for the centre coordinate
    double distance_cathode = std::abs(plane_depth - ScintPoint.X());
    // calculate hits on cathode plane via geometric acceptance
    double cathode_visibility_geo = std::exp(-1. * distance_cathode / fvuv_absorption_length) *
                                    (solid_angle_cathode / (4. * CLHEP::pi));

    // determine Gaisser-Hillas correction including border effects
    // use flat correction
    double r = std::hypot(ScintPoint.Y() - fcathode_centre[1], ScintPoint.Z() - fcathode_centre[2]);
    double pars_ini[4] = {0, 0, 0, 0};
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    if (fIsFlatPDCorr) {
      pars_ini[0] = fGHvuvpars_flat[0][0];
      pars_ini[1] = fGHvuvpars_flat[1][0];
      pars_ini[2] = fGHvuvpars_flat[2][0];
      pars_ini[3] = fGHvuvpars_flat[3][0];
      s1 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[0], 0, true);
      s2 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[1], 0, true);
      s3 = interpolate(fborder_corr_angulo_flat, fborder_corr_flat[2], 0, true);
    }
    else {
      throw cet::exception("SemiAnalyticalModel")
        << "Error: flat optical detector VUV correction required for reflected "
           "semi-analytic hits. - configuration error in semi-analytical model."
        << "\n";
    }

    // add border correction
    pars_ini[0] = pars_ini[0] + s1 * r;
    pars_ini[1] = pars_ini[1] + s2 * r;
    pars_ini[2] = pars_ini[2] + s3 * r;
    pars_ini[3] = pars_ini[3];

    // calculate corrected number of hits
    double GH_correction = Gaisser_Hillas(distance_cathode, pars_ini);
    const double cathode_visibility_rec = GH_correction * cathode_visibility_geo;

    // 2). detemine visibility of each PD
    const geo::Point_t hotspot = {plane_depth, ScintPoint.Y(), ScintPoint.Z()};
    ReflDetectedVisibilities.resize(fNOpDets);
    for (size_t const OpDet : util::counter(fNOpDets)) {
      if (!isOpDetInSameTPC(ScintPoint, fOpDetector[OpDet].center)) {
        ReflDetectedVisibilities[OpDet] = 0.;
        continue;
      }

      ReflDetectedVisibilities[OpDet] =
        VISVisibility(ScintPoint, fOpDetector[OpDet], cathode_visibility_rec, hotspot, AnodeMode);
    }
  }

  double SemiAnalyticalModel::VISVisibility(geo::Point_t const& ScintPoint,
                                            OpticalDetector const& opDet,
                                            const double cathode_visibility,
                                            geo::Point_t const& hotspot,
                                            bool AnodeMode) const
  {
    // set correct plane_depth
    double plane_depth;
    if (AnodeMode)
      plane_depth = fanode_plane_depth;
    else
      plane_depth = ScintPoint.X() < 0. ? -fplane_depth : fplane_depth;

    // calculate visibility of the optical
    // detector from the hotspot using solid angle:

    geo::Vector_t const emission_relative = hotspot - opDet.center;

    // calculate distances and angles for application of corrections
    // distance from hotspot to optical detector
    const double distance_vis = emission_relative.R();
    //  angle between hotspot and optical detector
    double cosine_vis;
    if (opDet.orientation == 2) { // lateral fixed at z
      cosine_vis = std::abs(emission_relative.Z()) / distance_vis;
    }
    else if (opDet.orientation == 1) { // lateral fixed at y
      cosine_vis = std::abs(emission_relative.Y()) / distance_vis;
    }
    else { // anode/cathode (default)
      cosine_vis = std::abs(emission_relative.X()) / distance_vis;
    }
    const double theta_vis = fast_acos(cosine_vis) * 180. / CLHEP::pi;

    // calculate solid angle of optical channel
    double solid_angle_detector = 0.;
    // ARAPUCAS/Bars (rectangle)
    if (opDet.type == 0) {
      // get hotspot coordinates relative to opDet
      geo::Vector_t const abs_emission_relative{std::abs(emission_relative.X()),
                                                std::abs(emission_relative.Y()),
                                                std::abs(emission_relative.Z())};
      solid_angle_detector =
        Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_emission_relative, opDet.orientation);
    }
    // PMTS (dome)
    else if (opDet.type == 1) {
      solid_angle_detector = Omega_Dome_Model(distance_vis, theta_vis);
    }
    // PMTs (disk)
    else if (opDet.type == 2) {
      const double zy_offset = std::sqrt(emission_relative.Y() * emission_relative.Y() +
                                         emission_relative.Z() * emission_relative.Z());
      const double x_distance = std::abs(emission_relative.X());
      solid_angle_detector = Disk_SolidAngle(zy_offset, x_distance, fradius);
    }
    else {
      throw cet::exception("SemiAnalyticalModel")
        << "Error: Invalid optical detector shape requested - configuration "
           "error in semi-analytical model, only rectangular, dome or disk "
           "optical detectors are supported."
        << "\n";
    }

    // calculate number of hits via geometeric acceptance
    double visibility_geo = (solid_angle_detector / (2. * CLHEP::pi)) *
                            cathode_visibility; // 2*pi due to presence of reflective foils

    // determine correction factor, depending on PD type
    const size_t k = (theta_vis / fdelta_angulo_vis); // off-set angle bin
    double r;
    if ((opDet.orientation == 1 && fIsFlatPDCorrLat) ||
        (opDet.orientation == 2 && fIsFlatPDCorrLat))
      r = std::hypot(ScintPoint.X() - fcathode_centre[0], ScintPoint.Z() - fcathode_centre[2]);
    else
      r = std::hypot(ScintPoint.Y() - fcathode_centre[1], ScintPoint.Z() - fcathode_centre[2]);
    double d_c = std::abs(ScintPoint.X() - plane_depth); // distance to cathode
    double border_correction = 0;
    // flat PDs
    if ((opDet.type == 0 || opDet.type == 2) && (fIsFlatPDCorr || fIsFlatPDCorrLat)) {
      // cathode/anode case
      if (opDet.orientation == 0 && fIsFlatPDCorr) {

        border_correction =
          interpolate2(fvis_distances_x_flat, fvis_distances_r_flat, fvispars_flat, d_c, r, k);
      }
      // laterals case
      else if ((opDet.orientation == 1 && fIsFlatPDCorrLat) ||
               (opDet.orientation == 2 && fIsFlatPDCorrLat)) {
        border_correction = interpolate2(fvis_distances_x_flat_lateral,
                                         fvis_distances_r_flat_lateral,
                                         fvispars_flat_lateral,
                                         d_c,
                                         r,
                                         k);
      }
      else {
        throw cet::exception("SemiAnalyticalModel")
          << "Error: Invalid optical detector shape requested or corrections "
             "are missing - configuration error in semi-analytical model."
          << "\n";
      }
    }
    // dome PDs
    else if (opDet.type == 1 && fIsDomePDCorr) {
      border_correction =
        interpolate2(fvis_distances_x_dome, fvis_distances_r_dome, fvispars_dome, d_c, r, k);
    }
    else {
      throw cet::exception("SemiAnalyticalModel")
        << "Error: Invalid optical detector shape requested or corrections are "
           "missing - configuration error in semi-analytical model."
        << "\n";
    }

    // apply anode reflectivity factor
    if (AnodeMode) border_correction *= fAnodeReflectivity;

    return border_correction * visibility_geo / cosine_vis;
  }

  //......................................................................
  // Gaisser-Hillas function definition
  double SemiAnalyticalModel::Gaisser_Hillas(const double x, const double* par) const
  {
    double X_mu_0 = par[3];
    double Normalization = par[0];
    double Diff = par[1] - X_mu_0;
    double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
    double Exponential = std::exp((par[1] - x) / par[2]);

    return (Normalization * Term * Exponential);
  }

  //......................................................................
  // solid angle of circular aperture
  double SemiAnalyticalModel::Disk_SolidAngle(const double d, const double h, const double b) const
  {
    if (b <= 0. || d < 0. || h <= 0.) return 0.;
    const double leg2 = (b + d) * (b + d);
    const double aa = std::sqrt(h * h / (h * h + leg2));
    if (isApproximatelyZero(d)) { return 2. * CLHEP::pi * (1. - aa); }
    double bb = 2. * std::sqrt(b * d / (h * h + leg2));
    double cc = 4. * b * d / leg2;

    if (isDefinitelyGreaterThan(d, b)) {
      try {
        return 2. * aa *
               (std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()) -
                boost::math::ellint_1(bb, noLDoublePromote()));
      }
      catch (std::domain_error& e) {
        if (isApproximatelyEqual(d, b, 1e-9)) {
          mf::LogWarning("SemiAnalyticalModel")
            << "Elliptic Integral in Disk_SolidAngle() given parameters "
               "outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("SemiAnalyticalModel")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters "
               "outside domain.\n"
            << "\nbb: " << bb << "\ncc: " << cc << "Exception message: " << e.what();
          return 0.;
        }
      }
    }
    if (isDefinitelyLessThan(d, b)) {
      try {
        return 2. * CLHEP::pi -
               2. * aa *
                 (boost::math::ellint_1(bb, noLDoublePromote()) +
                  std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()));
      }
      catch (std::domain_error& e) {
        if (isApproximatelyEqual(d, b, 1e-9)) {
          mf::LogWarning("SemiAnalyticalModel")
            << "Elliptic Integral in Disk_SolidAngle() given parameters "
               "outside domain."
            << "\nbb: " << bb << "\ncc: " << cc << "\nException message: " << e.what()
            << "\nRelax condition and carry on.";
          return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
        }
        else {
          mf::LogError("SemiAnalyticalModel")
            << "Elliptic Integral inside Disk_SolidAngle() given parameters "
               "outside domain.\n"
            << "\nbb: " << bb << "\ncc: " << cc << "Exception message: " << e.what();
          return 0.;
        }
      }
    }
    if (isApproximatelyEqual(d, b)) {
      return CLHEP::pi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
    }
    return 0.;
  }

  //......................................................................
  // solid angle of rectangular aperture
  double SemiAnalyticalModel::Rectangle_SolidAngle(const double a,
                                                   const double b,
                                                   const double d) const
  {
    double aa = a / (2. * d);
    double bb = b / (2. * d);
    double aux = (1. + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
    return 4. * fast_acos(std::sqrt(aux));
  }

  double SemiAnalyticalModel::Rectangle_SolidAngle(Dims const& o,
                                                   geo::Vector_t const& v,
                                                   int OpDetOrientation) const
  {
    // v is the position of the track segment with respect to
    // the center position of the arapuca window

    // solid angle calculation depends on orientation of PD, set correct distances
    // to use
    double d1;
    double d2;
    if (OpDetOrientation == 2) {
      // lateral PD, arapuca plane fixed in z direction
      d1 = std::abs(v.X());
      d2 = std::abs(v.Z());
    }
    else if (OpDetOrientation == 1) {
      // lateral PD, arapuca plane fixed in y direction
      d1 = std::abs(v.X());
      d2 = std::abs(v.Y());
    }
    else {
      // anode/cathode PD, arapuca plane fixed in x direction [default]
      d1 = std::abs(v.Y());
      d2 = std::abs(v.X());
    }
    // arapuca plane fixed in x direction
    if (isApproximatelyZero(d1) && isApproximatelyZero(v.Z())) {
      return Rectangle_SolidAngle(o.h, o.w, d2);
    }
    if (isDefinitelyGreaterThan(d1, o.h * .5) &&
        isDefinitelyGreaterThan(std::abs(v.Z()), o.w * .5)) {
      double A = d1 - o.h * .5;
      double B = std::abs(v.Z()) - o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (A + o.h), 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * A, 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * (A + o.h), 2. * B, d2) +
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    if ((d1 <= o.h * .5) && (std::abs(v.Z()) <= o.w * .5)) {
      double A = -d1 + o.h * .5;
      double B = -std::abs(v.Z()) + o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (o.h - A), 2. * (o.w - B), d2) +
                          Rectangle_SolidAngle(2. * A, 2. * (o.w - B), d2) +
                          Rectangle_SolidAngle(2. * (o.h - A), 2. * B, d2) +
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    if (isDefinitelyGreaterThan(d1, o.h * .5) && (std::abs(v.Z()) <= o.w * .5)) {
      double A = d1 - o.h * .5;
      double B = -std::abs(v.Z()) + o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (A + o.h), 2. * (o.w - B), d2) -
                          Rectangle_SolidAngle(2. * A, 2. * (o.w - B), d2) +
                          Rectangle_SolidAngle(2. * (A + o.h), 2. * B, d2) -
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }
    if ((d1 <= o.h * .5) && isDefinitelyGreaterThan(std::abs(v.Z()), o.w * .5)) {
      double A = -d1 + o.h * .5;
      double B = std::abs(v.Z()) - o.w * .5;
      double to_return = (Rectangle_SolidAngle(2. * (o.h - A), 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * (o.h - A), 2. * B, d2) +
                          Rectangle_SolidAngle(2. * A, 2. * (B + o.w), d2) -
                          Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
                         .25;
      return to_return;
    }

    return 0.;
  }

  //......................................................................
  // solid angle of dome aperture
  double SemiAnalyticalModel::Omega_Dome_Model(const double distance, const double theta) const
  {
    // this function calculates the solid angle of a semi-sphere of radius b,
    // as a correction to the analytic formula of the on-axix solid angle,
    // as we move off-axis an angle theta. We have used 9-angular bins
    // with delta_theta width.

    // par0 = Radius correction close
    // par1 = Radius correction far
    // par2 = breaking distance betwween "close" and "far"

    constexpr double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221};
    constexpr double par1[9] = {
      0., 0., 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
    const double delta_theta = 10.; // TODO: should this be fdelta_angulo_vuv?
    int j = int(theta / delta_theta);
    // PMT radius
    const double b = fradius; // cm
    // distance form which the model parameters break (empirical value)
    const double d_break = 5 * b; // par2

    if (distance >= d_break) {
      double R_apparent_far = b - par1[j];
      double ratio_square = (R_apparent_far * R_apparent_far) / (distance * distance);
      return (2 * CLHEP::pi * (1 - std::sqrt(1 - ratio_square)));
    }
    else {
      double R_apparent_close = b - par0[j];
      double ratio_square = (R_apparent_close * R_apparent_close) / (distance * distance);
      return (2 * CLHEP::pi * (1 - std::sqrt(1 - ratio_square)));
    }
  }

  //......................................................................
  // checks photo-detector is in same TPC/argon volume as scintillation
  bool SemiAnalyticalModel::isOpDetInSameTPC(geo::Point_t const& ScintPoint,
                                             geo::Point_t const& OpDetPoint) const
  {
    // method working for SBND, uBooNE, DUNE-HD 1x2x6 and DUNE-VD 1x8x6
    // will need to be replaced to work in full DUNE geometry, ICARUS geometry
    // check x coordinate has same sign or is close to zero, otherwise return
    // false
    if (((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) && std::abs(OpDetPoint.X()) > 10. &&
        fNTPC == 2) { // TODO: replace with geometry service
      return false;
    }
    return true;
  }

  std::vector<SemiAnalyticalModel::OpticalDetector> SemiAnalyticalModel::opticalDetectors() const
  {
    std::vector<SemiAnalyticalModel::OpticalDetector> opticalDetector;
    for (size_t const i : util::counter(fNOpDets)) {
      geo::OpDetGeo const& opDet = fGeom.OpDetGeoFromOpDet(i);
      geo::Point_t center = opDet.GetCenter();
      double height;
      double length;
      int type;
      int orientation;
      if (opDet.isSphere()) { // dome PMTs
        type = 1;             // dome
        orientation = 0;      // anode/cathode (default)
        length = -1;
        height = -1;
      }
      else if (opDet.isBar()) {
        type = 0; // (X)Arapucas/Bars
        // determine orientation to get correction OpDet dimensions
        length = opDet.Length();
        height = opDet.Height();
        if (opDet.Width() > opDet.Length()) { // laterals along x-y plane, Z dimension smallest
          orientation = 2;
          length = opDet.Width();
        }
        else if (opDet.Width() > opDet.Height()) { // laterals along x-z plane, Y dimension smallest
          orientation = 1;
          height = opDet.Width();
        }
        else { // anode/cathode (default), X dimension smallest
          orientation = 0;
        }
      }
      else {
        type = 2;        // disk PMTs
        orientation = 0; // anode/cathode (default)
        length = -1;
        height = -1;
      }
      opticalDetector.emplace_back(
        SemiAnalyticalModel::OpticalDetector{height, length, center, type, orientation});
    }
    return opticalDetector;
  }

} // namespace phot
