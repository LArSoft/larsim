#ifndef SEMIANALYTICALMODEL_H
#define SEMIANALYTICALMODEL_H

// SemiAnalyticalModel
//  - fast optical simulation using semi-analytical model
//  - contains functions to calculate the number of direct and reflected photons
//  incident
//    each photo-detector, along with the necessary ultility functions (geometry
//    calculations etc.)
//  - full description of model: Eur. Phys. J. C 81, 349 (2021)

// Nov 2021 by P. Green

// LArSoft Libraries
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larsim/IonizationScintillation/ISTPC.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"

// fhicl
#include "fhiclcpp/fwd.h"

#include "TVector3.h"

#include "boost/math/policies/policy.hpp"

#include <map>
#include <vector>
#include <memory>

// Define a new policy *not* internally promoting RealType to double:
typedef boost::math::policies::policy<boost::math::policies::promote_double<false>>
  noLDoublePromote;

namespace phot {
  class SemiAnalyticalModel {

  public:
    // constructor
    SemiAnalyticalModel(const fhicl::ParameterSet& VUVHits,
                        const fhicl::ParameterSet& VISHits,
                        const std::shared_ptr<OpticalPath> &OpticalPath,
                        const bool doReflectedLight = false,
                        const bool includeAnodeReflections = false,
                        const bool useXeAbsorption = false);

    // direct / VUV light
    void detectedDirectVisibilities(std::vector<double>& DetectedVisibilities,
                                    geo::Point_t const& ScintPoint) const;

    // reflected / visible light
    void detectedReflectedVisibilities(std::vector<double>& ReflDetectedVisibilities,
                                       geo::Point_t const& ScintPoint,
                                       bool AnodeMode = false) const;

  private:
    double VUVAbsorptionLength() const;

    // structure for rectangular solid angle calculation
    struct Dims {
      double h, w; // height, width
    };

    // structure for optical detector information
    struct OpticalDetector {
      double h; // height
      double w; // width
      geo::Point_t center;
      int type;
      int orientation;
    };

    // direct light photo-detector visibility calculation
    double VUVVisibility(geo::Point_t const& ScintPoint, OpticalDetector const& opDet) const;

    // reflected light photo-detector visibility calculation
    double VISVisibility(geo::Point_t const& ScintPoint,
                         OpticalDetector const& opDet,
                         const double cathode_visibility,
                         geo::Point_t const& hotspot,
                         bool AnodeMode = false) const;

    // Gaisser-Hillas
    double Gaisser_Hillas(const double x, const double* par) const;

    // solid angle calculations
    // rectangular aperture
    double Rectangle_SolidAngle(const double a, const double b, const double d) const;
    double Rectangle_SolidAngle(Dims const& o,
                                geo::Vector_t const& v,
                                const int OpDetOrientation) const;
    // circular aperture
    double Disk_SolidAngle(const double d, const double h, const double b) const;
    // dome aperture calculation
    double Omega_Dome_Model(const double distance, const double theta) const;

    bool isOpDetInSameTPC(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) const;
    std::vector<OpticalDetector> opticalDetectors() const;

    // geometry properties
    geo::GeometryCore const& fGeom;
    geo::WireReadoutGeom const& fChannelMap;
    const larg4::ISTPC fISTPC;
    const int fNTPC;
    const std::vector<geo::BoxBoundedGeo> fActiveVolumes;
    const TVector3 fcathode_centre, fanode_centre;
    double fplane_depth, fanode_plane_depth;
    double fDriftDistance;

    // photodetector geometry properties
    const size_t fNOpDets;
    const std::vector<OpticalDetector> fOpDetector;
    double fradius;
    Dims fcathode_plane;
    Dims fanode_plane;

    // For VUV semi-analytic hits
    double fdelta_angulo_vuv;
    // flat PDs
    bool fIsFlatPDCorr;
    std::vector<std::vector<double>> fGHvuvpars_flat;
    std::vector<double> fborder_corr_angulo_flat;
    std::vector<std::vector<double>> fborder_corr_flat;
    // lateral PDs
    bool fIsFlatPDCorrLat;
    std::vector<std::vector<double>> fGHvuvpars_flat_lateral;
    std::vector<double> fborder_corr_angulo_flat_lateral;
    std::vector<std::vector<double>> fborder_corr_flat_lateral;

    // dome PDs
    bool fIsDomePDCorr;
    std::vector<std::vector<double>> fGHvuvpars_dome;
    std::vector<double> fborder_corr_angulo_dome;
    std::vector<std::vector<double>> fborder_corr_dome;
    // Field cage scaling
    bool fApplyFieldCageTransparency;
    double fFieldCageTransparencyLateral;
    double fFieldCageTransparencyCathode;

    // For VIS semi-analytic hits
    const bool fDoReflectedLight;
    const bool fIncludeAnodeReflections;
    // correction parameters for VIS Nhits estimation
    double fdelta_angulo_vis;
    double fAnodeReflectivity;
    // flat PDs
    std::vector<double> fvis_distances_x_flat;
    std::vector<double> fvis_distances_r_flat;
    std::vector<std::vector<std::vector<double>>> fvispars_flat;
    // lateral PDs
    std::vector<double> fvis_distances_x_flat_lateral;
    std::vector<double> fvis_distances_r_flat_lateral;
    std::vector<std::vector<std::vector<double>>> fvispars_flat_lateral;
    // dome PDs
    std::vector<double> fvis_distances_x_dome;
    std::vector<double> fvis_distances_r_dome;
    std::vector<std::vector<std::vector<double>>> fvispars_dome;

    // absorption length
    const bool fUseXeAbsorption;
    double fvuv_absorption_length;

    // maximum distance
    double fMaxPDDistance;

    // flag to apply border corrections for vertical direction only
    bool fVerticalBorderCorrectionMode;

    // Tool to to determine visibility of optical detectors from scintillation emission points
    std::shared_ptr<OpticalPath> fOpticalPath;
  };

} // namespace phot

#endif
