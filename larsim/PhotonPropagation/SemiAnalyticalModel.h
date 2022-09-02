#ifndef SEMIANALYTICALMODEL_H
#define SEMIANALYTICALMODEL_H

// SemiAnalyticalModel
//  - fast optical simulation using semi-analytical model
//  - contains functions to calculate the number of direct and reflected photons incident
//    each photo-detector, along with the necessary ultility functions (geometry calculations etc.)
//  - full description of model: Eur. Phys. J. C 81, 349 (2021)

// Nov 2021 by P. Green

// LArSoft Libraries
#include "larsim/IonizationScintillation/ISTPC.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// fhicl
#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"

#include "boost/math/policies/policy.hpp"

#include <map>
#include <vector>


// Define a new policy *not* internally promoting RealType to double:
typedef boost::math::policies::policy<boost::math::policies::promote_double<false>> noLDoublePromote;


class SemiAnalyticalModel {

public:

    // constructor
    SemiAnalyticalModel(fhicl::ParameterSet VUVHits, fhicl::ParameterSet VISHits, bool doReflectedLight = false, bool includeAnodeReflections = false);

    // direct / VUV light
    void detectedDirectVisibilities(std::vector<double>& DetectedVisibilities,
                                    geo::Point_t const& ScintPoint) const;

    // reflected / visible light
    void detectedReflectedVisibilities(std::vector<double>& ReflDetectedVisibilities,
                                       geo::Point_t const& ScintPoint,
                                       bool AnodeMode = false) const;

private:

    // parameter and geometry initialization
    void Initialization();

    int VUVAbsorptionLength() const;

    // structure for rectangular solid angle calculation
    struct Dims {
      double h, w; // height, width
    };

    // structure for optical detector information
    struct OpticalDetector {
      double h; // height
      double w; // width
      geo::Point_t OpDetPoint;
      int type;
      int orientation;
    };

    // direct light photo-detector visibility calculation
    double VUVVisibility(geo::Point_t const& ScintPoint, OpticalDetector const& opDet) const;

    // reflected light photo-detector visibility calculation
    double VISVisibility(geo::Point_t const& ScintPoint, OpticalDetector const& opDet,
                         const double cathode_visibility, geo::Point_t const& hotspot,
                         bool AnodeMode = false) const;

    // Gaisser-Hillas
    double Gaisser_Hillas(const double x, const double* par) const;

    // solid angle calculations
    // rectangular aperture
    double Rectangle_SolidAngle(const double a, const double b, const double d) const;
    double Rectangle_SolidAngle(Dims const& o, geo::Vector_t const& v, const double OpDetOrientation) const;
    // circular aperture
    double Disk_SolidAngle(const double d, const double h, const double b) const;
    // dome aperture calculation
    double Omega_Dome_Model(const double distance, const double theta) const;

    // utility functions
    bool isOpDetInSameTPC(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) const;

    double fast_acos(double x) const;

    double interpolate(const std::vector<double>& xData,
                       const std::vector<double>& yData,
                       double x,
                       bool extrapolate,
                       size_t i = 0) const;

    double interpolate2(const std::vector<double>& xDistances,
                        const std::vector<double>& rDistances,
                        const std::vector<std::vector<std::vector<double>>>& parameters,
                        const double x,
                        const double r,
                        const size_t k) const;

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
    inline constexpr static bool
    isApproximatelyZero(TReal a, TReal tolerance = std::numeric_limits<TReal>::epsilon())
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

    // fhicl parameter sets
    const fhicl::ParameterSet fVUVHitsParams;
    const fhicl::ParameterSet fVISHitsParams;

    // geometry properties
    const larg4::ISTPC fISTPC;
    geo::GeometryCore const& fGeom;
    const int fNTPC;
    const std::vector<geo::BoxBoundedGeo> fActiveVolumes;
    const TVector3 fcathode_centre, fanode_centre;
    double fplane_depth, fanode_plane_depth;

    // photodetector geometry properties
    const size_t nOpDets;
    double fradius;
    Dims fcathode_plane;
    Dims fanode_plane;
    std::vector<geo::Point_t> fOpDetCenter;
    std::vector<int> fOpDetType;
    std::vector<int> fOpDetOrientation;
    std::vector<double> fOpDetLength;
    std::vector<double> fOpDetHeight;

    const int fvuv_absorption_length;

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

};

#endif