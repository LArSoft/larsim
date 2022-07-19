#ifndef PROPAGATIONTIMEMODEL_H
#define PROPAGATIONTIMEMODEL_H

// PropagationTimeModel
//  - parameterized fast optical simulation of photon propagation times
//  - contains functions to calculate the propagation times of direct and reflected photons incident
//  each photo-detector, along with the necessary ultility functions
//  - full description of model: Eur. Phys. J. C 81, 349 (2021)

// Nov 2021 by P. Green

// LArSoft libraries
#include "larsim/IonizationScintillation/ISTPC.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// fhicl
#include "fhiclcpp/ParameterSet.h"

namespace CLHEP {
  class HepRandomEngine;
  class RandFlat;
}

// other
#include "TVector3.h"
#include "TF1.h"

#include <array>
#include <memory>
#include <vector>


class PropagationTimeModel {

public:

  // constructor
  PropagationTimeModel(fhicl::ParameterSet VUVTimingParams,
                       fhicl::ParameterSet VISTimingParams,
                       CLHEP::HepRandomEngine& ScintTimeEngine,
                       bool doReflectedLight = false,
                       bool GeoPropTimeOnly = false);
  // destructor
  ~PropagationTimeModel();

  // propagation time
  void propagationTime(std::vector<double>& arrival_time_dist,
                         geo::Point_t const& x0,
                         const size_t OpChannel,
                         bool Reflected = false);

private:

  // parameter and geometry initialization
  void Initialization();

  // direct / VUV light
  void getVUVTimes(std::vector<double>& arrivalTimes,
                   const double distance_in_cm,
                   const size_t angle_bin);

  void getVUVTimesGeo(std::vector<double>& arrivalTimes,
                      const double distance_in_cm);

  void generateParam(const size_t index,
                     const size_t angle_bin);

  // reflected / visible light
  void getVISTimes(std::vector<double>& arrivalTimes,
                 const TVector3 &ScintPoint,
                 const TVector3 &OpDetPoint);

  // utility functions
  double fast_acos(double x) const;

  double interpolate(const std::vector<double>& xData,
                     const std::vector<double>& yData,
                     double x,
                     bool extrapolate,
                     size_t i = 0) const;

  void interpolate3(std::array<double, 3>& inter,
                    const std::vector<double>& xData,
                    const std::vector<double>& yData1,
                    const std::vector<double>& yData2,
                    const std::vector<double>& yData3,
                    double x,
                    bool extrapolate);

  static double finter_d(const double* x, const double* par);

  static double model_close(const double* x, const double* par);

  static double model_far(const double* x, const double* par);

  // fhicl parameter sets
  const fhicl::ParameterSet fVUVTimingParams;
  const fhicl::ParameterSet fVISTimingParams;

  // configuration
  const bool fdoReflectedLight;
  const bool fGeoPropTimeOnly;

  // ISTPC class
  larg4::ISTPC fISTPC;

  // random numbers
  CLHEP::HepRandomEngine& fScintTimeEngine;
  std::unique_ptr<CLHEP::RandFlat> fUniformGen;

  // geometry properties
  double fplane_depth;
  TVector3 fcathode_centre;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  // photodetector geometry properties
  size_t nOpDets;
  std::vector<geo::Point_t> fOpDetCenter;
  std::vector<int> fOpDetOrientation;

  // For VUV propagation time parametrization
  double fstep_size, fmax_d, fmin_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance, fangle_bin_timing_vuv;
  std::vector<std::vector<double>> fparameters[7];
  // vector containing generated VUV timing parameterisations
  std::vector<std::vector<TF1>> fVUV_timing;
  // vector containing min and max range VUV timing parameterisations are sampled to
  std::vector<std::vector<double>> fVUV_max;
  std::vector<std::vector<double>> fVUV_min;

  // For VIS propagation time parameterisation
  double fvis_vmean, fangle_bin_timing_vis;
  std::vector<double> fdistances_refl;
  std::vector<double> fradial_distances_refl;
  std::vector<std::vector<std::vector<double>>> fcut_off_pars;
  std::vector<std::vector<std::vector<double>>> ftau_pars;

};

#endif
