//////////////////////////////////////////////////////////////////////////////
/// \file ActiveVolumeVertexSampler.cxx
/// \brief Algorithm that samples vertex locations uniformly within the
/// active volume of a detector. It is fully experiment-agnostic and multi-TPC
/// aware.
///
/// \author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "ActiveVolumeVertexSampler.h"

namespace {
  constexpr int MAX_BOX_ITERATIONS = 10000;
}

//------------------------------------------------------------------------------
evgen::ActiveVolumeVertexSampler::ActiveVolumeVertexSampler(
  const fhicl::Table<evgen::ActiveVolumeVertexSampler::Config>& conf,
  rndm::NuRandomService& rand_service, const geo::Geometry& geom,
  const std::string& generator_name)
  : fVertexType(vertex_type_t::kSampled), fGeneratorName(generator_name),
  fTPCDist(nullptr)
{
  // Configure the algorithm using the FHiCL parameters
  this->reconfigure(conf, geom);

  // Register the TPC sampling engine with the seed service. If you need the
  // seed later, get it from the seed service using the value of the variable
  // generator_name as the instance name.
  rndm::NuRandomService::seed_t tpc_seed = rand_service.registerEngine(
    [this](rndm::NuRandomService::EngineId const& /* unused */,
      rndm::NuRandomService::seed_t lar_seed) -> void
    {
      auto seed = static_cast<uint_fast64_t>(lar_seed);
      // Use the obtained seed to prepare the random number engine.  This is
      // an attempt to do a decent job, but optimally accomplishing this can
      // be tricky (see, for example,
      // http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
      std::seed_seq seed_sequence{seed};
      fTPCEngine.seed(seed_sequence);
    },
    fGeneratorName, conf.get_PSet(), { "seed" }
  );

  // TODO: resolve the other workaround mentioned in the MARLEYHelper
  // class, then fix this as well
  uint_fast64_t tpc_cast_seed = static_cast<uint_fast64_t>(tpc_seed);
  std::seed_seq tpc_seed_sequence{tpc_cast_seed};
  fTPCEngine.seed(tpc_seed_sequence);
}

//------------------------------------------------------------------------------
TLorentzVector evgen::ActiveVolumeVertexSampler::sample_vertex_pos(
  const geo::Geometry& geom)
{
  // sample a new position if needed
  if (fVertexType == vertex_type_t::kSampled) {

    // Sample a TPC index using the active masses as weights
    size_t tpc_index = fTPCDist->operator()(fTPCEngine);

    // Get the dimensions of the chosen TPC's active volume
    const auto& tpc = geom.TPC(tpc_index);
    double minX = tpc.MinX();
    double maxX = tpc.MaxX();
    double minY = tpc.MinY();
    double maxY = tpc.MaxY();
    double minZ = tpc.MinZ();
    double maxZ = tpc.MaxZ();
    std::uniform_real_distribution<double>::param_type x_range(minX, maxX);
    std::uniform_real_distribution<double>::param_type y_range(minY, maxY);
    std::uniform_real_distribution<double>::param_type z_range(minZ, maxZ);

    // Sample a location uniformly over this volume
    std::uniform_real_distribution<double> uniform_dist;
    double x = uniform_dist(fTPCEngine, x_range);
    double y = uniform_dist(fTPCEngine, y_range);
    double z = uniform_dist(fTPCEngine, z_range);
    MF_LOG_INFO("ActiveVolumeVertexSampler " + fGeneratorName)
      << "Sampled primary vertex in TPC #" << tpc_index << ", x = " << x
      << ", y = " << y << ", z = " << z;

    // Update the vertex position 4-vector
    fVertexPosition.SetXYZT(x, y, z, 0.); // TODO: add time sampling
  }
  else if (fVertexType == vertex_type_t::kBox) {
    bool ok = false;
    int num_iterations = 0;

    // TODO: reduce code duplication here
    // Get the range of coordinates covered by the box
    std::uniform_real_distribution<double>::param_type x_range(fXmin, fXmax);
    std::uniform_real_distribution<double>::param_type y_range(fYmin, fYmax);
    std::uniform_real_distribution<double>::param_type z_range(fZmin, fZmax);

    // Sample a location uniformly over this volume
    std::uniform_real_distribution<double> uniform_dist;

    double x = 0.;
    double y = 0.;
    double z = 0.;
    while ( !ok && num_iterations < MAX_BOX_ITERATIONS ) {
      x = uniform_dist(fTPCEngine, x_range);
      y = uniform_dist(fTPCEngine, y_range);
      z = uniform_dist(fTPCEngine, z_range);

      // If we've been asked to check whether the vertex is within the
      // active volume of a TPC, verify that. Otherwise just accept the
      // first vertex position that was sampled unconditionally.
      if ( fCheckActive ) {
        size_t num_tpcs = geom.NTPC();
        for (size_t iTPC = 0; iTPC < num_tpcs; ++iTPC) {
          const auto& tpc = geom.TPC(iTPC);
          double minX = tpc.MinX();
          double maxX = tpc.MaxX();
          double minY = tpc.MinY();
          double maxY = tpc.MaxY();
          double minZ = tpc.MinZ();
          double maxZ = tpc.MaxZ();
          if ( x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ ) {
            ok = true;
            break;
          }
        }
      }
      else ok = true;
    }

    if ( !ok ) throw cet::exception("ActiveVolumeVertexSampler " + fGeneratorName)
      << "Failed to sample a vertex within a TPC active volume after " << MAX_BOX_ITERATIONS
      << " iterations";

    MF_LOG_INFO("ActiveVolumeVertexSampler " + fGeneratorName)
      << "Sampled primary vertex at x = " << x << ", y = " << y
      << ", z = " << z;

    // Update the vertex position 4-vector
    fVertexPosition.SetXYZT(x, y, z, 0.); // TODO: add time sampling
  }

  // if we're using a fixed vertex position, we don't need to do any sampling
  return fVertexPosition;
}

//------------------------------------------------------------------------------
void evgen::ActiveVolumeVertexSampler::reconfigure(
  const fhicl::Table<evgen::ActiveVolumeVertexSampler::Config>& conf,
  const geo::Geometry& geom)
{
  auto type = conf().type_();
  if (type == "sampled") {

    fVertexType = vertex_type_t::kSampled;

    // Get the active masses of each of the TPCs in the current geometry. Use
    // them as weights for sampling a TPC to use for the primary vertex.
    std::vector<double> tpc_masses;
    size_t num_tpcs = geom.NTPC();
    for (size_t iTPC = 0; iTPC < num_tpcs; ++iTPC) {
      // For each TPC, use its active mass (returned in kg) as its sampling
      // weight
      tpc_masses.push_back(geom.TPC(iTPC).ActiveMass());
    }

    // Load the discrete distribution used to sample TPCs with the up-to-date
    // set of weights
    fTPCDist.reset(new std::discrete_distribution<size_t>(tpc_masses.begin(),
      tpc_masses.end()));
  }
  else if (type == "fixed") {

    fVertexType = vertex_type_t::kFixed;

    auto vertex_pos = conf().position_();
    double Vx = vertex_pos.at(0);
    double Vy = vertex_pos.at(1);
    double Vz = vertex_pos.at(2);

    fVertexPosition.SetXYZT(Vx, Vy, Vz, 0.); // TODO: add time setting
  }
  else if (type == "box") {

    fVertexType = vertex_type_t::kBox;
    auto min_pos = conf().min_position_();
    auto max_pos = conf().max_position_();

    fXmin = min_pos.at(0);
    fYmin = min_pos.at(1);
    fZmin = min_pos.at(2);

    fXmax = max_pos.at(0);
    fYmax = max_pos.at(1);
    fZmax = max_pos.at(2);

    // By default, don't check whether each point is in a TPC active volume
    fCheckActive = false;
    // If the user specified this optional parameter, use that instead
    conf().check_active_( fCheckActive );
  }

  else throw cet::exception("ActiveVolumeVertexSampler " + fGeneratorName)
    << "Invalid vertex type '" << type << "' requested. Allowed values are"
    << " 'sampled' and 'fixed'";
}
