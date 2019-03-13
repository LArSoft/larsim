//////////////////////////////////////////////////////////////////////////////
/// \file ActiveVolumeVertexSampler.h
/// \brief Algorithm that samples vertex locations uniformly within the
/// active volume of a detector. It is fully experiment-agnostic and multi-TPC
/// aware.
///
/// \author Steven Gardiner <sjgardiner@ucdavis.edu>
//////////////////////////////////////////////////////////////////////////////

#ifndef LARSIM_ALGORITHMS_ACTIVEVOLUMEVERTEXSAMPLER_H
#define LARSIM_ALGORITHMS_ACTIVEVOLUMEVERTEXSAMPLER_H

// standard library includes
#include <memory>
#include <random>
#include <string>

// framework includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"

// ROOT includes
#include "TLorentzVector.h"

namespace evgen {

  class ActiveVolumeVertexSampler {

    public:

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      /// Collection of configuration parameters used to
      /// determine the vertex location for each event
      struct Config {
        fhicl::Atom<std::string> type_ {
          Name("type"),
          Comment("Technique used to choose vertex locations"),
          "sampled" // default value
        };

        fhicl::OptionalAtom<std::string> seed_ {
          Name("seed"),
          Comment("Seed used for sampling vertex locations"),
          [this]() -> bool { return type_() != "fixed"; }
        };

        fhicl::Sequence<double, 3> position_ {
          Name("position"),
          Comment("Coordinates of the fixed vertex position"),
          [this]() -> bool { return type_() == "fixed"; }
        };

        fhicl::Sequence<double, 3> min_position_ {
          Name("min_position"),
          Comment("The minimum allowed values for the x, y, and z coordinates"),
          [this]() -> bool { return type_() == "box"; }
        };

        fhicl::Sequence<double, 3> max_position_ {
          Name("max_position"),
          Comment("The maximum allowed values for the x, y, and z coordinates"),
          [this]() -> bool { return type_() == "box"; }
        };

        fhicl::OptionalAtom<bool> check_active_ {
          Name("check_active"),
          Comment("Whether to enforce that the sampled vertices are within a TPC"
            " active volume"),
          [this]() -> bool { return type_() == "box"; }
        };

      }; // struct Config

      enum class vertex_type_t { kSampled, kFixed, kBox };

      // Configuration-checking constructors
      ActiveVolumeVertexSampler(const fhicl::Table<Config>& conf,
        rndm::NuRandomService& rand_service, const geo::Geometry& geom,
        const std::string& generator_name);

      ActiveVolumeVertexSampler(const fhicl::ParameterSet& pset,
        rndm::NuRandomService& rand_service, const geo::Geometry& geom,
        const std::string& generator_name)
        : ActiveVolumeVertexSampler(fhicl::Table<Config>(pset, {}),
        rand_service, geom, generator_name) {}

      void reconfigure(const fhicl::Table<Config>& conf,
        const geo::Geometry& geom);

      // Function that selects a primary vertex location for each event.
      // TODO: add time sampling
      TLorentzVector sample_vertex_pos(const geo::Geometry& geom);

    protected:

      // Currently sampled vertex position (doesn't change value if the vertex
      // is fixed)
      TLorentzVector fVertexPosition;

      vertex_type_t fVertexType;

      std::string fGeneratorName;

      // Discrete distribution object used to sample TPCs based on their active
      // masses
      std::unique_ptr<std::discrete_distribution<size_t> > fTPCDist;

      // RNG object used to sample TPCs
      std::mt19937_64 fTPCEngine;

      // Helper variables used only for "box" sampling mode
      double fXmin;
      double fYmin;
      double fZmin;

      double fXmax;
      double fYmax;
      double fZmax;

      bool fCheckActive;

  }; // class evgen::ActiveVolumeVertexSampler

} // namespace evgen

#endif
