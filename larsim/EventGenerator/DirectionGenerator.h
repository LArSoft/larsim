/**
 * @author      : Daniele Guffanti (University and INFN Milano-Bicocca)
 * @file        : DirectionGenerator.h
 * @created     : Friday Oct 16, 2026 15:30:33 CEST
 * @description : Interface for algorithms that generate a particle direction.
 */

#ifndef LARSIM_DIRECTIONGENERATOR_H

#define LARSIM_DIRECTIONGENERATOR_H

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Vector_t
#include "CLHEP/Random/RandomEngine.h"

#include <sstream>
#include <type_traits>

namespace evgen {
  /// Interface for algorithms that generate a particle direction.
  class IDirectionGenerator {
    public:
      virtual ~IDirectionGenerator() noexcept = default;


      /// Dump the configuration of the generator to the output.
      template <typename Stream>
        void dumpSettings(Stream&& out) const 
        {
          if constexpr (std::is_base_of_v<std::ostream, std::remove_reference_t<Stream>>) {
            dumpSettingsImpl(out);
          } else {
            std::ostringstream oss;
            dumpSettingsImpl(oss);
            out << oss.str();
          }
        }

      /// Returns a unit direction vector in the detector frame (dimensionless).
      virtual geo::Vector_t generate(CLHEP::HepRandomEngine& engine) = 0;

    protected:
      virtual void dumpSettingsImpl(std::ostream& out) const = 0;
  };
} 

#endif /* end of include guard LARSIM_DIRECTIONGENERATOR_H */

