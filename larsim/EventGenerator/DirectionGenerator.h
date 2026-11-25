/**
 * @author      : Daniele Guffanti (University and INFN Milano-Bicocca)
 * @file        : DirectionGenerator.h
 * @created     : Friday Oct 16, 2026 15:30:33 CEST
 * @brief       : Interface for particle direction generators.
 */

#ifndef LARSIM_DIRECTIONGENERATOR_H

#define LARSIM_DIRECTIONGENERATOR_H

#include "CLHEP/Random/RandomEngine.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Vector_t

#include <sstream>
#include <type_traits>

namespace evgen {
  /**
   * @brief Interface for algorithms that generate particle directions.
   *
   * Implementations sample unit vectors in the detector coordinate frame
   * using the random engine supplied by the caller.
   */
  class IDirectionGenerator {
  public:
    virtual ~IDirectionGenerator() noexcept = default;

    /**
       * @brief Writes the generator configuration to an output stream.
       * @param out Output stream or stream-like object accepting a string.
       */
    template <typename Stream>
    void dumpSettings(Stream&& out) const
    {
      if constexpr (std::is_base_of_v<std::ostream, std::remove_reference_t<Stream>>) {
        dumpSettingsImpl(out);
      }
      else {
        std::ostringstream oss;
        dumpSettingsImpl(oss);
        out << oss.str();
      }
    }

    /// Returns a unit direction vector in the detector frame (dimensionless).
    virtual geo::Vector_t generate(CLHEP::HepRandomEngine& engine) = 0;

  protected:
    /**
       * @brief Writes implementation-specific settings to a standard stream.
       * @param out Output stream.
       */
    virtual void dumpSettingsImpl(std::ostream& out) const = 0;
  };
}

#endif /* end of include guard LARSIM_DIRECTIONGENERATOR_H */
