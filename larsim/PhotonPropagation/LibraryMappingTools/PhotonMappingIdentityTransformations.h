/**
 * @file   larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.h
 * @brief  A photon mapping identity transformation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 * @see    `larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations_tool.cc`
 *
 */

#ifndef LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_PHOTONMAPPINGIDENTITYTRANSFORMATIONS_H
#define LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_PHOTONMAPPINGIDENTITYTRANSFORMATIONS_H

// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/IPhotonMappingTransformations.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

// framework libraries
#include "art/Utilities/ToolConfigTable.h"

namespace geo { class GeometryCore; }

namespace phot {

  /**
   * @brief Trivial photon mapping transformation.
   *
   * This implementation of `phot::IPhotonMappingTransformations` interface
   * uses the identity mapping to have the library space matching the world one.
   *
   * This is the most generic case when there are no symmetries to be exploited:
   * neither locations nor the optical detectors are changed by the mapping.
   *
   * @note The documentation of this class describes the details of this
   *       implementation of `phot::IPhotonMappingTransformations` interface.
   *       For information about how to _use it_, see the documentation of that
   *       interface.
   */
  class PhotonMappingIdentityTransformations
    : public phot::IPhotonMappingTransformations
  {

      public:
    struct Config {

    //   using Name = fhicl::Name;
    //   using Comment = fhicl::Comment;

      // no configuration required

    }; // struct Config

    using Parameters = art::ToolConfigTable<Config>;


    /// Constructor: ignores the configuration.
    PhotonMappingIdentityTransformations(Config const&)
      : fGeom(lar::providerFrom<geo::Geometry>())
      , fOpDetsToLibraryIndicesMap(makeOpDetsToLibraryIndicesMap())
      , fLibraryIndicesToOpDetsMap(makeLibraryIndicesToOpDetsMap())
      {}

    /// Constructor: ignores the configuration.
    PhotonMappingIdentityTransformations(Parameters const& config)
      : PhotonMappingIdentityTransformations(config())
      {}


    // --- BEGIN Geometry mapping interface ------------------------------------
    /// @name Geometry mapping interface
    /// @{

    /**
     * @brief Returns the representation within the library of a detector
     *        location.
     * @param location position in world coordinates [cm]
     * @return a vector expressing `location` in the library space
     *
     * The returned vector is an exact copy of `location`.
     *
     * No exception is ever thrown.
     */
    virtual geo::Point_t detectorToLibrary
      (geo::Point_t const& location) const override
      { return location; }

    /// @}
    // --- END Geometry mapping interface --------------------------------------


    // --- BEGIN Optical detector mapping interface ----------------------------
    /// @name Optical detector mapping interface
    /// @{

    /**
     * @brief Maps an optical detector with its library index.
     * @param location _(unused)_ world reference point
     * @param opDetID ID of the optical detector to be mapped
     *              (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @return library index corresponding to the specified optical detector
     * @see `libraryIndexToOpDet()`, `opDetsToLibraryIndices()`
     *
     * The mapping is trivial: each library index has the same value as the
     * ID of the optical detector it corresponds to.
     */
    virtual LibraryIndex_t opDetToLibraryIndex
      (geo::Point_t const& location, OpDetID_t opDetID) const override
      { return LibraryIndex_t{ opDetID }; }

    /**
     * @brief Maps a library index with an optical detector.
     * @param location world reference point _(unused)_
     * @param libIndex the library index to be mapped
     * @return optical detector ID
     *         (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @see `opDetToLibraryIndex()`, `libraryIndicesToOpDets()`
     *
     * The mapping is trivial: each library index has the same value as the
     * ID of the optical detector it corresponds to.
     */
    virtual OpDetID_t libraryIndexToOpDet
      (geo::Point_t const& location, LibraryIndex_t libIndex) const override
      { return OpDetID_t{ libIndex }; }


    /**
     * @brief Returns a map of library indices as function of optical detectors.
     * @param location world reference point _(unused)_
     * @return library indices for all optical detectors
     * @see `opDetToLibraryIndex()`, `opDetsToLibraryIndices()`
     *
     * The mapping is trivial, as each library index has the same value as the
     * ID of the optical detector it corresponds to.
     * This mapping is global and does not depend on any location.
     */
    virtual OpDetToLibraryIndexMap const& opDetsToLibraryIndices
      (geo::Point_t const& location) const override
      { return fOpDetsToLibraryIndicesMap; }

    /**
     * @brief Expected number of mappings of optical detector into library
     *        index.
     * @return the expected size of the mapping of optical detectors
     * @see `opDetsToLibraryIndices()`
     *
     * This is effectively the number of available optical detectors, as well
     * as the size of the mapping as returned by `opDetsToLibraryIndices()`.
     */
    virtual std::size_t opDetMappingSize() const override
      { return fOpDetsToLibraryIndicesMap.size(); }


    /**
     * @brief Returns a map of optical detectors identifiers, one for each
     *        library index.
     * @param location a world reference point _(unused)_
     * @return optical detector identifiers for all library indices
     * @see `opDetsToLibraryIndices()`, `libraryIndexToOpDet()`
     *
     * The mapping is trivial, as each library index has the same value as the
     * ID of the optical detector it corresponds to.
     * This mapping is global and does not depend on any location.
     */
    virtual LibraryIndexToOpDetMap const& libraryIndicesToOpDets
      (geo::Point_t const& location) const override
      { return fLibraryIndicesToOpDetsMap; }

    /**
     * @brief Size of the mapping of library optical detectors.
     * @param location a world reference point _(unused)_
     * @return the expected size of the mapping of library optical detectors
     * @see `opDetsFromLibrary()`
     *
     * This is also the number of optical detectors.
     */
    virtual std::size_t libraryMappingSize
      (geo::Point_t const& location) const override
      { return fLibraryIndicesToOpDetsMap.size(); }


    /// @}
    // --- END Optical detector identifier mapping interface -------------------


      protected:

    /// Detector geometry service provider.
    geo::GeometryCore const* fGeom = nullptr;

    /// Complete optical detector identifier mapping world-to-library.
    OpDetToLibraryIndexMap fOpDetsToLibraryIndicesMap;

    /// Complete optical detector identifier mapping library-to-world.
    LibraryIndexToOpDetMap fLibraryIndicesToOpDetsMap;


    /// Return a trivial global optical detector identifier mapping.
    virtual LibraryIndexToOpDetMap makeLibraryIndicesToOpDetsMap() const;

    /// Return a trivial library index mapping.
    virtual OpDetToLibraryIndexMap makeOpDetsToLibraryIndicesMap() const;

  }; // class PhotonMappingIdentityTransformations


} // namespace phot


#endif // LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_PHOTONMAPPINGIDENTITYTRANSFORMATIONS_H
