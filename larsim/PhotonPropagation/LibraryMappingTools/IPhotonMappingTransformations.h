/**
 * @file   larsim/PhotonPropagation/LibraryMappingTools/IPhotonMappingTransformations.h
 * @brief  Interface for transformation of photon visibility maps.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 *
 * This is a pure interface.
 */

#ifndef LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_IPHOTONMAPPINGTRANSFORMATIONS_H
#define LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_IPHOTONMAPPINGTRANSFORMATIONS_H

// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/OpDetVisibilityData.h"
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t<>
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// C++ standard libraries
#include <vector>
#include <utility> // std::forward()
#include <functional> // std::cref()
#include <stdexcept> // std::out_of_range, std::domain_error


namespace phot {

  /**
   * @brief Collection of functions to transform photon mapping data.
   *
   * If the optical detectors of an experiment sport some symmetry, it may be
   * convenient to exploit that symmetry to reduce the actual size of the
   * internal mapping to a basic volume, and using the symmetry rules to convert
   * between that basic volume and the world volume where the entire detector
   * is represented.
   * It is common that the symmetries require also a logic rearrangement of the
   * actual optical detectors: a hook for that type of mapping is also provided,
   * bidirectional.
   *
   * This interface expresses the following transformations:
   * * `detectorToLibrary()` transform a point in the detector space
   *     to a point in the library space
   * * `opDetToLibraryIndex()` transforms a optical detector ID into the
   *     corresponding library index
   * * `libraryIndexToOpDet()` transforms a library index into the corresponding
   *     global optical detector ID
   *
   *
   * Glossary
   * =========
   * @anchor PhotonMappingTransformations_Glossary
   *
   * In this class and its documentation the following terms are in use:
   *
   * * _library_ is a source of data; our libraries mapping the real space into
   *      some parametrization, and store the data as function of the
   *      parametrization values. In particular, they usually parametrize the
   *      covered space into tiny boxes ("voxels") and they address the data
   *      with a voxel index; and for each voxel index they store a collection
   *      of data values, addressed by another independent index that we call...
   * * _detector space_ or _world_ is the physical space where the whole
   *      detector is represented; in LArSoft this space is described in the
   *      global (or world) coordinate system, and the geometry description
   *      (provided by `geo::GeometryCore`) works in this space
   * * _library space_ is a representation of the physical space covered by the
   *      library; it may effectively be a subset of the world space (e.g.
   *      a single drift volume) or it may be a completely local representation,
   *      that needs a mapping with the world space in order to be of any use
   * * _library identifier_ or _library ID_ or _library index_ is the index used
   *      to access a single value of the library (after having chosen the voxel
   *      to get data for); note that the "library identifier" does not actually
   *      identify a library, but rather some specific data within the a
   *      library. This mapping class uses this concept to identify and select
   *      library data values, and uses the type `LibraryIndex_t` to represent
   *      it in the code
   * * _optical detector_ is a LArSoft concept, representing a physical detector
   *      with a presence in the real world; in LArSoft, they may be represented
   *      by...
   * * _optical detector number_ or _optical detector ID_ is a unique identifier
   *      in the form of either an integral number or as a location structure,
   *      `geo::OpDetID`; this concept is also used, e.g., as argument of
   *      `geo::GeometryCore::OpDetGeoFromOpDet()`. This mapping class uses the
   *      concept and the single unique number representation to identify the
   *      optical detector, and uses the local type `OpDetID_t` to represent it
   *      in the code;
   * * _optical detector readout channel_ or just _channel_ is a LArSoft
   * *    concept, representing a data single source from optical detectors, and
   *      associated with optical detector readout. In most detectors, there is
   *      one channel per optical detector and vice versa. It may be described
   *      in LArSoft by...
   * * _optical detector channel number_ or _channel number_ or (again)
   *   _channel_ is a unique numeric identifier for a optical detector readout
   *      channel; this is used, e.g., as argument of
   *      `geo::GeometryCore::OpDetGeoFromOpDet()`)
   *
   *
   * @note While optical detector and optical detector channels are different
   *       concepts, nothing this class is currently dealing with relies on any
   *       difference between the two. Nevertheless, given the nature of the
   *       library of providing visibility information based on the geometry of
   *       the detectors, using the concept of physical optical detector is more
   *       appropriate than using the one of readout channel.
   *
   *
   * Technical notes
   * ================
   *
   * Mapping
   * --------
   *
   * Currently the mapping is expressed with a specific type. If more
   * flexibility is needed, the mapping might use run time polymorphism.
   * This case is not supported right now. If need arose, something might be
   * designed on top of what is here now.
   *
   */
  class IPhotonMappingTransformations {

      public:

    virtual ~IPhotonMappingTransformations() = default;

    // --- BEGIN Data types ----------------------------------------------------
    using OpDetID_t = int; ///< Type describing a optical detector ID.

    /// Type describing a library index. FIXME former LibraryOpDetID_t
    using LibraryIndex_t = OpDetID_t; // currently bound to channel ID

    /**
     * @brief Type describing the mapping of optical detectors into library
     *        indices.
     *
     * This type of mapping allows querying data values by optical detector ID,
     * fetching them from library data source.
     *
     * This type guarantees:
     *  * `operator[](OpDetID_t)` to fetch the library index associated to
     *    the specified optical detector
     *  * `size()` to report how many optical detectors are supported
     *
     */
    using LibraryIndexToOpDetMap = std::vector<OpDetID_t>;

    /**
     * @brief Type describing the mapping of library indices into optical
     *        detectors.
     *
     * This type of mapping allows querying data values by library index,
     * fetching them from an optical-detector-based data source.
     *
     * This type guarantees:
     *  * `operator[](LibraryIndex_t)` to access the optical detector ID of
     *    the specified library index
     *  * `size()` to report how many library optical detectors are supported
     *  * either conversion to `bool`, response to `empty()` or to `std::size()`
     *    to determine if the contained data is valid
     *    (actually, just response to `phot::isValidLibraryData()`)
     *
     */
    using OpDetToLibraryIndexMap = std::vector<LibraryIndex_t>; // FIXME former LibOpDetIDmap

    /**
     * @brief Type of optical detector data collection.
     * @tparam LibDataColl type of the original library data container
     *
     * This type behaves like a random access collection with as index the
     * optical detector ID and as value the data corresponding to that optical
     * detector.
     */
    template <typename LibDataColl>
    using MappedOpDetData_t = phot::OpDetVisibilityData
      <
        util::collection_reference_t<LibDataColl>,
        util::collection_reference_t<OpDetToLibraryIndexMap const>
      >;

    // --- END Data types ------------------------------------------------------


    /**
     * @brief Value used to identify an invalid optical detector.
     *
     * It may be used as value of the library-to-detector mapping.
     */
    static constexpr OpDetID_t InvalidOpDetID
      = std::numeric_limits<OpDetID_t>::max();

    /**
     * @brief Value used for an invalid library index.
     *
     * It may be used as value of the detector-to-library mapping.
     */
    static constexpr LibraryIndex_t InvalidLibraryIndex
      = std::numeric_limits<LibraryIndex_t>::max();


    // --- BEGIN Geometry mapping interface ------------------------------------
    /// @name Geometry mapping interface
    /// @{

    /**
     * @brief Returns the representation within the library of a detector
     *        location.
     * @param location position in world coordinates [cm]
     * @return a vector expressing `location` in the library space
     * @throw std::domain_error if the function can't handle `location`
     *
     * The units and meaning of the returned location are library-dependent.
     *
     * Example: assuming that `fVoxelDef` contains the description of the
     * voxelization of the library, the ID of the voxel at `location` can be
     * obtained from a `mapping` as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const voxelID
     *   = fVoxelDef.GetVoxelID(mapping->detectorToLibrary(location));
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where `mapping` points to an object implementing the
     * `phot::IPhotonMappingTransformations` interface.
     * Note that the caller here has no business with the translated location,
     * beside propagating it to the library: caller simply is not explained
     * the meaning of that translated location.
     */
    virtual geo::Point_t detectorToLibrary
      (geo::Point_t const& location) const = 0;

    /// @}
    // --- END Geometry mapping interface --------------------------------------


    // --- BEGIN Optical detector mapping interface ----------------------------
    /**
     * @name Optical detector mapping interface
     *
     * Here the mapping is provided between optical detectors (`OpDetID_t`) and
     * library indices (`LibraryIndex_t`).
     *
     * Terms and roles of the various concepts involved in the mapping may be
     * confusing: see the @ref PhotonMappingTransformations_Glossary "glossary"
     * in `phot::IPhotonMappingTransformations` class description for their
     * definition.
     */
    /// @{

    /**
     * @brief Returns the library index for the specified optical detector.
     * @param location a point in detector space to provide global context
     * @param opDetID optical detector identifier
     * @return index corresponding to `opDetID`, or `InvalidLibraryIndex`
     * @throw std::out_of_range if input optical detector can't be handled
     * @throw std::domain_error if location can't be handled
     * @see `libraryIndexToOpDet()`, `opDetsToLibraryIndices()`
     *
     * The specified global identifier `opDetID` of an optical detector is
     * converted into an index for the library.
     *
     * If the library does not include information for this optical detector,
     * the value `InvalidLibraryIndex` is returned instead.
     *
     * The specified `location` is used to inform the conversion, and it is
     * usually pointing to the source of scintillation.
     */
    virtual LibraryIndex_t opDetToLibraryIndex
      (geo::Point_t const& location, OpDetID_t opDetID) const = 0;

    /**
     * @brief Returns the optical detector ID for the specified library index.
     * @param location a point in detector space to provide global context
     * @param libIndex library index to be mapped
     * @return optical detector corresponding to `libIndex`, or `InvalidOpDetID`
     * @throw std::out_of_range if input index can't be handled
     * @throw std::domain_error if location can't be handled
     * @see `opDetToLibraryIndex()`, `libraryIndicesToOpDets()`
     *
     * The library index `libIndex` is converted into a optical detector ID.
     * If the library does not include information for this optical detector,
     * the value `InvalidLibraryIndex` is returned instead.
     *
     * The specified `location` is used to inform the conversion, and it is
     * usually pointing to the source of scintillation.
     */
    virtual OpDetID_t libraryIndexToOpDet
      (geo::Point_t const& location, LibraryIndex_t libIndex) const = 0;


    /**
     * @brief Returns a map of library indices as function of optical detectors.
     * @param location a world reference point to provide global context
     * @return library indices for all optical detectors
     * @throw std::domain_error if location can't be handled
     * @see `opDetToLibraryIndex()`, `opDetsToLibraryIndices()`
     *
     * The returned value is a mapping object (see `OpDetToLibraryIndexMap`
     * documentation for the interface). If an optical detector has no
     * representation in the library at `location`, the index value
     * corresponding to it is `InvalidLibraryIndex`.
     *
     * The specified `location` is used to inform the conversion, and it is
     * usually pointing to the source of scintillation. It may be left unused.
     */
    virtual OpDetToLibraryIndexMap const& opDetsToLibraryIndices
      (geo::Point_t const& location) const = 0;

    /**
     * @brief Expected number of mappings of optical detector into library
     *        index.
     * @return the expected size of the mapping of optical detectors
     * @see `opDetsToLibraryIndices()`
     *
     * This is effectively the number of available optical detectors, as well
     * as the size of the mapping as returned by `opDetsToLibraryIndices()`.
     */
    virtual std::size_t opDetMappingSize() const = 0;

    /**
     * @brief Returns a map of optical detectors identifiers, one for each
     *        library index.
     * @param location a world reference point to provide global context
     * @return optical detector identifiers for all library indices
     * @throw std::domain_error if location can't be handled
     * @see `opDetsToLibraryIndices()`, `libraryIndexToOpDet()`
     *
     * The returned value is a mapping object (see `LibraryIndexToOpDetMap`
     * documentation for the interface). If a library index does not map to any
     * optical detector in the library at `location` (which is unusual!),
     * the optical detector ID corresponding to it is `InvalidOpDetID`.
     *
     * The specified `location` is used to provide context in a similar
     * fashion as `detectorToLibrary()` does. It can be used to choose the
     * correct mapping among the available ones.
     */
    virtual LibraryIndexToOpDetMap const& libraryIndicesToOpDets
      (geo::Point_t const& location) const = 0;

    /**
     * @brief Expected size of the mapping from library to optical detectors.
     * @param location a world reference point to provide global context
     * @throw std::domain_error if location can't be handled
     * @return the expected size of the mapping from library indices to
     *         optical detectors
     * @see `libraryIndicesToOpDets()`
     *
     * This is effectively the size of the mapping returned by
     * `libraryIndicesToOpDets()`. It represents how many library indices are
     * provided by the library for the specified `location`. For example, in
     * an implementation exploiting the modularity of a detector with _N_
     * channels and made of _M_ modules, this value might be _N / M_.
     */
    virtual std::size_t libraryMappingSize
      (geo::Point_t const& location) const = 0;

    /// @}
    // --- END Optical detector mapping interface ------------------------------


    // --- BEGIN Utility methods -----------------------------------------------
    /// @name Utility methods
    /// @{

    /**
     * @brief Remaps a collection indexed by library index into one indexed by
     *        optical detector IDs.
     * @tparam Coll type of collection of data from the library, to be mapped
     * @param opDetToLibraryMap the mapping to be applied
     * @param source collection of library data to be mapped
     * @param defaultValue value returned for unmapped optical detector IDs
     * @return a random access collection providing data by optical detector ID
     *
     * This method returns a collection proxy that dynamically applies the
     * specified mapping. The original data is not copied.
     *
     * The format of the map (`opDetMap`) is an object with a indexing operator
     * where the index is _the one in the destination collection_ (i.e. the
     * optical detector ID) and the mapped value is the index in the source
     * collection (the libary index). Each valid optical detector ID must be
     * mapped. If the library index mapped to an optical detector is
     * `InvalidLibraryIndex`, that optical detector is not mapped anywhere and
     * its entry in the mapped collection is assigned the value `defaultValue`.
     * Otherwise, the library index must be valid (which means, up to the result
     * of `libraryMappingSize()` applied to the location where `opDetMap` is
     * valid).
     *
     */
    template <typename Coll>
    MappedOpDetData_t<Coll> applyOpDetMapping(
      OpDetToLibraryIndexMap const& opDetToLibraryMap, Coll&& source,
      util::collection_value_t<Coll> defaultValue = {}
      ) const;

    /**
     * @brief Remaps a collection indexed by library index into one indexed by
     *        optical detector IDs according to the mapping at `location`.
     * @tparam Coll type of collection of data from the library, to be mapped
     * @param location a world reference point to provide global context
     * @param source collection of library data to be mapped
     * @param defaultValue value returned for unmapped optical detector IDs
     * @return a random access collection providing data by optical detector ID
     * @see `applyOpDetMapping()`
     *
     * This method applies a mapping (see
     * `applyOpDetMapping(OpDetToLibraryIndexMap const&, Coll&&, util::collection_value_t<Coll>)`)
     * that is obtained via `libraryIndicesToOpDets()`.
     */
    template <typename Coll>
    auto applyOpDetMapping(
      geo::Point_t const& location, Coll&& source,
      util::collection_value_t<Coll> defaultValue = {}
      ) const
      {
        return applyOpDetMapping(
          opDetsToLibraryIndices(location), std::forward<Coll>(source),
          defaultValue
          );
      }

    /// @}
    // --- END Utility methods -------------------------------------------------


  }; // class IPhotonMappingTransformations


} // namespace phot



//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
template <typename Coll>
auto phot::IPhotonMappingTransformations::applyOpDetMapping(
  OpDetToLibraryIndexMap const& opDetToLibraryMap, Coll&& source,
  util::collection_value_t<Coll> defaultValue /* = {} */
) const
  -> MappedOpDetData_t<Coll>
{

  using std::size;

  // number of available destination slots (i.e. how many optical detectors)
  auto const n = size(opDetToLibraryMap);

  return MappedOpDetData_t<Coll>{
      util::make_collection_reference(std::forward<Coll>(source))
    , std::cref(opDetToLibraryMap) // mapping is referenced
    , n                            // size
    , defaultValue
    };

} // phot::IPhotonMappingTransformations::applyOpDetMapping()


//------------------------------------------------------------------------------

#endif // LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_IPHOTONMAPPINGTRANSFORMATIONS_H
