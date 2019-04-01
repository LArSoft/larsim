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
   * It is common that the symmetries require also a rearrangement of the actual
   * optical detector channels: a hook for that type of mapping is also
   * provided, bidirectional.
   * 
   * This interface expresses the following transformations:
   * * `locationToLibrary()` transform a point in the detector (world) space
   *     to a point in the library space
   * * `opDetToLibrary()` transforms a optical detector channel ID into the
   *     corresponding ID as seen by the library
   * * `opDetFromLibrary()` transforms a detector ID in the library into a
   *     global optical detector ID
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
   * In that case, `util::MappedContainer` will not be sufficient any more,
   * but a similar class may be created.
   * 
   */
  class IPhotonMappingTransformations {
    
      public:
    
    // --- BEGIN Data types ----------------------------------------------------
    using OpDetID_t = int; ///< Type describing a optical detector ID.
    
    /**
     * @brief Type to return the mapping of optical IDs into library IDs.
     * 
     * This type guarantees:
     *  * `operator[](OpDetID_t)` to access the library ID of the specified
     *    optical detector ID
     *  * `size()` to report how many optical detectors are supported
     * 
     */
    using OpDetIDmap = std::vector<int>;
    
    /**
     * @brief Type to return the mapping of library IDs into optical IDs.
     * 
     * This type guarantees:
     *  * `operator[](OpDetID_t)` to access the optical detector ID of the
     *    specified library detector ID
     *  * `size()` to report how many library optical detectors are supported
     *  * either conversion to `bool`, response to `empty()` or to `std::size()`
     *    to determine if the contained data is valid
     *    (actually, just response to `phot::isValidLibraryData()`)
     * 
     */
    using LibOpDetIDmap = std::vector<int>;
    
    template <typename LibDataColl>
    using MappedOpDetData_t = phot::OpDetVisibilityData
      <
        util::collection_reference_t<LibDataColl>,
        util::collection_reference_t<LibOpDetIDmap const>
      >;
    
    // --- END Data types ------------------------------------------------------
    
    
    /// Value used for an invalid index (for both detector and library).
    static constexpr OpDetID_t InvalidIndex
      = std::numeric_limits<OpDetID_t>::max();
    
    
    // --- BEGIN Geometry mapping interface ------------------------------------
    /// @name Geometry mapping interface
    /// @{
    
    /**
     * @brief Returns the representation within the library of a detector
     *        location.
     * @param location position in world coordinates [cm]
     * @return a vector expressing `location` in the library space
     * @see `()`
     * @throw std::domain_error if the function can't handle `location`
     * 
     * The units and meaning of the returned location are 
     * 
     * Example: assuming that `fVoxelDef` contains the description of the
     * voxelization of the library, the ID of the voxel at `location` can be
     * obtained from a `mapping` as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const voxelID
     *   = fVoxelDef.GetVoxelID(mapping->locationToLibrary(location));
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where `mapping` points to an object implementing the
     * `phot::IPhotonMappingTransformations` interface.
     */
    virtual geo::Point_t locationToLibrary
      (geo::Point_t const& location) const = 0;
    
    /// @}
    // --- END Geometry mapping interface --------------------------------------
    
    
    // --- BEGIN Optical detector identifier mapping interface -----------------
    /// @name Optical detector identifier mapping interface
    /// @{
    
    /**
     * @brief Returns the identifier within the library of the specified
     *        optical detector.
     * @param location a world reference point to provide global context
     * @param opDet optical detector ID
     *              (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @return identifier within the library of specified optical detector 
     * @throw std::out_of_range if input optical detector can't be handled
     * @see `opDetFromLibrary()`
     * 
     * The global identifier `opDet` of an optical detector is converted into
     * an identifier internal to the library.
     */
    virtual int opDetToLibrary
      (geo::Point_t const& location, OpDetID_t opDet) const = 0;
    
    /**
     * @brief Returns the identifier in the detector of the optical detector
     *        specified via a library-specific identifier.
     * @param location a world reference point to provide global context
     * @param libOpDet optical detector ID internal to the library
     * @return identifier in the detector of specified optical detector 
     *         (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @throw std::out_of_range if input identifier can't be handled
     * @see `opDetToLibrary()`
     * 
     * The library-specific identifier `libOpDet` of an optical detector is
     * converted into an identifier valid for detector geometry.
     * 
     * The specified `location` is used to inform the conversion, and it is
     * usually pointing to the source of scintillation.
     * 
     * Note that such a reference is not needed for the inverse operation
     * (`opDetToLibrary()`) because in the mapping each optical detector must
     * map to exactly one library optical detector.
     */
    virtual int opDetFromLibrary
      (geo::Point_t const& location, OpDetID_t libOpDet) const = 0;
    
    
    /**
     * @brief Returns a map of identifiers within the library for each optical
     *        detector.
     * @param location a world reference point to provide global context
     * @return identifiers within the library for all optical detectors
     * @see `opDetToLibrary()`, `opDetsToLibrary()`
     * 
     * The returned value is a mapping object (see `OpDetIDmap` documentation
     * for the interface).
     * 
     * The specified `location` is used to inform the conversion, and it is
     * usually pointing to the source of scintillation. It may be left unused.
     * 
     * This is a static global mapping that does not depend on any location:
     * each optical detector maps to exactly one library optical detector.
     * 
     */
    virtual OpDetIDmap const& opDetsToLibrary
      (geo::Point_t const& location) const = 0;
    
    /**
     * @brief Expected size of the mapping of optical detectors.
     * @return the expected size of the mapping of library optical detectors
     * @see `opDetsToLibrary()`
     * 
     * This is effectively the number of available optical detectors, as well
     * as the size of the mapping as returned by `opDetsToLibrary()`.
     */
    virtual std::size_t opDetMappingSize() const = 0;
    
    /**
     * @brief Returns a map of identifiers of optical detectors for each optical
     *        detector identifier within the library around the specified
     *        `location`.
     * @param location a world reference point to provide global context
     * @return identifiers of optical detectors for all library identifiers
     * @throw std::domain_error if input identifier can't be handled
     * @see `opDetToLibrary()`, `opDetsFromLibrary()`
     * 
     * The returned value is a mapping object (see `LibOpDetIDmap` documentation
     * for the interface).
     * 
     * The specified `location` is used to provide context for the mapping, in
     * a similar fashion as `locationToLibrary()` does.
     * 
     * Differently from `opDetsFromLibrary()`, this method returns a copy of the
     * mapping rather than a reference to an existing one.
     * 
     */
    virtual LibOpDetIDmap const& opDetsFromLibrary
      (geo::Point_t const& location) const = 0;
    
    /**
     * @brief Expected size of the mapping of library optical detectors.
     * @param location a world reference point to provide global context
     * @return the expected size of the mapping of library optical detectors
     * @see `opDetsFromLibrary()`
     * 
     * This is effectively the size of the mapping as returned by
     * `opDetsFromLibrary()`.
     */
    virtual std::size_t libraryOpDetMappingSize
      (geo::Point_t const& location) const = 0;
    
    /// @}
    // --- END Optical detector identifier mapping interface -------------------
    
    
    // --- BEGIN Utility methods -----------------------------------------------
    /// @name Utility methods
    /// @{
    
    /**
     * @brief Remaps a collection indexed by library optical detector IDs into
     *        one indexed by (global) optical detector IDs.
     * @tparam Coll type of collection of data to be mapped
     * @param opDetMap the mapping to be applied
     * @tparam source collection of data to be mapped
     * @return a random access collection remapping the original data
     * 
     * This method returns a collection proxy that dynamically applies the
     * specified mapping. The original data is not copied.
     * 
     * The format of the map is an object with a indicing operator where the
     * index is the one in the destination collection (i.e. the optical detector
     * index) and the mapped value is the index in the source collection (the
     * libarry internal optical detector index). Each valid input index must be
     * mapped. If the mapped output index value is the largest
     * possible index (`InvalidIndex`), the item is not mapped anywhere and the
     * corresponding entry in the mapped collection is assigned the value
     * `defaultValue`. Otherwise, the index must be within the range of input
     * indices (which is up to `libraryOpDetMappingSize()` applied to the
     * location where `opDetMap` is valid).
     * 
     */
    template <typename Coll>
    MappedOpDetData_t<Coll> applyOpDetMapping(
      LibOpDetIDmap const& opDetMap, Coll&& source,
      util::collection_value_t<Coll> defaultValue = {}
      ) const;
    
    /**
     * @brief Remaps a collection indexed by library optical detector IDs into
     *        one indexed by (global) optical detector IDs.
     * @tparam Coll type of collection of data to be mapped
     * @param location a world reference point to provide global context
     * @param source collection of data to be mapped
     * @param defaultValue value returned for unmapped indices
     * @return a vector with copies of the original data, reshuffled
     * @see `applyOpDetMappingPtr()`
     * 
     * This method applies a mapping (see
     * `applyOpDetMapping(LibOpDetIDmap const&, Coll&&)`) that is
     * obtained via `opDetsFromLibrary()`.
     */
    template <typename Coll>
    auto applyOpDetMapping(
      geo::Point_t const& location, Coll&& source,
      util::collection_value_t<Coll> defaultValue = {}
      ) const
      {
        return applyOpDetMapping(
          opDetsFromLibrary(location), std::forward<Coll>(source), defaultValue
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
  LibOpDetIDmap const& opDetMap, Coll&& source,
  util::collection_value_t<Coll> defaultValue /* = {} */
) const 
  -> MappedOpDetData_t<Coll>
{
  
  auto const n = opDetMappingSize(); // number of available destination slots
  
  return MappedOpDetData_t<Coll>{
      util::make_collection_reference(std::forward<Coll>(source))
    , std::cref(opDetMap) // mapping is referenced
    , n                   // size
    , defaultValue
    };
  
} // phot::IPhotonMappingTransformations::applyOpDetMapping()


//------------------------------------------------------------------------------

#endif // LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_IPHOTONMAPPINGTRANSFORMATIONS_H
