/**
 * @file   larsim/PhotonPropagation/PhotonMappingIdentityTransformations.h
 * @brief  A photon mapping identity transformation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 * @see    `larsim/PhotonPropagation/PhotonMappingIdentityTransformations.cxx`
 * 
 */

#ifndef LARSIM_PHOTONPROPAGATION_PHOTONMAPPINGIDENTITYTRANSFORMATIONS_H
#define LARSIM_PHOTONPROPAGATION_PHOTONMAPPINGIDENTITYTRANSFORMATIONS_H

// LArSoft libraries
#include "larsim/PhotonPropagation/IPhotonMappingTransformations.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"


namespace phot {
  
  /**
   * @brief Trivial photon mapping transformation.
   * 
   * This implementation of `phot::IPhotonMappingTransformation` interface
   * uses the identity mapping to have the library space matching the world one.
   * 
   * This is the most generic case when there are no symmetries to be exploited:
   * neither locations nor the optical detector identifiers are changed by the
   * mapping.
   * 
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
      , fDirectOpDetMap(makeDirectOpDetMap())
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
     * @see `locationFromLibrary()`
     * 
     * The returned vector is an exact copy of `location`.
     * 
     * No exception is ever thrown.
     */
    virtual geo::Point_t locationToLibrary
      (geo::Point_t const& location) const override
      { return location; }
    
    /// @}
    // --- END Geometry mapping interface --------------------------------------
    
    
    // --- BEGIN Optical detector identifier mapping interface -----------------
    /// @name Optical detector identifier mapping interface
    /// @{
    
    /**
     * @brief Maps an optical detector identifier with itself.
     * @param location world reference point _(unused)_
     * @param opDet optical detector ID
     *              (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @return identifier within the library of specified optical detector 
     * @see `opDetFromLibrary()`
     * 
     * The mapping is trivial: each optical detector has the same identifier
     * value in both the library and the detector.
     */
    virtual int opDetToLibrary
      (geo::Point_t const& location, OpDetID_t opDet) const override
      { return opDet; }
    
    /**
     * @brief Maps a library optical detector identifier with itself.
     * @param location world reference point _(unused)_
     * @param libOpDet optical detector ID internal to the library
     * @return identifier in the detector of specified optical detector 
     *         (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @see `opDetToLibrary()`
     * 
     * The mapping is trivial: each optical detector has the same identifier
     * value in both the library and the detector.
     */
    virtual int opDetFromLibrary
      (geo::Point_t const& location, OpDetID_t libOpDet) const override
      { return libOpDet; }
    
    
    /**
     * @brief Returns a map of identifiers within the library for each optical
     *        detector.
     * @param location world reference point _(unused)_
     * @return identifiers within the library for all optical detectors
     * @see `opDetToLibrary()`, `opDetsToLibrary()`
     * 
     * The mapping is trivial.
     * 
     * This is a static global mapping that does not depend on any location:
     * each optical detector maps to exactly one library optical detector.
     * 
     */
    virtual OpDetIDmap const& opDetsToLibrary
      (geo::Point_t const& location) const override
      { return fDirectOpDetMap; }
    
    /**
     * @brief Expected size of the mapping of optical detectors.
     * @return the expected size of the mapping of library optical detectors
     * @see `opDetsToLibrary()`
     * 
     * This is the number of available optical detectors.
     */
    virtual std::size_t opDetMappingSize() const override
      { return fDirectOpDetMap.size(); }
    
    
    /**
     * @brief Returns a map of identifiers of optical detectors for each optical
     *        detector identifier within the library around the specified
     *        `location`
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
    virtual LibOpDetIDmap opDetsFromLibrary
      (geo::Point_t const& location) const override
      { return fDirectOpDetMap; }
    
    /**
     * @brief Size of the mapping of library optical detectors.
     * @param location a world reference point to provide global context
     * @return the expected size of the mapping of library optical detectors
     * @see `opDetsFromLibrary()`
     */
    virtual std::size_t libraryOpDetMappingSize
      (geo::Point_t const& location) const override
      { return fDirectOpDetMap.size(); }
    
    
    /// @}
    // --- END Optical detector identifier mapping interface -------------------
    
    
      protected:
    
    /// Detector geometry service provider.
    geo::GeometryCore const* fGeom = nullptr;
    
    /// Complete optical detector identifier mapping world-to-library.
    OpDetIDmap fDirectOpDetMap;
    
    
    /// Return a global optical detector identifier mapping.
    virtual OpDetIDmap makeDirectOpDetMap() const;
    
  }; // class PhotonMappingIdentityTransformations
  
  
} // namespace phot



#endif // LARSIM_PHOTONPROPAGATION_PHOTONMAPPINGIDENTITYTRANSFORMATIONS_H
