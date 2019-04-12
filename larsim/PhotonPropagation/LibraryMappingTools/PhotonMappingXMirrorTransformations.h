/**
 * @file   larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingXMirrorTransformations.h
 * @brief  A photon mapping transformation with reflection at x = 0.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 * @see    `larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingXMirrorTransformations_tool.cc`
 * 
 */

#ifndef LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_PHOTONMAPPINGXMIRRORTRANSFORMATIONS_H
#define LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_PHOTONMAPPINGXMIRRORTRANSFORMATIONS_H

// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"

// C/C++ standard libraries
#include <cmath> // std::abs()

namespace phot {
  
  /**
   * @brief Photon mapping transformation with reflection at x = 0.
   * 
   * This mapping describes an infinite planar mirror placed at @f$ x = 0 @f$.
   * 
   */
  class PhotonMappingXMirrorTransformations
    : public phot::PhotonMappingIdentityTransformations
  {
    
      public:
    struct Config: phot::PhotonMappingIdentityTransformations::Config {
      
    //  using Name = fhicl::Name;
    //  using Comment = fhicl::Comment;
      
      // no additional configuration required
      
    }; // struct Config
    
    using Parameters = art::ToolConfigTable<Config>;
    
    
    /// Constructor: same configuration
    /// as `phot::PhotonMappingIdentityTransformations`.
    PhotonMappingXMirrorTransformations(Parameters const& config)
      : PhotonMappingXMirrorTransformations(config()) {}
    
    /// Constructor: same configuration
    /// as `phot::PhotonMappingIdentityTransformations`.
    PhotonMappingXMirrorTransformations(Config const& config)
      : PhotonMappingIdentityTransformations(config) {}
    
    
    /**
     * @brief Returns the representation within the library of a detector
     *        location.
     * @param location position in world coordinates [cm]
     * @return a vector expressing `location` in the library space
     * 
     * The returned vector is the same as `location`, but with the _x_ component
     * always positive..
     * 
     * No exception is ever thrown.
     */
    virtual geo::Point_t detectorToLibrary
      (geo::Point_t const& location) const override
      { return { std::abs(location.X()), location.Y(), location.Z() }; }
    
  }; // class PhotonMappingXMirrorTransformations
  
  
} // namespace phot


#endif // LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_PHOTONMAPPINGXMIRRORTRANSFORMATIONS_H
