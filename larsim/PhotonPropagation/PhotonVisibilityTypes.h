/**
 * @file   larsim/PhotonPropagation/PhotonVisibilityTypes.h
 * @brief  Declaration of types related to photon visibility.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 5, 2019
 * @see    larsim/PhotonPropagation/PhotonVisibilityService.h
 * 
 * The purpose of this header is to allow the inclusion of the types needed
 * for interacting with `phot::PhotonVisibilityService` without having to
 * include the full service.
 * 
 * This is a header-only library.
 */

#ifndef LARSIM_PHOTONPROPAGATION_PHOTONVISIBILITYTYPES_H
#define LARSIM_PHOTONPROPAGATION_PHOTONVISIBILITYTYPES_H


// LArSoft libraries
#include "larsim/PhotonPropagation/IPhotonLibrary.h"
#include "larsim/PhotonPropagation/LibraryMappingTools/IPhotonMappingTransformations.h"


namespace phot{
  
  /// Type of (global) optical detector ID.
  using OpDetID_t = phot::IPhotonMappingTransformations::OpDetID_t;
  
  /**
   * @brief Type of mapped visibility counts.
   *
   * No data storage is provided.
   * 
   * This is the type returned by `phot::PhotonVisibilityService` when asked
   * about visibility of a point from _all_ the optical detectors.
   * 
   * It used to be `float const*`.
   */
  using MappedCounts_t
    = phot::IPhotonMappingTransformations::MappedOpDetData_t
      <phot::IPhotonLibrary::Counts_t>
    ;
  
  /**
   * @brief Type of mapped arrival times.
   *
   * No data storage is provided.
   * 
   * This is the type returned by `phot::PhotonVisibilityService` when asked
   * about time of arrival from a point to _all_ the optical detectors.
   */
  using MappedT0s_t
    = phot::IPhotonMappingTransformations::MappedOpDetData_t
      <phot::IPhotonLibrary::T0s_t>
    ;
  
  /**
   * @brief Type of set of parameters for functions.
   *
   * No data storage is provided.
   * 
   * This is the type returned by `phot::PhotonVisibilityService` when asked
   * about any parametrized... parameter, from a point to _all_ the optical
   * detectors.
   */
  using MappedParams_t
    = phot::IPhotonMappingTransformations::MappedOpDetData_t
      <phot::IPhotonLibrary::Params_t>
    ;
  
  /**
   * @brief Type of mapped parametrization functions.
   *
   * No data storage is provided.
   * 
   * This is the type returned by `phot::PhotonVisibilityService` when asked
   * about any parameter directlyt described by functions (glorious ROOT `TF1`),
   * from a point to _all_ the optical detectors.
   */
  using MappedFunctions_t
    = phot::IPhotonMappingTransformations::MappedOpDetData_t
      <phot::IPhotonLibrary::Functions_t>
    ;
  
} // namespace phot


#endif // LARSIM_PHOTONPROPAGATION_PHOTONVISIBILITYTYPES_H
