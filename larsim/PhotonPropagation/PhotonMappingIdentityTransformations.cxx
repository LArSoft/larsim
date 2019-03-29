/**
 * @file   larsim/PhotonPropagation/PhotonMappingIdentityTransformations.cxx
 * @brief  A photon mapping identity transformation: implementation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 * @see    `larsim/PhotonPropagation/PhotonMappingIdentityTransformations.h`
 * 
 */


// LArSoft libraries
#include "larsim/PhotonPropagation/PhotonMappingIdentityTransformations.h"


//------------------------------------------------------------------------------
phot::IPhotonMappingTransformations::OpDetIDmap
phot::PhotonMappingIdentityTransformations::makeDirectOpDetMap() const {
  
  auto const nOpDets = fGeom->NOpDets();
  
  OpDetIDmap map(nOpDets);
  std::iota(map.begin(), map.end(), 0);
  
  return map;
} // phot::PhotonMappingIdentityTransformations::makeDirectOpDetMap()


//------------------------------------------------------------------------------

