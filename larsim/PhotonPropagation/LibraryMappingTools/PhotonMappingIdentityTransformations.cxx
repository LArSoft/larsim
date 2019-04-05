/**
 * @file   larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.cxx
 * @brief  A photon mapping identity transformation: implementation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 * @see    `larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.h`
 * @see    `larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations_tool.cc`
 * 
 */


// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.h"


//------------------------------------------------------------------------------
phot::IPhotonMappingTransformations::LibraryIndexToOpDetMap
phot::PhotonMappingIdentityTransformations::makeLibraryIndicesToOpDetsMap
  () const
{
  
  auto const nOpDets = fGeom->NOpDets();
  
  LibraryIndexToOpDetMap map(nOpDets);
  std::iota(map.begin(), map.end(), 0);
  
  return map;
} // phot::PhotonMappingIdentityTransformations::makeLibraryIndicesToOpDetsMap()


//------------------------------------------------------------------------------
phot::IPhotonMappingTransformations::OpDetToLibraryIndexMap
phot::PhotonMappingIdentityTransformations::makeOpDetsToLibraryIndicesMap
  () const
{
  
  auto const nOpDets = fGeom->NOpDets();
  
  OpDetToLibraryIndexMap map(nOpDets);
  std::iota(map.begin(), map.end(), 0);
  
  return map;
} // phot::PhotonMappingIdentityTransformations::makeOpDetsToLibraryIndicesMap()


//------------------------------------------------------------------------------

