cet_make_library(LIBRARY_NAME PhotonMappingTransformations INTERFACE
  SOURCE
  IPhotonMappingTransformations.h
  OpDetVisibilityData.h
  LIBRARIES INTERFACE
  larcorealg::headers
  lardataalg::UtilitiesHeaders
  larcoreobj::geo_vectors
)

cet_make_plugin_builder(lar::PhotonMappingTransformations art::tool Modules
  LIBRARIES CONDITIONAL
  larsim::PhotonMappingTransformations
)

include(lar::PhotonMappingTransformations)

cet_build_plugin(PhotonMappingIdentityTransformations
  lar::PhotonMappingTransformations
  IMPL_SOURCE PhotonMappingIdentityTransformations.cxx
  LIBRARIES
  CONDITIONAL
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcoreobj::geo_vectors
  PRIVATE
  larcorealg::Geometry
)

cet_build_plugin(PhotonMappingXMirrorTransformations
  lar::PhotonMappingTransformations
  IMPL_SOURCE PhotonMappingXMirrorTransformations.cxx
  LIBRARIES CONDITIONAL
  larsim::PhotonPropagation_LibraryMappingTools_PhotonMappingIdentityTransformations_tool
)

install_headers()
install_fhicl()
install_source()
