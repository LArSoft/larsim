# using `art_make()` instead of `art_make_library()` because too lazy to list all sources:
art_make(NO_PLUGINS
  EXCLUDE
    PhotonMappingIdentityTransformations_tool.cc
    PhotonMappingXMirrorTransformations_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    larcorealg_Geometry
    art::Framework_Services_Registry
    ${CANVAS}
    ${ART_UTILITIES}
  )

simple_plugin(PhotonMappingIdentityTransformations "tool"
              larsim_PhotonPropagation_LibraryMappingTools)

simple_plugin(PhotonMappingXMirrorTransformations "tool"
              larsim_PhotonPropagation_LibraryMappingTools)

install_headers()
install_fhicl()
install_source()
