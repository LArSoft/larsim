# using `art_make()` instead of `art_make_library()` because too lazy to list all sources:
art_make(NO_PLUGINS
  EXCLUDE
    PhotonMappingIdentityTransformations_tool.cc
    PhotonMappingXMirrorTransformations_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    larcorealg_Geometry
    art::Framework_Services_Registry
    canvas::canvas
    art::Utilities
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(PhotonMappingIdentityTransformations "tool"
              larsim_PhotonPropagation_LibraryMappingTools) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(PhotonMappingXMirrorTransformations "tool"
              larsim_PhotonPropagation_LibraryMappingTools) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
