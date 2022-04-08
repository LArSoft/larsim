art_make(NO_PLUGINS
  EXCLUDE
    ScintTimeLAr_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    larcorealg_Geometry
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Services_Registry
    canvas::canvas
    cetlib::cetlib
    cetlib_except::cetlib_except
    fhiclcpp
    art::Utilities
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


simple_plugin(ScintTimeLAr "tool"
              larsim_PhotonPropagation_ScintTimeTools
              LIBRARIES PRIVATE CLHEP::CLHEP
              ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types


install_headers()
install_fhicl()
install_source()
