art_make(NO_PLUGINS
  EXCLUDE
    ScintTimeLAr_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    larcorealg_Geometry
    art::Framework_Services_Registry
    ${CANVAS}
    cetlib
    cetlib_except
    fhiclcpp
    ${CLHEP_LIBRARIES}
    art::Utilities
  )


simple_plugin(ScintTimeLAr "tool"
              larsim_PhotonPropagation_ScintTimeTools
              )


install_headers()
install_fhicl()
install_source()