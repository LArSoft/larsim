art_make(NO_PLUGINS
  EXCLUDE
    ScintTimeLAr_tool.cc
  LIB_LIBRARIES
    ${ROOT_CORE}
    larcorealg_Geometry
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${CANVAS}
    cetlib
    cetlib_except
    fhiclcpp
    ${CLHEP}
    ${ART_UTILITIES}
  )


simple_plugin(ScintTimeLAr "tool"
              larsim_PhotonPropagation_ScintTimeTools
              )


install_headers()
install_fhicl()
install_source()