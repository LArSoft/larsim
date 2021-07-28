art_make(NO_PLUGINS
  EXCLUDE
    ScintTimeLAr_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    larcorealg_Geometry
    art::Framework_Services_Registry
    ${CANVAS}
    cetlib::cetlib
    cetlib_except::cetlib_except
    fhiclcpp
    ${CLHEP}
    ${ART_UTILITIES}
  )


simple_plugin(ScintTimeLAr "tool"
              larsim_PhotonPropagation_ScintTimeTools
              LIBRARIES PRIVATE CLHEP::CLHEP
              )


install_headers()
install_fhicl()
install_source()