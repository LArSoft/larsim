art_make(
  MODULE_LIBRARIES
  #  lardata_Utilities
    nusimdata_SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    art_Persistency_Common canvas_Persistency_Common
    art_Persistency_Provenance canvas_Persistency_Provenance
    art_Utilities canvas_Utilities
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${CETLIB}
    ${ROOT_EG}
    ${ROOT_BASIC_LIB_LIST}
  )

install_headers()
install_fhicl()
install_source()
