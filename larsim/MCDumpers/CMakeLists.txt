art_make(
  LIB_LIBRARIES
    ${ROOT_EG}
    ${ROOT_CORE}
  MODULE_LIBRARIES
  #  lardata_Utilities
    nusimdata_SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    cetlib cetlib_except
    ${ROOT_EG}
    ${ROOT_BASIC_LIB_LIST}
  )

install_headers()
install_fhicl()
install_source()
