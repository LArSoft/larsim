art_make(
  MODULE_LIBRARIES
  #  Utilities
    ${SIMULATIONBASE}
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_PERSISTENCY_COMMON}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${CETLIB}
    ${ROOT_BASIC_LIB_LIST}
  )

install_headers()
install_fhicl()
install_source()
