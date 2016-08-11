
add_subdirectory(Providers)

art_make(
  SERVICE_LIBRARIES
    art_Framework_Core
    art_Framework_Services_Registry
    art_Persistency_Provenance canvas_Persistency_Provenance
    art_Utilities canvas_Utilities
    art_Framework_Services_Optional_RandomNumberGenerator_service
    MF_MessageLogger
    MF_Utilities
    fhiclcpp
    cetlib
    ${ROOT_CORE}
)


install_source()
install_headers()
install_fhicl()
