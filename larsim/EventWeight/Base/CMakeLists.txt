art_make(LIB_LIBRARIES
         nutools_RandomUtils_NuRandomService_service
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         canvas
         cetlib_except
         ${FHICLCPP}
         ${CLHEP}
         ${ROOT_CORE}
         ${ROOT_MATRIX}
         ${Boost_SYSTEM_LIBRARY})

install_headers()
install_fhicl()
install_source()
