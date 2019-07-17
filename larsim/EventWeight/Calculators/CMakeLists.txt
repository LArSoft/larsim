include_directories($ENV{GENIE_INC}/GENIE)

art_make(LIB_LIBRARIES
         larsim_EventWeight_Base
         nugen_NuReweight_art
         nugen_NuReweight
         nurandom_RandomUtils_NuRandomService_service
         ${ART_FRAMEWORK_PRINCIPAL}
         ${ART_PERSISTENCY_PROVENANCE}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${FHICLCPP}
         cetlib_except)

install_headers()
install_fhicl()
install_source()
