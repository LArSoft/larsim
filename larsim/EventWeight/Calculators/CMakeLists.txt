include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(LIB_LIBRARIES larsim_EventWeight_Base
                       lardata_Utilities
                       nutools_NuReweight_art
                       nutools_NuReweight
                       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
                       nusimdata_SimulationBase
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${FHICLCPP}
                       ${CLHEP})

install_headers()
install_fhicl()
install_source()
