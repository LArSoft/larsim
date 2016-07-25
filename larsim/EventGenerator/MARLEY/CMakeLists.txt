
include_directories ( $ENV{MARLEY_INC} )
link_directories ( $ENV{MARLEY_LIB} )

art_make( BASENAME_ONLY
          MODULE_LIBRARIES 
                        larsim_Simulation
                        larsimobj_Simulation
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_RandomUtils_LArSeedService_service
                        MARLEY
                        MARLEY_ROOT
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
                        ${G4_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
