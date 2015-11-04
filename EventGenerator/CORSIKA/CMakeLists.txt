
cet_find_library( SQLITE3 NAMES sqlite3 PATHS ENV NO_DEFAULT_PATH )
art_make( BASENAME_ONLY
          MODULE_LIBRARIES 
                        SummaryData
                        Geometry
                        Geometry_service
                        SeedService_service # artextensions
                        Simulation
                        ${SIMULATIONBASE}
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
                        ${SQLITE3}
                        ${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()

