
include_directories ( $ENV{GaisserParamHOME}/src )

art_make( BASENAME_ONLY
          MODULE_LIBRARIES 
	  		PhotonPropagation
                        PhotonVisibilityService_service
                        Simulation
                        SummaryData
                        Geometry
                        Geometry_service
                        SeedService_service # artextensions
                        ${GaisserParam}
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
			${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()

