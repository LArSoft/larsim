
include_directories ( $ENV{MuonPropagationHOME}/src )

art_make( 
          MODULE_LIBRARIES 
	  		larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_Simulation
                        larcore_SummaryData
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        SeedService_service # artextensions
                        ${GaisserParam}
			${MUSUN}
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

