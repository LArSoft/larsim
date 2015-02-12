
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LIBXML2_FQ_DIR}/include/libxml2 )
include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{XERCES_C_INC} )
include_directories ( $ENV{CRYHOME}/src )

art_make( BASENAME_ONLY
          LIBRARY_NAME  LArG4
          LIB_LIBRARIES 
                        PhotonPropagation
                        PhotonVisibilityService_service
                        SpaceCharge_service
                        Utilities
			TimeService_service
                        LArProperties_service
                        DetectorProperties_service
	                Geometry
	                Geometry_service
                        ${SIMULATIONBASE}
                        ${G4BASE}
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
			${ROOT_EG}
			${ROOT_TREEPLAYER} 
			${ROOT_FFTW}
			${ROOT_REFLEX}
			${ROOT_EGPYTHIA6}
			${ROOT_GUI}
			${ROOT_BASIC_LIB_LIST}
	                ${CLHEP}
                        ${G4_LIB_LIST}
           MODULE_LIBRARIES LArG4
                        SpaceCharge_service
                        Simulation
	                Geometry
	                Geometry_service
                        ${SIMULATIONBASE}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
      G4Base      
      CLHEP
      ${G4_LIB_LIST}
		)

install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
