
art_make( 
          LIBRARY_NAME  PhotonPropagation
          LIB_LIBRARIES Filters
	                RawData
	                Geometry
	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
          SERVICE_LIBRARIES PhotonPropagation
                        Simulation
                        Filters
	                RawData
	                Geometry
	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                        Simulation
                  			${ROOT_BASIC_LIB_LIST}
                        Geometry
      	                larcore_Geometry_Geometry_service
                        
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
	 )

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)

