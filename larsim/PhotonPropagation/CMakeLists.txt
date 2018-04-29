
art_make( 
          
          LIB_LIBRARIES larevt_Filters
	                lardataobj_RawData
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			cetlib cetlib_except
			${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
          SERVICE_LIBRARIES larsim_PhotonPropagation
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larevt_Filters
	                lardataobj_RawData
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larsim_PhotonPropagation
                       larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_IonizationScintillation
                        lardataobj_Simulation                   
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larcorealg_Geometry
      	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			cetlib cetlib_except
                  	${ROOT_BASIC_LIB_LIST}
	 )

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)

