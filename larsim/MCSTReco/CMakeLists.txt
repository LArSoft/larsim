
add_subdirectory(job)

art_make( 
	  
          LIB_LIBRARIES    larcore_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation
			   SimulationBase
			   larcore_SummaryData
			   lardata_Utilities
			   lardata_MCBase
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
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
			   ${ROOT_GEOM}
			   ${ROOT_XMLIO}
			   ${ROOT_GDML}
			   ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larsim_MCSTReco
			   larcore_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation
			   SimulationBase
			   larcore_SummaryData
			   lardata_Utilities
			   lardata_MCBase
			   SeedService_service # artextensions
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
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
			   ${ROOT_GEOM}
			   ${ROOT_XMLIO}
			   ${ROOT_GDML}
			   ${ROOT_BASIC_LIB_LIST}
		  	   ${ROOT_EG}
			   ${CLHEP}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			)

install_headers()
install_fhicl()
install_source()
