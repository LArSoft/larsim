
add_subdirectory(job)

art_make( 
	  
          LIB_LIBRARIES    larcore_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			   nusimdata_SimulationBase
			   larcoreobj_SummaryData
			   lardata_Utilities
			   lardataobj_MCBase
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   art_Persistency_Common canvas_Persistency_Common
			   art_Persistency_Provenance canvas_Persistency_Provenance
			   art_Utilities canvas_Utilities
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
			   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			   nusimdata_SimulationBase
			   larcoreobj_SummaryData
			   lardata_Utilities
			   lardataobj_MCBase
			   nutools_RandomUtils_NuRandomService_service
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   art_Persistency_Common canvas_Persistency_Common
			   art_Persistency_Provenance canvas_Persistency_Provenance
			   art_Utilities canvas_Utilities
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
