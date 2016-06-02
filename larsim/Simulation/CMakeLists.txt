
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( NO_PLUGINS
          
          LIB_LIBRARIES larsim_Simulation_LArVoxelCalculator_service
	                larsim_Simulation_LArG4Parameters_service
	                lardata_Utilities
	                nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(LArVoxelCalculator "service"
              
)

simple_plugin(LArG4Parameters "service"
              
)

install_headers()
install_fhicl()
install_source()
