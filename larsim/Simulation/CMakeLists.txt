
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( NO_PLUGINS
          
          LIB_LIBRARIES larsim_Simulation_LArVoxelCalculator_service
	                larsim_Simulation_LArG4Parameters_service
			lardataobj_Simulation
	                lardata_Utilities
	                nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(LArVoxelCalculator "service"
			cetlib_except
)

simple_plugin(LArG4Parameters "service"
			cetlib_except
)

install_headers()
install_fhicl()
install_source()
