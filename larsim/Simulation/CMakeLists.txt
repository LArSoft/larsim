
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( NO_PLUGINS
          LIBRARY_NAME  Simulation
          LIB_LIBRARIES LArVoxelCalculator_service
	                LArG4Parameters_service
	                Utilities
	                ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(LArVoxelCalculator "service"
              BASENAME_ONLY
)

simple_plugin(LArG4Parameters "service"
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
