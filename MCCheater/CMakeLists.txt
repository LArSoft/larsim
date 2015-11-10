
simple_plugin(BackTracker "service"
			Utilities
			DetectorClocksServiceStandard_service
			Simulation
                        RecoBase
	                Geometry
	                Geometry_service
			${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(BackTrackerLoader "module"
              BackTracker_service
        			${ART_FRAMEWORK_SERVICES_REGISTRY}              
              BASENAME_ONLY
)

simple_plugin(CheckBackTracking "module"
                        BackTracker_service
                        RecoBase
                        Simulation
			${MF_MESSAGELOGGER}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
              BASENAME_ONLY
)

simple_plugin(RecoCheckAna "module"
                        BackTracker_service
                        RecoBase
      Simulation
			${SIMULATIONBASE}
			${MF_MESSAGELOGGER}
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
             BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
