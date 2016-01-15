
simple_plugin(BackTracker "service"
			lardata_Utilities
			larsim_Simulation
                        lardata_RecoBase
	                larcore_Geometry
	                larcore_Geometry_Geometry_service
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
              
)

simple_plugin(BackTrackerLoader "module"
              larsim_MCCheater_BackTracker_service
        			${ART_FRAMEWORK_SERVICES_REGISTRY}              
              
)

simple_plugin(CheckBackTracking "module"
                        larsim_MCCheater_BackTracker_service
                        lardata_RecoBase
                        larsim_Simulation
			${MF_MESSAGELOGGER}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
              
)

simple_plugin(RecoCheckAna "module"
                        larsim_MCCheater_BackTracker_service
                        lardata_RecoBase
      larsim_Simulation
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
             
)

install_headers()
install_fhicl()
install_source()
