
simple_plugin(BackTracker "service"
			lardata_Utilities
			larsim_Simulation
                        lardata_RecoBase
	                larcore_Geometry
	                larcore_Geometry_Geometry_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
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
			nusimdata_SimulationBase
			${MF_MESSAGELOGGER}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
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
			${ROOT_BASIC_LIB_LIST}
             
)

install_headers()
install_fhicl()
install_source()
