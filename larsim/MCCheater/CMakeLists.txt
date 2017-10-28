
simple_plugin(PhotonBackTracker "service"
			lardata_Utilities
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardataobj_RecoBase
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
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

simple_plugin(BackTracker "service"
			lardata_Utilities
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardataobj_RecoBase
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
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

simple_plugin(PhotonBackTrackerLoader "module"
              larsim_MCCheater_PhotonBackTracker_service
        			${ART_FRAMEWORK_SERVICES_REGISTRY}              
              )
              
simple_plugin(BackTrackerLoader "module"
              larsim_MCCheater_BackTracker_service
        			${ART_FRAMEWORK_SERVICES_REGISTRY}              
              
)

simple_plugin(CheckBackTracking "module"
                        larsim_MCCheater_BackTracker_service
                        lardataobj_RecoBase
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			${MF_MESSAGELOGGER}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			cetlib cetlib_except
              
)

simple_plugin(RecoCheckAna "module"
                        larsim_MCCheater_BackTracker_service
                        lardataobj_RecoBase
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			nusimdata_SimulationBase
			${MF_MESSAGELOGGER}
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
