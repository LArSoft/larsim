art_make( 
          
          LIB_LIBRARIES 
			larsim_Simulation
			lardataobj_Simulation
			lardataalg_DetectorInfo
			larevt_SpaceCharge
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
           MODULE_LIBRARIES larsim_IonizationScintillation
                        larsim_Simulation
			lardataobj_Simulation
			#lardata_DetectorInfoServices
			#larevt_SpaceChargeServices
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                )

install_headers()
install_fhicl()
install_source()
