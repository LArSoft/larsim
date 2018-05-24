art_make( 
          
          LIB_LIBRARIES 
			larsim_Simulation
			larsim_IonizationScintillation
			lardataobj_Simulation
			lardataalg_DetectorInfo
			larevt_SpaceCharge
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
           MODULE_LIBRARIES #larsim_ElectronDrift
                      larsim_Simulation lardataobj_Simulation larsim_LArG4
		      larsim_IonizationScintillation
                      lardataobj_RawData
		      larevt_SpaceCharge
                      larcorealg_Geometry
                      larcore_Geometry_Geometry_service
                      lardata_Utilities
                      nutools_RandomUtils_NuRandomService_service
                      ${ART_FRAMEWORK_CORE}
                      ${ART_FRAMEWORK_PRINCIPAL}
                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                      art_Persistency_Common
                      art_Persistency_Provenance
                      art_Utilities
#                      canvas_Utilities
                      canvas
#                      canvas_Persistency_Provenance
                      ${MF_MESSAGELOGGER}
                      ${MF_UTILITIES}
                      ${FHICLCPP}
		      ${CETLIB}
                      ${CLHEP}
                      ${ROOT_GEOM}
                      ${ROOT_GENVECTOR}
                      ${ROOT_BASIC_LIB_LIST}
                )

install_headers()
install_fhicl()
install_source()
