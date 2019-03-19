art_make(LIB_LIBRARIES larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       larsim_Simulation
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CLHEP}
                       cetlib_except
                       ${ROOT_EG}
                       ${ROOT_RIO}
                       ${ROOT_HIST}
                       ${ROOT_MATRIX}
                       ${ROOT_TREE}
          SERVICE_LIBRARIES larsim_PhotonPropagation
                       larsim_Simulation
                       nutools_ParticleNavigation
                       lardataobj_Simulation
                       larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       ${MF_MESSAGELOGGER}
          MODULE_LIBRARIES larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_IonizationScintillation
                        lardataobj_Simulation
                        larsim_Simulation
                        nutools_ParticleNavigation
                        lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        nutools_RandomUtils_NuRandomService_service
                        art_Persistency_Provenance
                        ${MF_MESSAGELOGGER}
                        ${CLHEP}
                        ${ROOT_GPAD}
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)
