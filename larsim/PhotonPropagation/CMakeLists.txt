art_make(LIB_LIBRARIES larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       larsim_Simulation
                       larsim_IonizationScintillation
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
                       ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                       canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CLHEP}
                       cetlib_except
                       ROOT::EG
                       ROOT::RIO
                       ROOT::Hist
                       ROOT::Matrix
                       ROOT::Tree
                       ROOT::GenVector
          SERVICE_LIBRARIES larsim_PhotonPropagation
                       larsim_Simulation
                       nug4_ParticleNavigation
                       lardataobj_Simulation
                       larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
                       ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                       ${MF_MESSAGELOGGER}
                       ROOT::GenVector
          MODULE_LIBRARIES larsim_PhotonPropagation
                        larsim_LArG4
                        larsim_ElectronDrift
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_IonizationScintillation
                        larsim_Simulation
                        lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nug4_ParticleNavigation
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ROOT::Core
                        ${MF_MESSAGELOGGER}
                        nurandom_RandomUtils_NuRandomService_service
                        art_Persistency_Provenance
                        ${CLHEP}
                        ROOT::GenVector
                        ROOT::Gpad
         )


install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(LibraryBuildTools)
add_subdirectory(LibraryMappingTools)
add_subdirectory(ScintTimeTools)