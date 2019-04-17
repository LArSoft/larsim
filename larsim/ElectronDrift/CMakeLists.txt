art_make(MODULE_LIBRARIES
         larsim_Simulation
         lardataobj_Simulation
         larsim_LArG4
         larsim_IonizationScintillation
         lardataobj_RawData
         larevt_SpaceCharge
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         lardata_Utilities
         nutools_RandomUtils_NuRandomService_service
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${ROOT_GEOM}
         ${ROOT_GENVECTOR})

install_headers()
install_fhicl()
install_source()
