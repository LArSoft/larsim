art_make(LIB_LIBRARIES
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           larsim_IonizationScintillation
           canvas
           cetlib_except
	   ${MF_MESSAGELOGGER}
           ${FHICLCPP}
         MODULE_LIBRARIES
           larsim_ElectronDrift
           larsim_Simulation
	   larsim_Utils
           larsim_IonizationScintillation
           larcorealg_Geometry
           nurandom_RandomUtils_NuRandomService_service
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${CLHEP}
           ${MF_MESSAGELOGGER}
           ROOT::Core
           ROOT::Tree
         )


install_headers()
install_fhicl()
install_source()
