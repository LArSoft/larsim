art_make(LIB_LIBRARIES
           ${ART_UTILITIES}
	   ${MF_MESSAGELOGGER}
         MODULE_LIBRARIES
           larsim_IonizationScintillation
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ROOT_TREE}
           ${ROOT_CORE}
         )

install_headers()
install_fhicl()
install_source()
