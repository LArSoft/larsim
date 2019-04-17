art_make(LIB_LIBRARIES
           larsim_Simulation
           lardataobj_Simulation
           lardataalg_DetectorInfo
           larevt_SpaceCharge
           nusimdata_SimulationBase
           ${FHICLCPP}
           ${MF_MESSAGELOGGER}
         MODULE_LIBRARIES
           larsim_IonizationScintillation
           larsim_Simulation
           lardataobj_Simulation
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ROOT_TREE}
           ${MF_MESSAGELOGGER})

install_headers()
install_fhicl()
install_source()
