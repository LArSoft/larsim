add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           nusimdata_SimulationBase
           lardataobj_MCBase
           canvas
           cetlib_except
           ${MF_MESSAGELOGGER}
           ${FHICLCPP}
           ${ROOT_CORE}
           ${ROOT_PHYSICS}
         MODULE_LIBRARIES
           larsim_MCSTReco
           ${MF_MESSAGELOGGER}
           ${ROOT_CORE}
           ${ROOT_PHYSICS}
         )

install_headers()
install_fhicl()
install_source()
