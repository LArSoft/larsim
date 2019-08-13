# for c2: the problem is in fhiclcpp/types/detail/Sequence.icc
cet_add_compiler_flags( CXX -Wno-missing-braces )

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${MARLEY}
           ${MARLEY_ROOT}
           ${MF_MESSAGELOGGER}
           cetlib
           cetlib_except
           ${ROOT_CORE}
           ${ROOT_GEOM}
           ${ROOT_PHYSICS}
         MODULE_LIBRARIES
           larsim_EventGenerator_MARLEY
           larcoreobj_SummaryData
           lardataobj_Simulation
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ROOT_HIST}
           ${ROOT_TREE}
           ${ROOT_RIO}
         )

install_headers()
install_fhicl()
install_source()
