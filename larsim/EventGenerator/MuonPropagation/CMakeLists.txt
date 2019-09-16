include_directories($ENV{MuonPropagationHOME}/src)

art_make(MODULE_LIBRARIES
           larcoreobj_SummaryData
           larcorealg_Geometry
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${MF_MESSAGELOGGER}
           ${CLHEP}
           ${ROOT_CORE}
           ${ROOT_EG}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
           ${ROOT_RIO}
           ${ROOT_TREE}
         )

install_headers()
install_fhicl()
install_source()
