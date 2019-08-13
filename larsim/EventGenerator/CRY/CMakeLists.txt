include_directories($ENV{CRYHOME}/src)

art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           nutools_EventGeneratorBase_CRY
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${MF_MESSAGELOGGER}
           ${ROOT_CORE}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
         )

install_headers()
install_fhicl()
install_source()
