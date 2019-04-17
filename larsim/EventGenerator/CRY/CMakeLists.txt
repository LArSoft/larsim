include_directories ( $ENV{CRYHOME}/src )

art_make(MODULE_LIBRARIES
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nutools_RandomUtils_NuRandomService_service
         nutools_EventGeneratorBase_CRY
         ${CRY}
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${CLHEP}
         ${MF_MESSAGELOGGER}
         ${ROOT_HIST}
         ${ROOT_PHYSICS})

install_headers()
install_fhicl()
install_source()
