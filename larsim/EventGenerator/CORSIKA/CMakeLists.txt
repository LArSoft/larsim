art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${CLHEP}
           ${IFDH}
           ${MF_MESSAGELOGGER}
           ROOT::Core
           ROOT::EG
           ROOT::Physics
           ${SQLITE3}
         )

install_headers()
install_fhicl()
install_source()
