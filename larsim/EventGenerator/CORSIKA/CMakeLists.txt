art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           art::Framework_Services_Registry
           ${CLHEP}
           ${IFDH}
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::EG
           ROOT::Physics
           ${SQLITE3}
         )

install_headers()
install_fhicl()
install_source()
