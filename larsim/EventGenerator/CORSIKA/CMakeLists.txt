cet_find_simple_package( ifdh INCPATH_SUFFIXES inc INCPATH_VAR IFDHC_INC )
include_directories(${IFDHC_INC})
art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${CLHEP_LIBRARIES}
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
