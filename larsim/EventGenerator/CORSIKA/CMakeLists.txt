find_library(SQLITE NAMES sqlite3 PATHS ENV SQLITE_LIB)
find_library( IFDH    NAMES ifdh    PATHS ENV IFDHC_LIB)
get_filename_component(IFDH_LIB    ${IFDH}    DIRECTORY)
include_directories("${IFDH_LIB}/../inc")

art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           art::Framework_Services_Registry
           ${CLHEP_LIBRARIES}
           ${IFDH}
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::EG
           ROOT::Physics
           ${SQLITE}
         )

install_headers()
install_fhicl()
install_source()
