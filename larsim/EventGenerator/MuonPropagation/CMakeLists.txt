include_directories($ENV{MuonPropagationHOME}/src)

art_make(MODULE_LIBRARIES
           larcoreobj_SummaryData
           larcorealg_Geometry
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ${CLHEP}
           ROOT::Core
           ROOT::EG
           ROOT::Hist
           ROOT::Physics
           ROOT::RIO
           ROOT::Tree
         )

install_headers()
install_fhicl()
install_source()
