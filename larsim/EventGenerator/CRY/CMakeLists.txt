include_directories($ENV{CRYHOME}/src)

art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           nutools_EventGeneratorBase_CRY
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
         )

install_headers()
install_fhicl()
install_source()
