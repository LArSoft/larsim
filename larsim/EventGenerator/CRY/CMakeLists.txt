art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom::RandomUtils_NuRandomService_service
           nusimdata::SimulationBase
           nutools::EventGeneratorBase_CRY
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
