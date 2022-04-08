art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           nurandom::RandomUtils_NuRandomService_service
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           CLHEP::CLHEP
           ifdh::ifdh
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::EG
           ROOT::Physics
           SQLite::SQLite3
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
