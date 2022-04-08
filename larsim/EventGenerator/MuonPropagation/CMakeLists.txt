include_directories($ENV{MuonPropagationHOME}/src) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make(MODULE_LIBRARIES
           larcoreobj_SummaryData
           larcorealg_Geometry
           nurandom::RandomUtils_NuRandomService_service
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::EG
           ROOT::Hist
           ROOT::Physics
           ROOT::RIO
           ROOT::Tree
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
