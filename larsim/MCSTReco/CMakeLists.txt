add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           art::Framework_Services_Registry
           nusimdata::SimulationBase
           lardataobj_MCBase
           canvas::canvas
           cetlib_except::cetlib_except
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           ROOT::Core
           ROOT::Physics
           art::Utilities
         MODULE_LIBRARIES
           larsim_MCSTReco
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Physics
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
