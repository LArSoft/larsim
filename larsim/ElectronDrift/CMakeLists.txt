art_make(LIB_LIBRARIES
           art::Framework_Services_Registry
           larsim_IonizationScintillation
           canvas::canvas
           cetlib_except::cetlib_except
	   messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
         MODULE_LIBRARIES
           larsim_ElectronDrift
           larsim_Simulation
	   larsim_Utils
           larsim_IonizationScintillation
           larcorealg_Geometry
           nurandom::RandomUtils_NuRandomService_service
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Tree
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


install_headers()
install_fhicl()
install_source()
