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
           nurandom_RandomUtils_NuRandomService_service
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           ${CLHEP}
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Tree
         )


install_headers()
install_fhicl()
install_source()
