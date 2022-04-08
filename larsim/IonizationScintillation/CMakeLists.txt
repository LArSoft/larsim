art_make(LIB_LIBRARIES
           larsim_Simulation
           lardataobj_Simulation
           lardataalg_DetectorInfo
           larevt_SpaceCharge
           nusimdata::SimulationBase
           larcorealg_Geometry
           larcore_Geometry_Geometry_service
           fhiclcpp::fhiclcpp
           messagefacility::MF_MessageLogger
         MODULE_LIBRARIES
           larsim_IonizationScintillation
           larsim_Simulation
           lardataobj_Simulation
           nurandom::RandomUtils_NuRandomService_service
           larcorealg_Geometry
           larcore_Geometry_Geometry_service
           art::Framework_Services_Registry
           art_root_io::tfile_support
           art_root_io::TFileService_service
           ROOT::Core
           ROOT::Tree
           messagefacility::MF_MessageLogger) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
