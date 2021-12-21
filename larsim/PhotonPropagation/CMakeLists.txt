art_make(LIB_LIBRARIES larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       larsim_Simulation
                       larsim_IonizationScintillation
                       art::Framework_Services_Registry
                       art_root_io::tfile_support ROOT::Core
                       art_root_io::TFileService_service
                       canvas::canvas
                       messagefacility::MF_MessageLogger
                       fhiclcpp::fhiclcpp
                       cetlib_except::cetlib_except
                       ROOT::EG
                       ROOT::RIO
                       ROOT::Hist
                       ROOT::Matrix
                       ROOT::Tree
                       ROOT::GenVector
                       ROOT::RooFit
          SERVICE_LIBRARIES larsim_PhotonPropagation
                       larsim_Simulation
                       nug4::ParticleNavigation
                       lardataobj_Simulation
                       larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       art_root_io::tfile_support ROOT::Core
                       art_root_io::TFileService_service
                       messagefacility::MF_MessageLogger
                       ROOT::GenVector
          MODULE_LIBRARIES larsim_PhotonPropagation
                        larsim_LegacyLArG4
                        larsim_ElectronDrift
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_IonizationScintillation
                        larsim_Simulation
                        lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nug4::ParticleNavigation
                        art_root_io::TFileService_service
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        messagefacility::MF_MessageLogger
                        nurandom::RandomUtils_NuRandomService_service
                        art::Persistency_Provenance
                        ROOT::GenVector
                        ROOT::Gpad
         )


install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(LibraryBuildTools)
add_subdirectory(LibraryMappingTools)
add_subdirectory(ScintTimeTools)
