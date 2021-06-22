add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           art::Framework_Services_Registry
           nusimdata_SimulationBase
           lardataobj_MCBase
           canvas
           cetlib_except
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp fhiclcpp::types
           ROOT::Core
           ROOT::Physics
           art::Utilities
         MODULE_LIBRARIES
           larsim_MCSTReco
           messagefacility::headers
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Physics
         )

install_headers()
install_fhicl()
install_source()
