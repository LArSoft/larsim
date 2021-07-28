add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           art::Framework_Services_Registry
           nusimdata_SimulationBase
           lardataobj_MCBase
           canvas::canvas
           cetlib_except::cetlib_except
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           ROOT::Core
           ROOT::Physics
           ${ART_UTILITIES}
         MODULE_LIBRARIES
           larsim_MCSTReco
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Physics
         )

install_headers()
install_fhicl()
install_source()
