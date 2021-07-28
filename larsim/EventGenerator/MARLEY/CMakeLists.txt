# for c2: the problem is in fhiclcpp/types/detail/Sequence.icc
cet_add_compiler_flags( CXX -Wno-missing-braces )

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${MARLEY}
           ${MARLEY_ROOT}
           messagefacility::MF_MessageLogger
           cetlib::cetlib
           cetlib_except::cetlib_except
           ROOT::Core
           ROOT::Geom
           ROOT::Physics
           ${ART_UTILITIES}
           fhiclcpp::fhiclcpp
         MODULE_LIBRARIES
           larsim_EventGenerator_MARLEY
           larcoreobj_SummaryData
           lardataobj_Simulation
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           ROOT::Hist
           ROOT::Tree
           ROOT::RIO
         )

install_headers()
install_fhicl()
install_source()
