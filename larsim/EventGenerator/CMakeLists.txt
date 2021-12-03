include_directories ( $ENV{NUGEN_INC} )

art_make(
         EXCLUDE
           "POTaccumulator_module.cc"
         LIB_LIBRARIES
           nusimdata::SimulationBase
           cetlib::cetlib
           fhiclcpp::fhiclcpp
           CLHEP::CLHEP
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ROOT::RIO
         MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           larsim_EventGenerator
           larsim_PhotonPropagation_PhotonVisibilityService_service
           nurandom::RandomUtils_NuRandomService_service
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::EG
           ROOT::Geom
           ROOT::Hist
           ROOT::Physics
           ROOT::RIO
        )

simple_plugin(POTaccumulator "module"
  LIBRARIES
  larcoreobj_SummaryData
  messagefacility::MF_MessageLogger
  range-v3::range-v3
  )

install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(MARLEY)
add_subdirectory(MuonPropagation)
add_subdirectory(CORSIKA)

# See ubooneoffline repository for a functioning NuWro implementation.
# Leave it commented out here, as below.
##add_subdirectory(NuWro)
