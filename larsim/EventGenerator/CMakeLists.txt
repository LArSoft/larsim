art_make(
         EXCLUDE
           "POTaccumulator_module.cc"
         LIB_LIBRARIES
           nusimdata_SimulationBase
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
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           ${CLHEP}
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::EG
           ROOT::Geom
           ROOT::Hist
           ROOT::Physics
           ROOT::RIO
        )

simple_plugin(POTaccumulator "module"
  larcoreobj_SummaryData
  messagefacility::MF_MessageLogger
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
