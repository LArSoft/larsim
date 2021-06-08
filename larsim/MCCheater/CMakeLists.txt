cet_make_library(
  LIBRARY_NAME larsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc ParticleInventory.tcc
  LIBRARIES
    canvas
    nug4_ParticleNavigation
    messagefacility::MFMessageLogger
    fhiclcpp::fhiclcpp fhiclcpp::types
    cetlib_except
    ROOT::Core
  )

cet_make_library(
  LIBRARY_NAME larsim_MCCheater_PhotonBackTracker
  SOURCE PhotonBackTracker.cc PhotonBackTracker.tcc
  LIBRARIES
    cetlib_except
    fhiclcpp::fhiclcpp fhiclcpp::types
    messagefacility::MFMessageLogger
    ROOT::Core
    canvas
    larsim_MCCheater_ParticleInventory
    larcorealg_Geometry
    lardataobj_Simulation
  )

cet_make_library(
  LIBRARY_NAME larsim_MCCheater_BackTracker
  SOURCE BackTracker.cc BackTracker.tcc
  LIBRARIES
    canvas
    cetlib_except
    fhiclcpp::fhiclcpp fhiclcpp::types
    messagefacility::MFMessageLogger
    ROOT::Core
    larcorealg_Geometry
    lardataobj_Simulation
    nug4_ParticleNavigation
  )

simple_plugin(ParticleInventoryService "service"
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_PERSISTENCY_PROVENANCE}
  larsim_MCCheater_ParticleInventory
  nug4_ParticleNavigation
  ROOT::Core
  ROOT::Physics
  messagefacility::MFMessageLogger)

simple_plugin(PhotonBackTrackerService "service"
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_PERSISTENCY_PROVENANCE}
  messagefacility::MFMessageLogger
  larsim_MCCheater_PhotonBackTracker)

simple_plugin(BackTrackerService "service"
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_PERSISTENCY_PROVENANCE}
  ROOT::Core
  messagefacility::MFMessageLogger
  larsim_MCCheater_BackTracker
  larsim_MCCheater_ParticleInventory)

simple_plugin(CheckBackTracking "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  messagefacility::MFMessageLogger
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ROOT::Core)

simple_plugin(RecoCheckAna "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  messagefacility::MFMessageLogger
  ROOT::Core
  ROOT::Hist
  ROOT::Tree)

simple_plugin(BackTrackerLoader "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY})

simple_plugin(PhotonBackTrackerLoader "module"
  larsim_MCCheater_ParticleInventoryService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY})

install_headers()
install_fhicl()
install_source()
