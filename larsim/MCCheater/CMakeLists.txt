cet_make_library(
  LIBRARY_NAME larsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc ParticleInventory.tcc
  LIBRARIES
    canvas
    nug4_ParticleNavigation
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib_except
    ${ROOT_CORE}
  )

cet_make_library(
  LIBRARY_NAME larsim_MCCheater_PhotonBackTracker
  SOURCE PhotonBackTracker.cc PhotonBackTracker.tcc
  LIBRARIES
    cetlib_except
    ${FHICLCPP}
    ${MF_MESSAGELOGGER}
    ${ROOT_CORE}
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
    ${FHICLCPP}
    ${MF_MESSAGELOGGER}
    ${ROOT_CORE}
    larcorealg_Geometry
    lardataobj_Simulation
    nug4_ParticleNavigation
  )

simple_plugin(ParticleInventoryService "service"
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_PERSISTENCY_PROVENANCE}
  larsim_MCCheater_ParticleInventory
  nug4_ParticleNavigation
  ${ROOT_CORE}
  ${ROOT_PHYSICS}
  ${MF_MESSAGELOGGER})

simple_plugin(PhotonBackTrackerService "service"
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_PERSISTENCY_PROVENANCE}
  ${MF_MESSAGELOGGER}
  larsim_MCCheater_PhotonBackTracker)

simple_plugin(BackTrackerService "service"
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_PERSISTENCY_PROVENANCE}
  ${ROOT_CORE}
  ${MF_MESSAGELOGGER}
  larsim_MCCheater_BackTracker
  larsim_MCCheater_ParticleInventory)

simple_plugin(CheckBackTracking "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ROOT_CORE})

simple_plugin(RecoCheckAna "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${MF_MESSAGELOGGER}
  ${ROOT_CORE}
  ${ROOT_HIST}
  ${ROOT_TREE})

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
