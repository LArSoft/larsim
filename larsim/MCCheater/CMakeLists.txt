cet_make_library(
  LIBRARY_NAME larsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc ParticleInventory.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation
  nug4_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

cet_make_library(
  LIBRARY_NAME larsim_MCCheater_PhotonBackTracker
  SOURCE PhotonBackTracker.cc PhotonBackTracker.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation
  larsim_MCCheater_ParticleInventory
  larcorealg_Geometry
  nug4_ParticleNavigation
  lardataobj_Simulation
  lardataobj_RawData
  nusimdata_SimulationBase
  )

cet_make_library(
  LIBRARY_NAME larsim_MCCheater_BackTracker
  SOURCE BackTracker.cc BackTracker.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation
  larsim_MCCheater_ParticleInventory
  larcorealg_Geometry
  nug4_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  )

simple_plugin(ParticleInventoryService "service"
  larsim_MCCheater_ParticleInventory
  larsim_Simulation
  lardata_Utilities
  nug4_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

simple_plugin(PhotonBackTrackerService "service"
  larsim_MCCheater_PhotonBackTracker
  larsim_MCCheater_ParticleInventory
  larsim_MCCheater_ParticleInventoryService_service
  lardata_Utilities
  larsim_Simulation
  nug4_ParticleNavigation
  lardataobj_Simulation
  lardataobj_RecoBase
  lardataobj_RawData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

simple_plugin(BackTrackerService "service"
  larsim_MCCheater_BackTracker
  larsim_MCCheater_ParticleInventory
  larsim_MCCheater_ParticleInventoryService_service
  lardata_Utilities
  larsim_Simulation
  nug4_ParticleNavigation
  lardataobj_Simulation
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

simple_plugin(CheckBackTracking "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  larsim_Simulation
  nug4_ParticleNavigation
  lardataobj_Simulation
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  )

simple_plugin(RecoCheckAna "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  larsim_Simulation
  nug4_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Provenance
  ${ROOT_TREE}
  ${ROOT_HIST}
  )

simple_plugin(BackTrackerLoader "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  )
simple_plugin(PhotonBackTrackerLoader "module"
  larsim_MCCheater_ParticleInventoryService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  )


install_headers()
install_fhicl()
install_source()
