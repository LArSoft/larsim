cet_make_library(
  LIBRARY_NAME larsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc ParticleInventory.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation 
  nutools_ParticleNavigation 
  lardataobj_Simulation
  nusimdata_SimulationBase
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
  nutools_ParticleNavigation 
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
  nutools_ParticleNavigation 
  lardataobj_Simulation
  nusimdata_SimulationBase
  ) 

simple_plugin(ParticleInventoryService "service"
  larsim_MCCheater_ParticleInventory
  larsim_Simulation 
  lardata_Utilities
  nutools_ParticleNavigation 
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  )

simple_plugin(PhotonBackTrackerService "service"
  larsim_MCCheater_PhotonBackTracker
  larsim_MCCheater_ParticleInventory
  larsim_MCCheater_ParticleInventoryService_service
  lardata_Utilities
  larsim_Simulation 
  nutools_ParticleNavigation 
  lardataobj_Simulation
  lardataobj_RecoBase
  lardataobj_RawData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  )

simple_plugin(BackTrackerService "service"
  larsim_MCCheater_BackTracker
  larsim_MCCheater_ParticleInventory
  larsim_MCCheater_ParticleInventoryService_service
  lardata_Utilities
  larsim_Simulation 
  nutools_ParticleNavigation 
  lardataobj_Simulation
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  )

simple_plugin(CheckBackTracking "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  cetlib cetlib_except

  )

simple_plugin(RecoCheckAna "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}

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
