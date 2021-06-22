cet_make_library(
  LIBRARY_NAME larsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc ParticleInventory.tcc
  LIBRARIES
    canvas
    nug4_ParticleNavigation
    messagefacility::MF_MessageLogger
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
    messagefacility::MF_MessageLogger
    ROOT::Core
    canvas
    larsim_MCCheater_ParticleInventory
    larcorealg_Geometry
    lardataobj_Simulation
    nug4_ParticleNavigation
  )

cet_make_library(
  LIBRARY_NAME larsim_MCCheater_BackTracker
  SOURCE BackTracker.cc BackTracker.tcc
  LIBRARIES
    canvas
    cetlib_except
    fhiclcpp::fhiclcpp fhiclcpp::types
    messagefacility::MF_MessageLogger
    ROOT::Core
    larcorealg_Geometry
    lardataobj_Simulation
    nug4_ParticleNavigation
  )

simple_plugin(ParticleInventoryService "service"
  art::Framework_Principal
  art::Persistency_Provenance
  larsim_MCCheater_ParticleInventory
  nug4_ParticleNavigation
  ROOT::Core
  ROOT::Physics
  messagefacility::MF_MessageLogger)

simple_plugin(PhotonBackTrackerService "service"
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  larsim_MCCheater_PhotonBackTracker)

simple_plugin(BackTrackerService "service"
  art::Framework_Principal
  art::Persistency_Provenance
  ROOT::Core
  messagefacility::MF_MessageLogger
  larsim_MCCheater_BackTracker
  nug4_ParticleNavigation
  larsim_MCCheater_ParticleInventory)

simple_plugin(CheckBackTracking "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  messagefacility::MF_MessageLogger
  art::Framework_Services_Registry
  nug4_ParticleNavigation
  ROOT::Core)

simple_plugin(RecoCheckAna "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nusimdata_SimulationBase
  nug4_ParticleNavigation
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  art_root_io::tfile_support
  messagefacility::MF_MessageLogger
  ROOT::Core
  ROOT::Hist
  ROOT::Tree)

simple_plugin(BackTrackerLoader "module"
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  art::Framework_Services_Registry
  nug4_ParticleNavigation
  ROOT::Core)

simple_plugin(PhotonBackTrackerLoader "module"
  larsim_MCCheater_ParticleInventoryService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  nug4_ParticleNavigation
  ROOT::Core
  art::Framework_Services_Registry)

install_headers()
install_fhicl()
install_source()
