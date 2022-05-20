cet_make_library(LIBRARY_NAME larsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc
  LIBRARIES
  PUBLIC
  lardataobj::Simulation
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  canvas::canvas
  fhiclcpp::types
  PRIVATE
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_make_library(LIBRARY_NAME larsim_MCCheater_PhotonBackTracker
  SOURCE PhotonBackTracker.cc
  LIBRARIES
  PUBLIC
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  lardataobj::Simulation
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::types
  PRIVATE
  larsim::MCCheater_ParticleInventory
  larcorealg::Geometry
)

cet_make_library(LIBRARY_NAME larsim_MCCheater_BackTracker
  SOURCE BackTracker.cc
  LIBRARIES
  PUBLIC
  larsim::MCCheater_ParticleInventory
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  lardataobj::Simulation
  canvas::canvas
  fhiclcpp::types
  PRIVATE
  larcorealg::Geometry
  lardataalg::DetectorInfo
)

cet_build_plugin(ParticleInventoryService art::service
  LIBRARIES
  PUBLIC
  larsim::MCCheater_ParticleInventory
  nusimdata::SimulationBase
  art::Persistency_Provenance
  fhiclcpp::types
  PRIVATE
  nug4::ParticleNavigation
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(PhotonBackTrackerService art::service
  LIBRARIES
  PUBLIC
  larsim::MCCheater_PhotonBackTracker
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  lardataobj::Simulation
  art::Persistency_Provenance
  canvas::canvas
  fhiclcpp::types
  PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larsim::MCCheater_ParticleInventory
  larcore::ServiceUtil
  larcorealg::Geometry
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(BackTrackerService art::service
  LIBRARIES
  PUBLIC
  larsim::MCCheater_BackTracker
  larcore::ServiceUtil
  lardataobj::RecoBase
  lardataobj::Simulation
  art::Persistency_Provenance
  canvas::canvas
  fhiclcpp::types
  PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larsim::MCCheater_ParticleInventory
  lardata::DetectorClocksService
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(CheckBackTracking art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  nug4::ParticleNavigation
  art::Framework_Services_Registry
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  canvas::canvas
)

cet_build_plugin(RecoCheckAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  lardataobj::Simulation
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  canvas::canvas
  messagefacility::MF_MessageLogger
  ROOT::Core
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(BackTrackerLoader art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  art::Framework_Services_Registry
)

cet_build_plugin(PhotonBackTrackerLoader art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larsim::MCCheater_PhotonBackTrackerService_service
  art::Framework_Services_Registry
)

install_headers()
install_fhicl()
install_source()
