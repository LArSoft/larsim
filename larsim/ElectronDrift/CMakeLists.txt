cet_make_library(SOURCE ISCalculationSeparate.cc
  LIBRARIES
  PUBLIC
  larcoreobj::geo_vectors
  PRIVATE
  larevt::SpaceCharge
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardataalg::DetectorInfo
  lardataobj::Simulation
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(DriftElectronstoPlane art::EDProducer
  LIBRARIES PRIVATE
  larsim::ElectronDrift
  larsim::Utils
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorPropertiesService
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  CLHEP::Random
  ROOT::MathCore
)

cet_build_plugin(ShiftEdepSCE art::EDProducer
  LIBRARIES PRIVATE
  larsim::Utils
  larsim::IonizationScintillation
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardataalg::DetectorInfo
  lardataobj::Simulation
  art_root_io::TFileService_service
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::MathCore
  ROOT::Tree
)

cet_build_plugin(SimDriftElectrons art::EDProducer
  LIBRARIES PRIVATE
  larsim::Simulation_LArG4Parameters_service
  larsim::Utils
  larsim::IonizationScintillation
  larevt::SpaceCharge
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardataalg::DetectorInfo
  larcoreobj::SimpleTypesAndConstants
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  CLHEP::Random
  ROOT::MathCore
)

cet_build_plugin(SimDriftedElectronClusterAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::Simulation
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::MathCore
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
