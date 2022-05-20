cet_make_library(SOURCE
  ISCalcCorrelated.cxx
  ISCalcNESTLAr.cxx
  ISCalcSeparate.cxx
  ISTPC.cxx
  LIBRARIES
  PUBLIC
  larcorealg::Geometry
  larcoreobj::geo_vectors
  PRIVATE
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::LArPropertiesService
  lardataalg::DetectorInfo
  lardataobj::Simulation
  messagefacility::MF_MessageLogger
  CLHEP::Random
  CLHEP::Vector
)

cet_build_plugin(IonAndScint art::EDProducer
  LIBRARIES PRIVATE
  larsim::IonizationScintillation
  larevt::SpaceChargeService
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  fhiclcpp::fhiclcpp
  canvas::canvas
  CLHEP::Random
)

cet_build_plugin(ISCalcAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::Simulation_LArG4Parameters_service
  larevt::SpaceCharge
  larsim::IonizationScintillation
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
  CLHEP::Random
)

install_headers()
install_fhicl()
install_source()
