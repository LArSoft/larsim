cet_build_plugin(SimWireAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::RawData
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(SimWire art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::Utilities_LArFFT_service
  larcorealg::Geometry
  lardataalg::DetectorInfo
  lardataobj::RawData
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  CLHEP::Random
  ROOT::MathCore
)

cet_build_plugin(WienerFilterAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::Utilities_LArFFT_service
  larcore::Geometry_Geometry_service
  lardata::Utilities_LArFFT_service
  lardataalg::DetectorInfo
  lardataobj::RawData
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

install_headers()
install_fhicl()
install_source()
