cet_build_plugin(GaisserParam art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::EG
  ROOT::Hist
  ROOT::Physics
  ROOT::Tree
  CLHEP::Random
)

cet_build_plugin(MUSUN art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::EG
  ROOT::Tree
  CLHEP::Random
)

install_headers()
install_fhicl()
install_source()
