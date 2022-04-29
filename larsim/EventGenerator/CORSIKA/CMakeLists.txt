cet_build_plugin(CORSIKAGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  larcoreobj::SummaryData
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ifdh::ifdh
  messagefacility::MF_MessageLogger
  ROOT::Core
  ROOT::EG
  CLHEP::Random
  SQLite::SQLite3
)

install_headers()
install_fhicl()
install_source()
