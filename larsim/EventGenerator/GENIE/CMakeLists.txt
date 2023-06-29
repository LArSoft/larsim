cet_build_plugin(GENIEGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::Utilities
  larcoreobj::SummaryData
  lardataalg::MCDumpers
  lardataobj::Simulation
  nugen::EventGeneratorBase_GENIE
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  GENIE::GFwEG
  GENIE::GFwAlg
  GENIE::GFwMsg
  GENIE::GFwGHEP
  GENIE::GTlFlx
  dk2nu::Genie
  dk2nu::Tree
  ROOT::Core
  ROOT::EG
  ROOT::Hist
)

cet_build_plugin(NeutronOsc art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nugen::EventGeneratorBase_GENIE
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  CLHEP::Random
  ${GENIE_LIB_LIST}
  log4cpp::log4cpp
)

cet_build_plugin(NucleonDecay art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nugen::EventGeneratorBase_GENIE
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  CLHEP::Random
  ${GENIE_LIB_LIST}
  log4cpp::log4cpp
)

install_headers()
install_fhicl()
install_source()
