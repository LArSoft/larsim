cet_build_plugin(EventWeight art::EDProducer
  LIBRARIES PRIVATE
  larsim::EventWeight_Calculators
  larsim::EventWeight_Base
  nugen::EventGeneratorBase_GENIE
  nusimdata::SimulationBase
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
