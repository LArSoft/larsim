cet_make_library(SOURCE MergeSimSources.cxx
  LIBRARIES
  PUBLIC
  lardataobj::Simulation
  nusimdata::SimulationBase
)
cet_build_plugin(MergeSimSources art::EDProducer
  LIBRARIES PRIVATE
  larsim::MergeSimSources
  larsim::Simulation
  lardataobj::Simulation
  larcorealg::headers
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::types
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
