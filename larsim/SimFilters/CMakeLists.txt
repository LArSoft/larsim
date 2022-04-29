cet_build_plugin(FilterCryostatNus art::EDFilter
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcorealg::Geometry
  nusimdata::SimulationBase
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(FilterGenInTime art::EDFilter
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  larcorealg::Geometry
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(FilterNoDirtNeutrinos art::EDFilter
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(FilterNoMCParticles art::SharedFilter
  LIBRARIES PRIVATE
  nusimdata::SimulationBase
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(FilterPrimaryPDG art::EDFilter
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  nug4::ParticleNavigation
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(FilterSimPhotonLiteTime art::SharedFilter
  LIBRARIES PRIVATE
  lardataobj::Simulation
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(FilterSimPhotonTime art::SharedFilter
  LIBRARIES PRIVATE
  lardataobj::Simulation
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(FilterStoppingMuon art::SharedFilter
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_fhicl()
install_source()
