add_subdirectory(job)

cet_make_library(SOURCE
  MCRecoEdep.cxx
  MCRecoPart.cxx
  MCShowerRecoAlg.cxx
  MCShowerRecoPart.cxx
  MCTrackCollectionAnaAlg.cxx
  MCTrackRecoAlg.cxx
  LIBRARIES
  PUBLIC
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::Simulation
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  cetlib_except::cetlib_except
  ROOT::Physics
  PRIVATE
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Core
  ROOT::MathCore
  ROOT::Tree
)

cet_build_plugin(MCReco art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCSTReco
  lardataobj::MCBase
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
