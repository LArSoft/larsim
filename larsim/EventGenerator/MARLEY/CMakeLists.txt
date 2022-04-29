cet_make_library(SOURCE
  ActiveVolumeVertexSampler.cxx
  MARLEYHelper.cxx
  MarleyParameterSetWalker.h
  LIBRARIES
  PUBLIC
  nusimdata::SimulationBase
  fhiclcpp::types
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Physics
  MARLEY::MARLEY
  PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  MARLEY::MARLEY_ROOT
  ROOT::Core
  nurandom::RandomUtils_NuRandomService_service
  cetlib_except::cetlib_except
)

cet_build_plugin(MARLEYGen art::EDProducer
  LIBRARIES PRIVATE
  larsim::EventGenerator_MARLEY
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Tree
)

cet_build_plugin(MARLEYTimeGen art::EDProducer
  LIBRARIES PRIVATE
  larsim::EventGenerator_MARLEY
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  lardataobj::Simulation
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  cetlib_except::cetlib_except
  MARLEY::MARLEY
  MARLEY::MARLEY_ROOT
  ROOT::Hist
  ROOT::RIO
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
