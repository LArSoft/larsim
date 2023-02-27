add_subdirectory(LibraryMappingTools)
add_subdirectory(ScintTimeTools)

cet_make_library(LIBRARY_NAME PhotonLibraryInterface INTERFACE
  SOURCE IPhotonLibrary.h
  LIBRARIES INTERFACE
  ROOT::Hist
)

cet_make_library(LIBRARY_NAME PhotonVisibilityTypes INTERFACE
  SOURCE PhotonVisibilityTypes.h
  LIBRARIES INTERFACE
  larsim::PhotonLibraryInterface
  larsim::PhotonMappingTransformations
)

cet_make_library(SOURCE
  PhotonLibrary.cxx
  PhotonLibraryHybrid.cxx
  PhotonPropagationUtils.cxx
  PropagationTimeModel.cxx
  SemiAnalyticalModel.cxx
  LIBRARIES
  PUBLIC
  larsim::PhotonLibraryInterface
  larcorealg::Geometry
  larcoreobj::geo_vectors
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Physics
  Boost::headers
  PRIVATE
  larsim::IonizationScintillation
  larsim::Simulation
  lardata::LArPropertiesService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  art::Framework_Services_Registry
  art_root_io::tfile_support
  messagefacility::MF_MessageLogger
  canvas::canvas
  cetlib::cetlib
  ROOT::Tree
  ROOT::MathCore
  ROOT::RooFitCore
  CLHEP::Random
  Boost::math_tr1
)

cet_build_plugin(PhotonVisibilityService art::service
  LIBRARIES
  PUBLIC
  larsim::PhotonVisibilityTypes
  larsim::PhotonLibraryInterface
  larsim::PhotonMappingTransformations
  larsim::Simulation
  larcorealg::Geometry
  larcoreobj::geo_vectors
  PRIVATE
  larsim::PhotonPropagation
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
)

cet_build_plugin(PDFastSimPAR art::EDProducer
  LIBRARIES PRIVATE
  larsim::PhotonPropagation
  larsim::IonizationScintillation
  larcore::ServiceUtil
  larcorealg::Geometry
  larcoreobj::geo_vectors
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  CLHEP::Random
)

cet_build_plugin(PDFastSimPVS art::EDProducer
  LIBRARIES PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::PhotonPropagation
  larsim::Simulation_LArG4Parameters_service
  lardata::LArPropertiesService
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art_plugin_support::toolMaker
  canvas::canvas
  fhiclcpp::fhiclcpp
  CLHEP::Random
)

cet_build_plugin(PhotonLibraryPropagation art::EDProducer
  LIBRARIES PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::IonizationScintillation
  larsim::Simulation_LArG4Parameters_service
  larcore::ServiceUtil
  lardata::LArPropertiesService
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art_plugin_support::toolMaker
  canvas::canvas
  fhiclcpp::fhiclcpp
  CLHEP::Random
)

cet_build_plugin(CreateHybridLibrary art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Simulation
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  art::Framework_Services_Registry
  ROOT::Physics
  ROOT::Tree
  ROOT::Hist
  ROOT::RIO
  ROOT::Gpad
)

cet_build_plugin(PhotonLibraryAnalyzer art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Simulation
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

add_subdirectory(LibraryBuildTools)

install_headers()
install_fhicl()
install_source()
install_scripts()
