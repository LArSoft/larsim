cet_make_library(SOURCE NueAr40CCGenerator.cxx
  LIBRARIES PRIVATE
  nusimdata::SimulationBase
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::Hist
  ROOT::MathCore
  ROOT::Physics
  ROOT::Tree
  CLHEP::Random
)

cet_build_plugin(FileMuons art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::EG
  ROOT::Hist
  ROOT::Physics
  ROOT::RIO
)

cet_build_plugin(LightSource art::EDProducer
  LIBRARIES PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Simulation
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcorealg::Geometry
  larcoreobj::geo_vectors
  larcoreobj::SummaryData
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::Core
  ROOT::Geom
  ROOT::Physics
  ROOT::Tree
  CLHEP::Random
)

cet_build_plugin(NDKGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::EG
  ROOT::Tree
  CLHEP::Random
)

cet_build_plugin(POTaccumulator art::EDAnalyzer
  LIBRARIES PRIVATE
  larcoreobj::SummaryData
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  range-v3::range-v3
)

cet_build_plugin(RadioGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataalg::DetectorInfo
  larcorealg::Geometry
  larcorealg::headers
  larcoreobj::SummaryData
  larcoreobj::geo_vectors
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)

# For evgenbase.h:
if (TARGET nugen::EventGeneratorBase)
  target_link_libraries(larsim_EventGenerator_RadioGen_module
    PRIVATE
    nugen::EventGeneratorBase
  )
else()
  target_link_libraries(larsim_EventGenerator_RadioGen_module
    PRIVATE
    ROOT::MathCore
  )
  target_include_directories(larsim_EventGenerator_RadioGen_module
    BEFORE PRIVATE ${nugen_INCLUDE_DIR}
  )
endif()

cet_build_plugin(SNNueAr40CCGen art::EDProducer
  LIBRARIES PRIVATE
  larsim::EventGenerator
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
)

cet_build_plugin(SingleGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcorealg::Geometry
  larcoreobj::SummaryData
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::EG
  ROOT::Hist
  ROOT::Physics
  CLHEP::Random
)

cet_build_plugin(TextFileGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Physics
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(MARLEY)
add_subdirectory(MuonPropagation)
add_subdirectory(CORSIKA)

# See ubooneoffline repository for a functioning NuWro implementation.
# Leave it commented out here, as below.
##add_subdirectory(NuWro)
