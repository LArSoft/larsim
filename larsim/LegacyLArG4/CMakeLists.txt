cet_make_library(SOURCE
  AllPhysicsLists.cc
  AuxDetReadout.cxx
  AuxDetReadoutGeometry.cxx
  CustomPhysicsTable.cxx
  FastOpticalPhysics.cxx
  G4BadIdeaAction.cxx
  G4ThermalElectron.cxx
  ISCalculation.cxx
  ISCalculationCorrelated.cxx
  ISCalculationNEST.cxx
  ISCalculationSeparate.cxx
  IonizationAndScintillation.cxx
  IonizationAndScintillationAction.cxx
  LArStackingAction.cxx
  LArVoxelReadout.cxx
  LArVoxelReadoutGeometry.cxx
  MaterialPropertyLoader.cxx
  MuNuclearSplittingProcess.cxx
  MuNuclearSplittingProcessXSecBias.cxx
  NestAlg.cxx
  NeutronHPphysics.cc
  OpBoundaryProcessSimple.cxx
  OpDetLookup.cxx
  OpDetPhotonTable.cxx
  OpDetReadoutGeometry.cxx
  OpDetSensitiveDetector.cxx
  OpFastScintillation.cxx
  OpParamAction.cxx
  OpParamSD.cxx
  OpticalPhysics.cxx
  ParticleListAction.cxx
  PhysicsList.cxx
  VisualizationAction.cxx
  LIBRARIES
  PUBLIC
  larsim::PhotonVisibilityTypes
  larsim::Simulation_LArG4Parameters_service
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcorealg::Geometry
  larcoreobj::geo_vectors
  lardataobj::Simulation
  nug4::G4Base
  nusimdata::SimulationBase
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  cetlib::cetlib
  cetlib_except::cetlib_except
  Geant4::G4digits_hits
  Geant4::G4event
  Geant4::G4geometry
  Geant4::G4global
  Geant4::G4materials
  Geant4::G4particles
  Geant4::G4physicslists
  Geant4::G4processes
  Geant4::G4run
  Geant4::G4track
  nug4::ParticleNavigation
  ROOT::Hist
  ROOT::Physics
  PRIVATE
  larsim::Simulation_LArVoxelCalculator_service
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Utils
  larevt::SpaceCharge
  larevt::SpaceChargeService
  lardata::LArPropertiesService
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::Simulation
  art_root_io::TFileService_service
  ROOT::Core
  ROOT::Geom
  ROOT::MathCore
  CLHEP::Geometry
  CLHEP::Random
  CLHEP::Vector
  Boost::math_tr1
)

cet_build_plugin(LArG4 art::EDProducer
  LIBRARIES PRIVATE
  larsim::LegacyLArG4
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Simulation_LArG4Parameters_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  larcorealg::headers
  lardataalg::MCDumpers
  lardataobj::Simulation
  nug4::G4Base
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  Geant4::G4digits_hits
  Geant4::G4geometry
  Geant4::G4run
)

cet_build_plugin(LArG4Ana art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  nug4::ParticleNavigation
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Physics
  ROOT::Tree
)

cet_build_plugin(LArSimChannelAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::Simulation
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Hist
)

install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
