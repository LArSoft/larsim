cet_build_plugin(LArVoxelCalculator art::service
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
)

cet_build_plugin(LArG4Parameters art::service
  LIBRARIES
  PUBLIC
  larcoreobj::SimpleTypesAndConstants
  PRIVATE
  fhiclcpp::fhiclcpp
)

cet_make_library(SOURCE
	PhotonVoxels.cxx
  LArVoxelData.cxx
  LArVoxelID.cxx
  LArVoxelList.cxx
  SimListUtils.cxx
  LIBRARIES
  PUBLIC
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  lardataobj::Simulation
  ROOT::Physics
  PRIVATE
  larsim::Simulation_LArVoxelCalculator_service
  larsim::Simulation_LArG4Parameters_service
  lardata::DetectorClocksService
  art::Framework_Principal
  art::Framework_Services_Registry
)

install_headers()
install_fhicl()
install_source()
