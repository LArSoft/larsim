foreach (dumper IN ITEMS
    GTruth
    MCParticles
    MCShowers
    MCTracks
    MCTruth
    OpDetBacktrackerRecords
    SimChannels
    SimPhotons
  )
  cet_build_plugin(Dump${dumper} art::EDAnalyzer
    LIBRARIES PRIVATE # Common dependencies.
    art::Framework_Principal
    messagefacility::MF_MessageLogger
    canvas::canvas
    fhiclcpp::types
  )
endforeach()

# Less common dependencies.
set_property(TARGET
  larsim_MCDumpers_DumpSimPhotons_module
  APPEND PROPERTY LINK_LIBRARIES
  larcorealg::headers
)

set_property(TARGET
  larsim_MCDumpers_DumpSimChannels_module
  APPEND PROPERTY LINK_LIBRARIES
  lardataobj::Simulation
)

set_property(TARGET
  larsim_MCDumpers_DumpGTruth_module
  larsim_MCDumpers_DumpMCParticles_module
  larsim_MCDumpers_DumpMCTruth_module
  APPEND PROPERTY LINK_LIBRARIES
  lardataalg::MCDumpers
)

set_property(TARGET
  larsim_MCDumpers_DumpMCParticles_module
  larsim_MCDumpers_DumpOpDetBacktrackerRecords_module
  larsim_MCDumpers_DumpSimPhotons_module
  APPEND PROPERTY LINK_LIBRARIES
  lardataobj::Simulation
)

set_property(TARGET
  larsim_MCDumpers_DumpMCShowers_module
  larsim_MCDumpers_DumpMCTracks_module
  APPEND PROPERTY LINK_LIBRARIES
  lardataobj::MCBase
)

set_property(TARGET
  larsim_MCDumpers_DumpGTruth_module
  larsim_MCDumpers_DumpMCParticles_module
  larsim_MCDumpers_DumpMCTruth_module
  APPEND PROPERTY LINK_LIBRARIES
  nusimdata::SimulationBase
)

install_headers()
install_fhicl()
install_source()
