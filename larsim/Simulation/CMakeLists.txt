art_make(NO_PLUGINS
         LIB_LIBRARIES
           larsim_Simulation_LArVoxelCalculator_service
           larsim_Simulation_LArG4Parameters_service
           lardataobj_Simulation
           lardata_Utilities
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           art::Framework_Principal
           art::Persistency_Provenance
           ROOT::Core
           ROOT::Physics
           ROOT::GenVector
           ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(LArVoxelCalculator "service") ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(LArG4Parameters "service") ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
