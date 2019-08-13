art_make(NO_PLUGINS
         LIB_LIBRARIES
           larsim_Simulation_LArVoxelCalculator_service
           larsim_Simulation_LArG4Parameters_service
           lardataobj_Simulation
           lardata_Utilities
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_PERSISTENCY_PROVENANCE}
           ${ROOT_CORE}
           ${ROOT_PHYSICS})

simple_plugin(LArVoxelCalculator "service")
simple_plugin(LArG4Parameters "service")

install_headers()
install_fhicl()
install_source()
