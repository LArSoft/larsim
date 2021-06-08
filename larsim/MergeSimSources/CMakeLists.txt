art_make(LIB_LIBRARIES  lardataobj_Simulation
                        nusimdata_SimulationBase
                        ROOT::Physics
                        ROOT::Core
           MODULE_LIBRARIES larsim_MergeSimSources
                        larsim_Simulation
                        lardataobj_Simulation
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                )

install_headers()
install_fhicl()
install_source()
