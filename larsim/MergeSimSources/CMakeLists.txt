art_make(LIB_LIBRARIES  lardataobj_Simulation
                        nusimdata_SimulationBase
                        ROOT::Physics
                        ROOT::Core
           MODULE_LIBRARIES larsim_MergeSimSources
                        larsim_Simulation
                        lardataobj_Simulation
                        nusimdata_SimulationBase
                        art::Framework_Services_Registry
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp fhiclcpp::types
                )

install_headers()
install_fhicl()
install_source()
