art_make(LIB_LIBRARIES  lardataobj_Simulation
                        nusimdata::SimulationBase
                        ROOT::Physics
                        ROOT::Core
           MODULE_LIBRARIES larsim_MergeSimSources
                        larsim_Simulation
                        lardataobj_Simulation
                        nusimdata::SimulationBase
                        art::Framework_Services_Registry
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                )

install_headers()
install_fhicl()
install_source()
