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
                ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
