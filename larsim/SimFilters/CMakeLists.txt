art_make(MODULE_LIBRARIES
         larcorealg_Geometry
         larsim_MCCheater_ParticleInventoryService_service
         nusimdata::SimulationBase
         art::Framework_Services_Registry
         messagefacility::MF_MessageLogger
         ROOT::Core
         ROOT::Physics) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

# install_headers()
install_fhicl()
install_source()
