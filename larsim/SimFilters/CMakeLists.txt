art_make(MODULE_LIBRARIES
         larcorealg_Geometry
         larsim_MCCheater_ParticleInventoryService_service
         nusimdata_SimulationBase
         art::Framework_Services_Registry
         messagefacility::MF_MessageLogger
         ROOT::Core
         ROOT::Physics)

# install_headers()
install_fhicl()
install_source()
