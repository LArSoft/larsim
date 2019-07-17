art_make(MODULE_LIBRARIES
         larcorealg_Geometry
         larsim_MCCheater_ParticleInventoryService_service
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${MF_MESSAGELOGGER}
         ${ROOT_CORE}
         ${ROOT_PHYSICS})

# install_headers()
install_fhicl()
install_source()
