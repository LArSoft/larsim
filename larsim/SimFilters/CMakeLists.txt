art_make(MODULE_LIBRARIES
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         larsim_Simulation
         nug4_ParticleNavigation
         lardataobj_Simulation
         larsim_MCCheater_BackTrackerService_service
         larsim_MCCheater_ParticleInventoryService_service
         lardata_Utilities
         larevt_Filters
         lardataobj_RawData
         lardataobj_RecoBase
         lardata_RecoObjects
         nusimdata_SimulationBase
         ${LARRECO_LIB}
         ${LARDATA_LIB}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${ROOT_GEOM}
         ${ROOT_XMLIO}
         ${ROOT_GDML}
         ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()
