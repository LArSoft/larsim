add_subdirectory(job)

art_make(LIB_LIBRARIES larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       larsim_Simulation
                       nug4_ParticleNavigation
                       lardataobj_Simulation
                       nusimdata_SimulationBase
                       larcoreobj_SummaryData
                       lardata_Utilities
                       lardataobj_MCBase
                       canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${ROOT_GEOM}
                       ${ROOT_XMLIO}
                       ${ROOT_GDML}
          MODULE_LIBRARIES larsim_MCSTReco
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation
                           nug4_ParticleNavigation
                           lardataobj_Simulation
                           nusimdata_SimulationBase
                           larcoreobj_SummaryData
                           lardata_Utilities
                           lardataobj_MCBase
                           nurandom_RandomUtils_NuRandomService_service
                           art_Persistency_Provenance
                           ${MF_MESSAGELOGGER}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_EG}
                           ${ROOT_HIST}
                           ${CLHEP}
                        )

install_headers()
install_fhicl()
install_source()
