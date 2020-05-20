include_directories ( ${BxDecay0_INCLUDE_DIR} )

art_make(MODULE_LIBRARIES
         lardataalg_MCDumpers
         larsim_Simulation lardataobj_Simulation
         nusimdata_SimulationBase
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nurandom_RandomUtils_NuRandomService_service
         nugen_EventGeneratorBase_GENIE
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${CLHEP}
         ${BxDecay0_LIBRARIES}
         ${ROOT_EG}
         ${ROOT_HIST}
         ${ROOT_MATHCORE}
        )

install_headers()
install_fhicl()
install_source()
