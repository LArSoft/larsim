# for c2: the problem is in fhiclcpp/types/detail/Sequence.icc
cet_add_compiler_flags( CXX -Wno-missing-braces )

art_make(
      LIB_LIBRARIES
                        larsim_Simulation
                        lardataobj_Simulation
                        lardataobj_Simulation_dict
                        nusimdata_SimulationBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        ${MARLEY}
                        ${MARLEY_ROOT}
                        ${FHICLCPP}
                        ${MF_MESSAGELOGGER}
                        cetlib_except
                        ${ROOT_GEOM}
      MODULE_LIBRARIES
                        larsim_EventGenerator_MARLEY
                        larsim_Simulation
                        lardata_Utilities
                        lardataobj_Simulation
                        lardataobj_Simulation_dict
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        ${MARLEY}
                        ${MARLEY_ROOT}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Provenance
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_HIST}
                        ${ROOT_TREE}
                        ${ROOT_RIO}
                        ${G4_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
