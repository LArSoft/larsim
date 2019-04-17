include_directories ( $ENV{MuonPropagationHOME}/src )

art_make(MODULE_LIBRARIES
                        larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_Simulation
                        lardataobj_Simulation
                        larcoreobj_SummaryData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        ${GaisserParam}
                        ${MUSUN}
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Provenance
                        ${MF_MESSAGELOGGER}
                        ${CLHEP}
                        ${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()
