include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LIBXML2_FQ_DIR}/include/libxml2 )
include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{XERCES_C_INC} )

art_make(LIB_LIBRARIES  larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        lardata_Utilities
                        larsim_Simulation  lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib_except
        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GENVECTOR}
                        ${ROOT_EG}
                        ${ROOT_TREEPLAYER}
                        ${ROOT_FFTW}
                        ${ROOT_REFLEX}
                        ${ROOT_EGPYTHIA6}
                        ${ROOT_GUI}
                        ${ROOT_GEOM}
                        ${ROOT_MATHCORE}
                        ${CLHEP}
                        ${G4_LIB_LIST}
           MODULE_LIBRARIES larsim_MergeSimSources
                        larsim_MCCheater_ParticleInventoryService_service
                        lardataalg_MCDumpers
                        larsim_Simulation  lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        ${ROOT_GENVECTOR}
        ${ROOT_BASIC_LIB_LIST}
                        ${CLHEP}
                        ${G4_LIB_LIST}
                )

install_headers()
install_fhicl()
install_source()
