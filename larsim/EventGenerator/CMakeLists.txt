art_make(LIB_LIBRARIES
           nusimdata_SimulationBase
           cetlib
           cetlib_except
           ${FHICLCPP}
           ${CLHEP}
           ${ROOT_CORE}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
           ${ROOT_RIO}
         MODULE_LIBRARIES
           larcorealg_Geometry
           larcoreobj_SummaryData
           larsim_EventGenerator
           larsim_PhotonPropagation_PhotonVisibilityService_service
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${CLHEP}
           ${MF_MESSAGELOGGER}
           ${ROOT_CORE}
           ${ROOT_EG}
           ${ROOT_GEOM}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
           ${ROOT_RIO}
        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(MARLEY)
add_subdirectory(MuonPropagation)
add_subdirectory(CORSIKA)

# See ubooneoffline repository for a functioning NuWro implementation.
# Leave it commented out here, as below.
##add_subdirectory(NuWro)
