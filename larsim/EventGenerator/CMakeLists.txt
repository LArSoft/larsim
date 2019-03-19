art_make(LIB_LIBRARIES nusimdata_SimulationBase
                       cetlib
                       cetlib_except
                       ${FHICLCPP}
                       ${CLHEP}
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${ROOT_CORE}
                       ${ROOT_RIO}
                       ${ROOT_HIST}
                       ${ROOT_PHYSICS}
          MODULE_LIBRARIES larsim_PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larsim_Simulation
                           lardataobj_Simulation
                           larcoreobj_SummaryData
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           nutools_RandomUtils_NuRandomService_service
                           larsim_EventGenerator
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Provenance
                           ${MF_MESSAGELOGGER}
                           ${CLHEP}
                           ${ROOT_EG}
                           ${ROOT_GEOM}
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
