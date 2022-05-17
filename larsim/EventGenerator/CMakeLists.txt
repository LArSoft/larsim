
art_make( 
          LIB_LIBRARIES
                        nusimdata_SimulationBase
                        art_Persistency_Provenance
                        canvas
                        ${ART_FRAMEWORK_CORE}
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ${FHICLCPP}
                        ${CLHEP}
			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}
                        ${ROOT_GEOM}
          MODULE_LIBRARIES larsim_PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larsim_Simulation lardataobj_Simulation
                           larcoreobj_SummaryData
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           nutools_RandomUtils_NuRandomService_service
                           larsim_EventGenerator
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${CLHEP}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_EG}
                           curl
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
