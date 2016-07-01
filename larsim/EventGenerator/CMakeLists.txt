
art_make( 
          LIB_LIBRARIES
                        nusimdata_SimulationBase
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        ${ART_FRAMEWORK_CORE}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${FHICLCPP}
                        ${CLHEP}
			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES larsim_PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larsim_Simulation larsimobj_Simulation
                           larcoreobj_SummaryData
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_RandomUtils_LArSeedService_service
                           larsim_EventGenerator
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           art_Persistency_Common canvas_Persistency_Common
                           art_Persistency_Provenance canvas_Persistency_Provenance
                           art_Utilities canvas_Utilities
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${CLHEP}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_EG}
        )

simple_plugin( GeneratedEventTimestamp "plugin"
                                       ${ART_FRAMEWORK_CORE}
                                       art_Persistency_Provenance canvas_Persistency_Provenance
                                       ${MF_MESSAGELOGGER}
                                       ${MF_UTILITIES}
                                       ${CETLIB}
				       ${Boost_FILESYSTEM_LIBRARY}
				       ${Boost_SYSTEM_LIBRARY}
               
             )

install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(MuonPropagation) 
add_subdirectory(CORSIKA)

# See ubooneoffline repository for a functioning NuWro implementation. 
# Leave it commented out here, as below.
##add_subdirectory(NuWro)
