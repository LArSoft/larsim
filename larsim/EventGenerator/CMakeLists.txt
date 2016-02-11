
art_make( 
          LIB_LIBRARIES
                        ${SIMULATIONBASE}
                        ${ART_PERSISTENCY_PROVENANCE}
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
                           larsim_Simulation
                           larcore_SummaryData
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           SeedService_service # artextensions
                           larsim_EventGenerator
                           ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_PERSISTENCY_PROVENANCE}
                           ${ART_UTILITIES}
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
                                       ${ART_PERSISTENCY_PROVENANCE}
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
