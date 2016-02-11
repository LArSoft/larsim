	
art_make( NO_PLUGINS
          
          LIB_LIBRARIES larevt_Filters
                        lardata_RawData
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
         )


simple_plugin(SimWireAna "module"
                         lardata_RawData
                         ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                         ${ART_PERSISTENCY_COMMON}
                         ${ART_PERSISTENCY_PROVENANCE}
                         ${ART_UTILITIES}
                         ${MF_MESSAGELOGGER}
                         ${MF_UTILITIES}
                         ${FHICLCPP}
                         ${CETLIB}
                         ${ROOT_BASIC_LIB_LIST}
              
)

simple_plugin(SimWire "module"
                      larsim_Simulation
                      lardata_RawData
                      larcore_Geometry
                      larcore_Geometry_Geometry_service
                      lardata_Utilities_LArFFT_service
                      lardata_Utilities_LArProperties_service
                      lardata_Utilities_DetectorProperties_service
                      lardata_Utilities
                      SeedService_service # artextensions
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
              
)

simple_plugin(SimWireT962 "module"
                          larsim_Simulation
                          lardata_RawData
                          larcore_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities_LArFFT_service
                          lardata_Utilities_LArProperties_service
                          lardata_Utilities_DetectorProperties_service
                          lardata_Utilities
                          SeedService_service # artextensions
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
              
)

simple_plugin(WienerFilterAna "module"
                              lardata_RawData
                              larcore_Geometry
                              larcore_Geometry_Geometry_service
                              ${ART_FRAMEWORK_CORE}
                              ${ART_FRAMEWORK_PRINCIPAL}
                              ${ART_FRAMEWORK_SERVICES_REGISTRY}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                              ${ART_PERSISTENCY_COMMON}
                              ${ART_PERSISTENCY_PROVENANCE}
                              ${ART_UTILITIES}
                              ${MF_MESSAGELOGGER}
                              ${MF_UTILITIES}
                              ${FHICLCPP}
                              ${CETLIB}
                              ${ROOT_BASIC_LIB_LIST}
              
)

simple_plugin(DumpRawDigits "module"
                            larevt_Filters
                            lardata_RawData
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${MF_MESSAGELOGGER}
                            ${FHICLCPP}
                            ${CETLIB}
              
)

install_headers()
install_fhicl()
install_source()
