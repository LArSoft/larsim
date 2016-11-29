	
art_make( NO_PLUGINS
          
          LIB_LIBRARIES larevt_Filters
                        lardataobj_RawData
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
         )


simple_plugin(SimWireAna "module"
                         lardataobj_RawData
                         ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                         art_Persistency_Common canvas_Persistency_Common
                         art_Persistency_Provenance canvas_Persistency_Provenance
                         art_Utilities canvas_Utilities
                         ${MF_MESSAGELOGGER}
                         ${MF_UTILITIES}
                         ${FHICLCPP}
                         ${CETLIB}
                         ${ROOT_BASIC_LIB_LIST}
              
)

simple_plugin(SimWire "module"
                      larsim_Simulation lardataobj_Simulation
                      lardataobj_RawData
                      larcore_Geometry
                      larcore_Geometry_Geometry_service
                      lardata_Utilities_LArFFT_service
                      lardata_Utilities
                      nutools_RandomUtils_NuRandomService_service
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
              
)

simple_plugin(SimWireT962 "module"
                          larsim_Simulation lardataobj_Simulation
                          lardataobj_RawData
                          larcore_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities_LArFFT_service
                          lardata_Utilities
                          nutools_RandomUtils_NuRandomService_service
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
              
)

simple_plugin(WienerFilterAna "module"
                              lardataobj_RawData
                              larcore_Geometry
                              larcore_Geometry_Geometry_service
                              ${ART_FRAMEWORK_CORE}
                              ${ART_FRAMEWORK_PRINCIPAL}
                              ${ART_FRAMEWORK_SERVICES_REGISTRY}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                              art_Persistency_Common canvas_Persistency_Common
                              art_Persistency_Provenance canvas_Persistency_Provenance
                              art_Utilities canvas_Utilities
                              ${MF_MESSAGELOGGER}
                              ${MF_UTILITIES}
                              ${FHICLCPP}
                              ${CETLIB}
                              ${ROOT_BASIC_LIB_LIST}
              
)

simple_plugin(DumpRawDigits "module"
                            larevt_Filters
                            lardataobj_RawData
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${MF_MESSAGELOGGER}
                            ${FHICLCPP}
                            ${CETLIB}
              
)

simple_plugin(RecompressRawDigits "module"
                lardataobj_RawData
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${MF_MESSAGELOGGER}
                ${FHICLCPP}
                ${CETLIB}
)

install_headers()
install_fhicl()
install_source()
