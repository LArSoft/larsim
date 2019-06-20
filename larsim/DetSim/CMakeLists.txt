simple_plugin(SimWireAna "module"
                         lardataobj_RawData
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                         art_Persistency_Provenance
                         ${MF_MESSAGELOGGER}
                         ${ROOT_HIST})

simple_plugin(SimWire "module"
                      larsim_Simulation lardataobj_Simulation
                      lardataobj_RawData
                      larcorealg_Geometry
                      larcore_Geometry_Geometry_service
                      lardata_Utilities_LArFFT_service
                      lardata_Utilities
                      nurandom_RandomUtils_NuRandomService_service
                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
                      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                      art_Persistency_Provenance
                      ${MF_MESSAGELOGGER}
                      ${CLHEP}
                      ${ROOT_MATHCORE}
                      ${ROOT_HIST})


simple_plugin(WienerFilterAna "module"
                              lardataobj_RawData
                              larcorealg_Geometry
                              larcore_Geometry_Geometry_service
                              ${ART_FRAMEWORK_SERVICES_REGISTRY}
                              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                              art_Persistency_Provenance
                              ${MF_MESSAGELOGGER}
                              ${ROOT_HIST})

install_headers()
install_fhicl()
install_source()
