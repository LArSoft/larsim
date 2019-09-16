simple_plugin(SimWireAna "module"
              lardataobj_RawData
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${ROOT_CORE}
              ${ROOT_HIST})

simple_plugin(SimWire "module"
              lardataobj_RawData
              lardataobj_Simulation
              larcorealg_Geometry
              lardata_Utilities_LArFFT_service
              nurandom_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_ROOT_IO_TFILE_SUPPORT}
              ${CLHEP}
              ${MF_MESSAGELOGGER}
              ${ROOT_CORE}
              ${ROOT_HIST}
              ${ROOT_MATHCORE})

simple_plugin(WienerFilterAna "module"
              larcorealg_Geometry
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${ROOT_CORE}
              ${ROOT_HIST})

install_headers()
install_fhicl()
install_source()
