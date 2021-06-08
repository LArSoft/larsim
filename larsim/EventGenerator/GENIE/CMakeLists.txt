include_directories ( ${GENIE_INC}/GENIE )

art_make(MODULE_LIBRARIES
         lardataalg_MCDumpers
         larsim_Simulation lardataobj_Simulation
         nusimdata_SimulationBase
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nurandom_RandomUtils_NuRandomService_service
         nugen_EventGeneratorBase_GENIE
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         messagefacility::MFMessageLogger
         ${CLHEP_LIBRARIES}
         ${GENIE_LIB_LIST}
         ROOT::EGPythia6    # FIXME!!! - resolving genie run time reference
         ROOT::EG
         ROOT::Hist
         ROOT::MathCore
        )

install_headers()
install_fhicl()
install_source()
