
#include_directories ( $ENV{MARLEY_INC} )

art_make( 
          MODULE_LIBRARIES 
                        larsim_Simulation lardataobj_Simulation
                        lardataobj_Simulation
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        ${MARLEY}
                        ${MARLEY_ROOT}
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
                        ${ROOT_GEOM}
                        ${G4_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
