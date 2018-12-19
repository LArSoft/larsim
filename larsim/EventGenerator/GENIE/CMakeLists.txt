
cet_add_compiler_flags( CXX -DHIDE_GENIE_LOG_XXX )

include_directories ( $ENV{GENIE_INC}/GENIE $ENV{LOG4CPP_INC} )

art_make( 
          MODULE_LIBRARIES 
                        lardataalg_MCDumpers
                        larsim_Simulation lardataobj_Simulation
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        nutools_EventGeneratorBase_GENIE
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
                        ${GENIE_LIB_LIST}
                        ${ROOT_EGPYTHIA6}    # FIXME!!! - resolving genie run time reference
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
                        ${LOG4CPP}
        )

install_headers()
install_fhicl()
install_source()

