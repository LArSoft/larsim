
cet_add_compiler_flags( CXX -DHIDE_GENIE_LOG_XXX )

include_directories ( $ENV{LOG4CPP_INC} )
include_directories ( $ENV{DK2NUDATA_INC} )
include_directories ( $ENV{DK2NUGENIE_INC} )

# genie
cet_find_library( DK2NUDATA  NAMES dk2nuTree  PATHS ENV DK2NUDATA_LIB )
cet_find_library( DK2NUGENIE NAMES dk2nuGenie PATHS ENV DK2NUGENIE_LIB )

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
			${DK2NUDATA}
                        ${DK2NUGENIE}
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

