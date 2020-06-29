include_directories ( $ENV{GENIE_INC}/GENIE $ENV{PPFX_DIR}/include $ENV{DK2NU_INC})
include_directories ( $ENV{LOG4CPP_INC} )

# Build against extra reweighting functionality in GENIE v3.0.4 MicroBooNE patch 01
add_definitions(-DGENIE_UB_PATCH)

art_make(LIBRARY_NAME  larsim_EventWeight_Calculators
         LIB_LIBRARIES
                       larsim_EventWeight_Base
		       lardata_Utilities
                       nutools_EventGeneratorBase_GENIE
		       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
		       nusimdata_SimulationBase
                       nutools_EventGeneratorBase_GENIE
		       ppfx
                       ${DK2NU_TREE}
                       ${DK2NU_GENIE}
              	       ${ART_FRAMEWORK_CORE}
		       ${ART_FRAMEWORK_PRINCIPAL}
		       ${ART_FRAMEWORK_SERVICES_REGISTRY}
		       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		       art_Persistency_Common canvas canvas
                       art_Persistency_Provenance canvas canvas
                       art_Utilities canvas canvas
                       cetlib cetlib_except
                       ${GENIE_LIB_LIST}
		       ${FHICLCPP}
		       ${CLHEP}
		       ${LOG4CPP}
                       ${ROOT_BASIC_LIB_LIST}
                       ${Boost_SYSTEM_LIBRARY}
	) 

install_headers()
install_fhicl()
install_source()

