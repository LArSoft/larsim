include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(LIBRARY_NAME  larsim_EventWeight_Calculators
         LIB_LIBRARIES
                       larsim_EventWeight_Base
		       lardata_Utilities
                       nutools_NuReweight_art   
                       nutools_NuReweight       
		       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
		       nusimdata_SimulationBase
              	       ${ART_FRAMEWORK_CORE}
		       ${ART_FRAMEWORK_PRINCIPAL}
		       ${ART_FRAMEWORK_SERVICES_REGISTRY}
		       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		       art_Persistency_Common canvas canvas
                       art_Persistency_Provenance canvas canvas
                       art_Utilities canvas canvas
                       cetlib cetlib_except
		       ${FHICLCPP}
		       ${CLHEP}
                       ${ROOT_BASIC_LIB_LIST}
                       ${Boost_SYSTEM_LIBRARY}
	) 

install_headers()
install_fhicl()
install_source()

