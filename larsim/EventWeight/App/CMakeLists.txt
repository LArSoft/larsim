art_make(LIBRARY_NAME     larsim_EventWeight
	 LIB_LIBRARIES
                          larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          ${Boost_SYSTEM_LIBRARY}
         MODULE_LIBRARIES larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
			  ${MF_MESSAGELOGGER}
			  ${GENIE_LIB_LIST}
                     	  
	)


install_headers()
install_fhicl()
install_source()
