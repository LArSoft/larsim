art_make(LIBRARY_NAME     larsim_EventWeight
         LIB_LIBRARIES    larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
         MODULE_LIBRARIES larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          nugen_EventGeneratorBase_GENIE
                          ${MF_MESSAGELOGGER})

install_headers()
install_fhicl()
install_source()
