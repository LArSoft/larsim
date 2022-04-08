art_make(LIBRARY_NAME     larsim_EventWeight
         LIB_LIBRARIES    larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
         MODULE_LIBRARIES larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          nugen::EventGeneratorBase_GENIE
                          messagefacility::MF_MessageLogger) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
