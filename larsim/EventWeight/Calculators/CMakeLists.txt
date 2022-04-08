
# Builds of MicroBooNE-specific larsim releases should uncomment the
# "add_definitions" line below. It enables extra reweighting functionality
# which is present only in "MicroBooNE patch" releases of GENIE v3.0.4. These
# require building against a special ups product.
#add_definitions(-DGENIE_UB_PATCH)

art_make(LIB_LIBRARIES
         larsim_EventWeight_Base
         nurandom::RandomUtils_NuRandomService_service
         nugen::EventGeneratorBase_GENIE
         nusimdata::SimulationBase
         art::Framework_Principal
         art::Persistency_Provenance
         art::Framework_Services_Registry
         fhiclcpp::fhiclcpp
         ${GENIE_LIB_LIST}
         log4cpp::log4cpp
         messagefacility::MF_MessageLogger
         cetlib_except::cetlib_except) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
