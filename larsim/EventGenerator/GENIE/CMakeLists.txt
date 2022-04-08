
art_make(MODULE_LIBRARIES
         lardataalg_MCDumpers
         larsim_Simulation lardataobj_Simulation
         nusimdata::SimulationBase
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nurandom::RandomUtils_NuRandomService_service
         nugen::EventGeneratorBase_GENIE
         art::Framework_Services_Registry
         art_root_io::tfile_support ROOT::Core
         art_root_io::TFileService_service
         art::Persistency_Provenance
         messagefacility::MF_MessageLogger
         ${GENIE_LIB_LIST}
         log4cpp::log4cpp
         ROOT::EGPythia6    # FIXME!!! - resolving genie run time reference
         ROOT::EG
         ROOT::Hist
         ROOT::MathCore
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
