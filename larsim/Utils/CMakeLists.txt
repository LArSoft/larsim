art_make(
             LIB_LIBRARIES 
                           larsim_MCCheater_BackTrackerService_service
                           lardataobj_RecoBase
                           art::Persistency_Common canvas::canvas
                           messagefacility::MF_MessageLogger
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
