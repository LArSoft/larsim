art_make(LIB_LIBRARIES
         nurandom::RandomUtils_NuRandomService_service
         canvas::canvas
         cetlib_except::cetlib_except
         ROOT::Core
         ROOT::Matrix
         art::Utilities
         CLHEP::CLHEP) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
