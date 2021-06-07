art_make(LIB_LIBRARIES
         canvas
         cetlib_except
         ${CLHEP_LIBRARIES}
         ROOT::Core
         ROOT::Matrix
         ${ART_UTILITIES})

install_headers()
install_fhicl()
install_source()
