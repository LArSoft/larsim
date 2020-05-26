art_make(MODULE_LIBRARIES
         lardataalg_MCDumpers
         nusimdata_SimulationBase
         ${MF_MESSAGELOGGER}
         ROOT::GenVector
         ROOT::Core
         ROOT::Physics)

install_headers()
install_fhicl()
install_source()
