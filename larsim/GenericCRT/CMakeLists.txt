#include_directories ( $ENV{GENIE_INC}/GENIE )
#include_directories ( $ENV{LIBXML2_FQ_DIR}/include/libxml2 )
#include_directories ( $ENV{GEANT4_FQ_DIR}/include )
#include_directories ( $ENV{XERCES_C_INC} )

art_make(LIB_LIBRARIES  larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        lardata_Utilities
                        larsim_Simulation  lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata::SimulationBase
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib_except::cetlib_except
                        ROOT::GenVector
                        ROOT::EG
                  #      ROOT::TreePlayer
                  #      ROOT::FFTW
                        
                  #      ROOT::EGPythia6
                  #      ROOT::Gui
                        ROOT::Geom
                        ROOT::MathCore
                        CLHEP::CLHEP
                        ${G4_LIB_LIST}
           MODULE_LIBRARIES larsim_GenericCRT 
                        larsim_MCCheater_ParticleInventoryService_service
                        lardataalg_MCDumpers
                        larsim_Simulation  lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata::SimulationBase
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        messagefacility::MF_MessageLogger
                        ROOT::GenVector
                )

install_headers()
install_fhicl()
install_source()
