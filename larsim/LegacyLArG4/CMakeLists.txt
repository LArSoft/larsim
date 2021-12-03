include_directories($ENV{GENIE_INC}/GENIE)
include_directories($ENV{LIBXML2_FQ_DIR}/include/libxml2)
#include_directories($ENV{GEANT4_FQ_DIR}/include)
#include_directories($ENV{XERCES_C_INC})

art_make(LIB_LIBRARIES
           larsim_PhotonPropagation_PhotonVisibilityService_service
	   larsim_Utils
           lardataobj_Simulation
           larcorealg_Geometry
           nug4::G4Base
           nug4::ParticleNavigation
           nug4::AdditionalG4Physics
           nurandom::RandomUtils_NuRandomService_service
           art::Utilities
           art::Framework_Services_Registry
           art_root_io::tfile_support
           art_root_io::TFileService_service
           canvas::canvas
           messagefacility::MF_MessageLogger
           cetlib_except::cetlib_except
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ROOT::Geom
           ROOT::MathCore
           Geant4::G4digits_hits
           Geant4::G4event
           Geant4::G4geometry
           Geant4::G4global
           Geant4::G4graphics_reps
           Geant4::G4materials
           Geant4::G4intercoms
           Geant4::G4particles
           Geant4::G4physicslists
           Geant4::G4processes
           Geant4::G4track
           Geant4::G4run
         MODULE_LIBRARIES
           larsim_LegacyLArG4
           larsim_MCCheater_ParticleInventoryService_service
           lardataalg_MCDumpers
           larcorealg_Geometry
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           ROOT::Hist
           ROOT::Core
           ROOT::Tree
         )

install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
