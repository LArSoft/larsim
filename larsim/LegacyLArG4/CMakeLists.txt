include_directories($ENV{GENIE_INC}/GENIE)
include_directories($ENV{LIBXML2_FQ_DIR}/include/libxml2)
include_directories($ENV{GEANT4_FQ_DIR}/include)
include_directories($ENV{XERCES_C_INC})

art_make(LIB_LIBRARIES
           larsim_PhotonPropagation_PhotonVisibilityService_service
	   larsim_Utils
           lardataobj_Simulation
           larcorealg_Geometry
           nug4_G4Base
           nug4_ParticleNavigation
           nurandom_RandomUtils_NuRandomService_service
           art::Utilities
           art::Framework_Services_Registry
           art_root_io::tfile_support
           art_root_io::TFileService_service
           canvas
           messagefacility::MF_MessageLogger
           cetlib_except
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ROOT::Geom
           ROOT::MathCore
           ${CLHEP_LIBRARIES}
           G4digits_hits
           G4event
           G4geometry
           G4global
           G4graphics_reps
           G4materials
           G4intercoms
           G4particles
           G4physicslists
           G4processes
           G4track
           G4run

         MODULE_LIBRARIES
           larsim_LegacyLArG4
           larsim_MCCheater_ParticleInventoryService_service
           lardataalg_MCDumpers
           larcorealg_Geometry
           nusimdata_SimulationBase
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
