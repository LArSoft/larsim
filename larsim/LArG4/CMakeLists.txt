include_directories($ENV{GENIE_INC}/GENIE)
include_directories($ENV{LIBXML2_FQ_DIR}/include/libxml2)
include_directories($ENV{GEANT4_FQ_DIR}/include)
include_directories($ENV{XERCES_C_INC})

art_make(LIB_LIBRARIES
           larsim_PhotonPropagation_PhotonVisibilityService_service
           lardataobj_Simulation
           larcorealg_Geometry
           nug4_G4Base
           nug4_ParticleNavigation
           nurandom_RandomUtils_NuRandomService_service
           ${ART_UTILITIES}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           canvas
           ${MF_MESSAGELOGGER}
           cetlib_except
           ${ROOT_CORE}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
           ${ROOT_GEOM}
           ${ROOT_MATHCORE}
           ${CLHEP}
           ${G4DIGITS_HITS}
           ${G4EVENT}
           ${G4GEOMETRY}
           ${G4GLOBAL}
           ${G4GRAPHICS_REPS}
           ${G4MATERIALS}
           ${G4INTERCOMS}
           ${G4PARTICLES}
           ${G4PHYSICSLISTS}
           ${G4PROCESSES}
           ${G4TRACK}
           ${G4RUN}
         MODULE_LIBRARIES
           larsim_LArG4
           larsim_MCCheater_ParticleInventoryService_service
           lardataalg_MCDumpers
           larcorealg_Geometry
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ROOT_HIST}
           ${ROOT_CORE}
           ${ROOT_TREE}
         )

install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
