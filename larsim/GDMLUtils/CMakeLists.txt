include_directories($ENV{GEANT4_FQ_DIR}/include)
include_directories($ENV{XERCES_C_INC})

art_make_exec(NAME simplifyGDML
              LIBRARIES
              Geant4::G4global
              Geant4::G4persistency
              Geant4::G4geometry
              ${XERCESC})
