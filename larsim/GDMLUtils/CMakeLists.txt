include_directories($ENV{GEANT4_FQ_DIR}/include)
include_directories($ENV{XERCES_C_INC})

art_make_exec(NAME simplifyGDML
              LIBRARIES
              ${G4PERSISTENCY}
	      ${G4GEOMETRY}
              ${XERCESC})
