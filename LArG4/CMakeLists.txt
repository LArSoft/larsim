
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LIBXML2_FQ_DIR}/include/libxml2 )
include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{XERCES_C_INC} )
include_directories ( $ENV{CRYHOME}/src )

art_make( BASENAME_ONLY
          LIBRARY_NAME  LArG4
          LIB_LIBRARIES 
                        PhotonPropagation
                        PhotonVisibilityService_service
                        Simulation
                        Utilities
                        LArProperties_service
                        DetectorProperties_service
	                Geometry
	                Geometry_service
                        ${SIMULATIONBASE}
                        ${G4BASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_EG}
			${ROOT_TREEPLAYER} 
			${ROOT_FFTW}
			${ROOT_REFLEX}
			${ROOTSYS}/lib/libEGPythia6.so
			${ROOT_GUI}
			${ROOT_BASIC_LIB_LIST}
	                ${CLHEP}
                        ${G4_LIB_LIST}
           MODULE_LIBRARIES LArG4
                        Simulation
	                Geometry
	                Geometry_service
                        ${SIMULATIONBASE}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
      G4Base      
      CLHEP
      ${G4_LIB_LIST}
		)

install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install( FILES ${g4_macros} DESTINATION ${${product}_lib_dir} )
# copy into build for tests
# mrb will add g4macropath to FW_SEARCH_PATH
set( mrbdir $ENV{MRB_BUILDDIR} )
if( mrbdir )
  set( g4macropath ${mrbdir}/${product}/lib )
else()
  set( g4macropath ${CETPKG_BUILD}/lib )
endif()
foreach( g4m ${g4_macros} )
  get_filename_component( g4nm ${g4m} NAME )
  configure_file( ${g4m} ${g4macropath}/${g4nm} COPYONLY )
endforeach( g4m )
