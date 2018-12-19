# ======================================================================
#  larsim main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsim/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

project(larsim)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( marley )
find_ups_product( nutools )
find_ups_product( nusimdata )
find_ups_product( art )
find_ups_product( cetbuildtools )
find_ups_geant4( geant4 )
find_ups_product( postgresql )
find_ups_product( ifdhc )
find_ups_product( genie )

# Wes put this here to use TRACE for debugging...
#find_ups_product( TRACE )

#  Find all the libraries needed by our dependent CMakeList.txt files
cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
#
find_library( CRY NAMES CRY PATHS $ENV{CRYHOME}/lib NO_DEFAULT_PATH )
# genie

# marley libraries
cet_find_library( MARLEY NAMES MARLEY PATHS ENV MARLEY_LIB NO_DEFAULT_PATH )
cet_find_library( MARLEY_ROOT NAMES MARLEY_ROOT PATHS ENV MARLEY_LIB NO_DEFAULT_PATH )

# temporarily needed since this is an unexpected header path
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# Determine the major version number of GENIE based on the GENIE_VERSION
# variable set by the genie ups product. Define the GENIE_PRE_R3 preprocessor
# macro (needed to set the right GENIE header locations, etc.) if needed.
#
# The top-level nutools CMakeLists.txt does this as well, but I don't think
# that we can access the result of that version check in the larsim build.
parse_ups_version( ${GENIE_VERSION} )
if( ${VMAJ} LESS 3 )

  add_definitions(-DGENIE_PRE_R3)

  # GENIE v2 library list
  set(GENIE_LIB_LIST  ${GALGORITHM}
                    ${GBARYONRESONANCE}
                    ${GBASE}
                    ${GBODEKYANG}
                    ${GCHARM}
                    ${GCOH}
                    ${GDFRC}
                    ${GDIS}
                    ${GCROSSSECTIONS}
                    ${GDECAY}
                    ${GELAS}
                    ${GELFF}
                    ${GHEP}
                    ${GEVGCORE}
                    ${GEVGMODULES}
                    ${GEVGDRIVERS}
                    ${GGIBUU}
                    ${GHADRONTRANSP}
                    ${GFRAGMENTATION}
                    ${GINTERACTION}
                    ${GLLEWELLYNSMITH}
                    ${GMEC}
                    ${GMESSENGER}
                    ${GNUGAMMA}
                    ${GNUE}
                    ${GNTUPLE}
                    ${GNUCLEAR}
                    ${GNUMERICAL}
                    ${GQPM}
                    ${GPDG}
                    ${GPDF}
                    ${GQEL}
                    ${GRES}
                    ${GREGISTRY}
                    ${GREINSEHGAL}
                    ${GUTILS}
                    ${GGEO}
                    ${GFLUXDRIVERS}
                    ${GMUELOSS}
                    ${GREWEIGHT}
                    ${GNUCLEONDECAY} )

  #MESSAGE("--LARSIM-- GENIE_VERSION=${VMAJ} ${VMIN} -DGENIE_PRE_R3")
else()

  # GENIE v3 library list
  set (GENIE_LIB_LIST ${GFWMSG}
                      ${GFWREG}
                      ${GFWALG}
                      ${GFWINT}
                      ${GFWGHEP}
                      ${GFWNUM}
                      ${GSL}  # FWNUM relies on GSL
                      ${GFWUTL}
                      ${GFWPARDAT}
                      ${GFWEG}
                      ${GFWNTP}
                      ${GPHXSIG}
                      ${GPHPDF}
                      ${GPHNUCLST}
                      ${GPHCMN}
                      ${GPHDCY}
                      ${GPHHADTRANSP}
                      ${GPHHADNZ}
                      ${GPHDEEX}
                      ${GPHAMNGXS}
                      ${GPHAMNGEG}
                      ${GPHCHMXS}
                      ${GPHCOHXS}
                      ${GPHCOHEG}
                      ${GPHDISXS}
                      ${GPHDISEG}
                      ${GPHDFRCXS}
                      ${GPHDFRCEG}
                      ${GPHGLWRESXS}
                      ${GPHGLWRESEG}
                      ${GPHIBDXS}
                      ${GPHIBDEG}
                      ${GPHMNUCXS}
                      ${GPHMNUCEG}
                      ${GPHMEL}
                      ${GPHNUELXS}
                      ${GPHNUELEG}
                      ${GPHQELXS}
                      ${GPHQELEG}
                      ${GPHRESXS}
                      ${GPHRESEG}
                      ${GPHSTRXS}
                      ${GPHSTREG}
                      ${GPHNDCY}
                      ${GTLGEO}
                      ${GTLFLX}
                      ${GRWFWK}
                      ${GRWIO}
                      ${GRWCLC})

  #MESSAGE("--LARSIM-- GENIE_VERSION=${VMAJ} ${VMIN} not defining -DGENIE_PRE_R3")
endif()

#MESSAGE("--LARSIM-- GENIE ${VMAJ} GENIE_LIB_LIST=${GENIE_LIB_LIST}")

# source
add_subdirectory(larsim)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
