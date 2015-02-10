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

set(CMAKE_MODULE_PATH $ENV{ART_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcore v1_00_00 )
find_ups_product( lardata v1_00_00 )
find_ups_product( larevt v1_00_00 )
find_ups_product( nutools v1_00_00 )
find_ups_product( art v1_09_00 )
find_ups_product( cetbuildtools v3_10_00 )
find_ups_geant4( geant4 v4_9_6 )
find_ups_product( postgresql v9_1_5 )

#  Find all the libraries needed by our dependent CMakeList.txt files
cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
#
find_library( CRY NAMES CRY PATHS $ENV{CRYHOME}/lib NO_DEFAULT_PATH )
# genie
cet_find_library( GALGORITHM NAMES GAlgorithm PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GBARYONRESONANCE NAMES GBaryonResonance PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GBASE NAMES GBase PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GBODEKYANG NAMES GBodekYang PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GCHARM NAMES GCharm PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GCOH NAMES GCoh PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GDFRC NAMES GDfrc PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GDIS NAMES GDIS PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GCROSSSECTIONS NAMES GCrossSections PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GDECAY NAMES GDecay PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GELAS NAMES GElas PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GELFF NAMES GElFF PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GHEP NAMES GHEP PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GEVGCORE NAMES GEVGCore  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GEVGMODULES NAMES GEVGModules PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GEVGDRIVERS NAMES GEVGDrivers PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GGIBUU NAMES GGiBUU PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GHADRONTRANSP NAMES GHadronTransp PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFRAGMENTATION NAMES GFragmentation PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GINTERACTION NAMES GInteraction PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GLLEWELLYNSMITH NAMES GLlewellynSmith  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GMEC NAMES GMEC PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GMESSENGER NAMES GMessenger PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUGAMMA NAMES GNuGamma PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUE NAMES GNuE PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNTUPLE NAMES GNtuple PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUCLEAR NAMES GNuclear PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUMERICAL NAMES GNumerical PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GQPM NAMES GQPM PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GPDG NAMES GPDG PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GPDF NAMES GPDF PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GQEL NAMES GQEL PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GRES NAMES GRES PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GREGISTRY NAMES GRegistry PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GREINSEGHAL NAMES GReinSeghal PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GUTILS NAMES GUtils PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GGEO NAMES GGeo PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFLUXDRIVERS NAMES GFluxDrivers PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GMUELOSS NAMES GMuELoss PATHS ENV GENIE_LIB NO_DEFAULT_PATH )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(DetSim)
add_subdirectory(EventGenerator)
add_subdirectory(LArG4)
add_subdirectory(MCCheater)
add_subdirectory(PhotonPropagation)
add_subdirectory(Simulation)
add_subdirectory(TriggerAlgo)
add_subdirectory(SimFilters)
add_subdirectory(MCSTReco)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
