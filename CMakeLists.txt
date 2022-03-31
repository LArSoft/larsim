# ======================================================================
#  larsim main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsim/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larsim VERSION 09.22.07 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)

cet_report_compiler_flags()

find_package(Range-v3 REQUIRED PUBLIC)
find_package(cetlib_except REQUIRED PUBLIC)
find_package(cetlib REQUIRED PUBLIC)
find_package(fhiclcpp REQUIRED PUBLIC)
find_package(messagefacility REQUIRED PUBLIC)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package(nusimdata REQUIRED PUBLIC)
find_package(nug4 REQUIRED PUBLIC)
find_package(nugen REQUIRED PUBLIC)
find_package(nurandom REQUIRED PUBLIC)
find_package(nutools REQUIRED PUBLIC)
find_package(PostgreSQL REQUIRED PUBLIC)
find_package(CLHEP COMPONENTS Vector REQUIRED PUBLIC)
find_package(ROOT COMPONENTS Core RIO GenVector EG Hist Physics Geom MathCore Matrix Tree RooFit Gpad REQUIRED PUBLIC)
find_package(Boost REQUIRED)
find_package(Geant4 REQUIRED PUBLIC)
find_package(GENIE REQUIRED PUBLIC) # find module in nugen
find_package(log4cpp REQUIRED PUBLIC)
find_package(ifdhc REQUIRED PUBLIC)
find_package(LibXml2 REQUIRED PUBLIC)
find_package(CRY REQUIRED PUBLIC) # find module in nutools

find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( marley )

# Wes put this here to use TRACE for debugging...
#find_ups_product( TRACE )

# marley libraries
cet_find_library( MARLEY NAMES MARLEY PATHS ENV MARLEY_LIB NO_DEFAULT_PATH )
cet_find_library( MARLEY_ROOT NAMES MARLEY_ROOT PATHS ENV MARLEY_LIB NO_DEFAULT_PATH )

# allow the use of #include "Geant4/.."
include_directories("${Geant4_INCLUDE_DIRS}/../")

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(larsim)

# tests
add_subdirectory(test)

# example
add_subdirectory(example)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
