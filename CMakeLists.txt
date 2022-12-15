# ======================================================================
#  larsim main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsim/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.11.02 REQUIRED)
project(larsim VERSION 09.30.01 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

cet_cmake_module_directories(Modules BINARY)

find_package(art REQUIRED EXPORT)
find_package(art_root_io REQUIRED EXPORT)
find_package(canvas_root_io REQUIRED)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)

find_package(nug4 REQUIRED EXPORT)
find_package(nugen REQUIRED EXPORT)
find_package(nurandom REQUIRED EXPORT)
find_package(nusimdata REQUIRED EXPORT)
find_package(nutools REQUIRED EXPORT)

find_package(Boost REQUIRED COMPONENTS headers math_tr1 REQUIRED EXPORT) # boost::math::policies, boost::math::ellint*
find_package(CLHEP COMPONENTS Evaluator Geometry Random Vector REQUIRED EXPORT)
find_package(CRY REQUIRED EXPORT) # find module in nutools
find_package(GENIE REQUIRED EXPORT) # find module in nugen
find_package(Geant4 REQUIRED EXPORT)
find_package(LibXml2 REQUIRED EXPORT)
find_package(MARLEY REQUIRED EXPORT)
find_package(PostgreSQL REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core EG GenVector Geom Gpad Hist MathCore Matrix Physics RIO RooFitCore Tree REQUIRED EXPORT)
find_package(Range-v3 REQUIRED EXPORT)
find_package(SQLite3 REQUIRED EXPORT)
find_package(ifdhc REQUIRED EXPORT)
find_package(log4cpp REQUIRED EXPORT)

find_package(larcore REQUIRED EXPORT)
find_package(larcorealg REQUIRED  EXPORT)
find_package(larcoreobj REQUIRED EXPORT)
find_package(lardata REQUIRED EXPORT)
find_package(lardataalg REQUIRED EXPORT)
find_package(lardataobj REQUIRED EXPORT)
find_package(larevt REQUIRED EXPORT)

# Wes put this here to use TRACE for debugging...
#find_ups_product( TRACE )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BuildPlugins)

# source
add_subdirectory(larsim)

add_subdirectory(Modules)

# tests
add_subdirectory(test)

# example
add_subdirectory(example)

# ups - table and config files

# packaging utility
cet_cmake_config()
