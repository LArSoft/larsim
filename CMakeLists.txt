# ======================================================================
#  larsim main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsim/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19) ### MIGRATE-ACTION-COMPATIBILITY (migrate-3.10.00) - Uncomment to preserve compatibility with older CMake versions

find_package(cetmodules REQUIRED)
project(larsim VERSION 09.22.10 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - add args: REPORT_THRESHOLD VERBOSE

find_package(Range-v3 REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED EXPORT)
find_package(art_root_io REQUIRED EXPORT)
find_package(nusimdata REQUIRED EXPORT)
find_package(nug4 REQUIRED EXPORT)
find_package(nugen REQUIRED EXPORT)
find_package(nurandom REQUIRED EXPORT)
find_package(nutools REQUIRED EXPORT)
find_package(PostgreSQL REQUIRED EXPORT)
find_package(CLHEP COMPONENTS Vector REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core RIO GenVector EG Hist Physics Geom MathCore Matrix Tree RooFit Gpad REQUIRED EXPORT)
find_package(Boost REQUIRED)
find_package(Geant4 REQUIRED EXPORT)
find_package(GENIE REQUIRED EXPORT) # find module in nugen
find_package(log4cpp REQUIRED EXPORT)
find_package(ifdhc REQUIRED EXPORT)
find_package(LibXml2 REQUIRED EXPORT)
find_package(CRY REQUIRED EXPORT) # find module in nutools

find_package( larcoreobj REQUIRED )
find_package( larcorealg REQUIRED )
find_package( larcore REQUIRED )
find_package( lardataobj REQUIRED )
find_package( lardataalg REQUIRED )
find_package( lardata REQUIRED )
find_package( larevt REQUIRED )
find_package( marley REQUIRED )

# Wes put this here to use TRACE for debugging...
#find_ups_product( TRACE )

# marley libraries
cet_find_library( MARLEY NAMES MARLEY PATHS ENV MARLEY_LIB NO_DEFAULT_PATH ) ### MIGRATE-ACTION-REQUIRED (migrate-3.10.00) - avoid: use find_package() with custom Find<pkg>.cmake for Spack compatibility
cet_find_library( MARLEY_ROOT NAMES MARLEY_ROOT PATHS ENV MARLEY_LIB NO_DEFAULT_PATH ) ### MIGRATE-ACTION-REQUIRED (migrate-3.10.00) - avoid: use find_package() with custom Find<pkg>.cmake for Spack compatibility

# allow the use of #include "Geant4/.."
include_directories("${Geant4_INCLUDE_DIRS}/../") ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

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

# packaging utility
cet_cmake_config()
