#[============================================================[.rst:
FindMARLEY
============
#]============================================================]
find_program(MARLEY_CONFIG_EXECUTABLE NAMES marley-config HINTS ENV MARLEY_FQ_DIR)
mark_as_advanced(MARLEY_CONFIG_EXECUTABLE)

set(_fm_libs MARLEY)
set(_fm_progs marley)
set(_fm_transitive_deps_MARLEY)
set(_fm_transitive_deps_MARLEY_ROOT MARLEY::MARLEY ROOT::Tree ROOT::Hist ROOT::RIO)
set(_fm_deps GSL)
set(_fm_fp_ROOT_args COMPONENTS Core Tree Hist RIO EXPORT)

unset(_fm_fphsa_extra_required_vars)
unset(_fm_fphsa_extra_args)
if (MARLEY_CONFIG_EXECUTABLE)
  if (MARLEY_FOUND)
    if (MARLEY_ROOT)
      list(PREPEND _fm_libs MARLEY_ROOT)
      list(PREPEND _fm_deps ROOT)
      list(APPEND _fm_progs marsum mroot)
    endif()
  else()
    execute_process(COMMAND ${MARLEY_CONFIG_EXECUTABLE}
      --use-root --cxx-std --version --incdir
      OUTPUT_VARIABLE _fm_config_output
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    string(REGEX MATCH
      "^(yes|no)[ \\t]+c\\+\\+([0-9]+)[ \\t]+([^ \\t]+)[ \\t]+(.*)$"
      _fm_config_output
      "${_fm_config_output}")
    unset(_fm_config_output)
    if (CMAKE_MATCH_1 STREQUAL "yes")
      list(PREPEND _fm_libs MARLEY_ROOT)
      list(PREPEND _fm_deps ROOT)
      list(APPEND _fm_progs marsum mroot)
      set(MARLEY_DEFINITIONS "USE_ROOT")
    else()
      unset(MARLEY_DEFINITIONS)
    endif()
    set(MARLEY_CXX_STANDARD "${CMAKE_MATCH_2}")
    set(MARLEY_VERSION "${CMAKE_MATCH_3}")
    set(MARLEY_INCLUDE_DIR "${CMAKE_MATCH_4}")
    set(MARLEY_LIBRARIES)
    foreach (_fm_lib IN LISTS _fm_libs)
      find_library(${_fm_lib}_LIBRARY NAMES ${_fm_lib} HINTS ENV MARLEY_LIB)
      mark_as_advanced(${_fm_lib}_LIBRARY)
      if (${_fm_lib}_LIBRARY)
        list(APPEND MARLEY_LIBRARIES "${${_fm_lib}_LIBRARY}")
      endif()
    endforeach()
  endif()
else()
  set(_fm_fphsa_extra_args "could not find executable \"marley-config\"")
endif()

list(TRANSFORM _fm_deps APPEND _FOUND
  OUTPUT_VARIABLE _fm_fphsa_extra_required_vars)

if (MARLEY_FOUND OR (MARLEY_LIBRARIES AND MARLEY_INCLUDE_DIR))
  foreach (_fm_prog IN LISTS _fm_progs)
    string(MAKE_C_IDENTIFIER "${_fm_prog}" _fm_prog_var)
    string(TOUPPER "${_fm_prog_var}_EXECUTABLE" _fm_prog_var)
    find_program(${_fm_prog_var} NAMES ${_fm_prog})
    mark_as_advanced(${_fm_prog_var})
    list(APPEND _fm_fphsa_extra_required_vars ${_fm_prog_var})
  endforeach()
  unset(_fm_prog_var)
  unset(_fm_missing_deps)
  foreach (_fm_dep IN LISTS _fm_deps)
    list(APPEND _fm_fphsa_extra_required_vars "${_fm_dep}_FOUND")
    get_property(_fm_${_fm_dep}_alreadyTransitive GLOBAL PROPERTY
      _CMAKE_${_fm_dep}_TRANSITIVE_DEPENDENCY)
    find_package(${_fm_dep} ${_fm_fp_${_fm_dep}_args} QUIET)
    if (NOT DEFINED cet_${_fm_dep}_alreadyTransitive OR cet_${_fm_dep}_alreadyTransitive)
      set_property(GLOBAL PROPERTY _CMAKE_${_fm_dep}_TRANSITIVE_DEPENDENCY TRUE)
    endif()
    unset(_fm_${_fm_dep}_alreadyTransitive)
    if (NOT ${_fm_dep}_FOUND)
      list(APPEND _fm_missing_deps ${_fm_dep})
    endif()
  endforeach()
  if (NOT "${_fm_missing_deps}" STREQUAL "")
    set(_fm_fphsa_extra_args
      REASON_FAILURE_MESSAGE "missing dependencies: ${_fm_missing_deps}"
    )
    unset(_fm_missing_deps)
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MARLEY
  REQUIRED_VARS MARLEY_CONFIG_EXECUTABLE MARLEY_INCLUDE_DIR MARLEY_LIBRARIES
  ${_fm_fphsa_extra_required_vars}
  VERSION_VAR MARLEY_VERSION
  ${_fm_fphsa_extra_args}
)
unset(_fm_fphsa_extra_required_vars)
unset(_fm_fphsa_extra_args)

if (MARLEY_FOUND)
  foreach (_fm_lib IN LISTS _fm_libs)
    if (${_fm_lib}_LIBRARY AND NOT TARGET MARLEY::${_fm_lib})
      add_library(MARLEY::${_fm_lib} UNKNOWN IMPORTED)
      set_target_properties(MARLEY::${_fm_lib} PROPERTIES
        IMPORTED_LOCATION "${${_fm_lib}_LIBRARY}"
        INTERFACE_COMPILE_FEATURES "cxx_std_${MARLEY_CXX_STANDARD}"
        INTERFACE_INCLUDE_DIRECTORIES "${MARLEY_INCLUDE_DIR}")
      if (_fm_lib STREQUAL "MARLEY_ROOT")
        set_property(TARGET MARLEY::${_fm_lib}
          APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
          "${MARLEY_DEFINITIONS}")
      endif()
      if (_fm_transitive_deps_${fm_lib})
        set_property(TARGET MARLEY::${_fm_lib}
          APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${_fm_transitive_deps_${fm_lib}}")
      endif()
      unset(_fm_transitive_deps_${_fm_lib})
    endif()
  endforeach()
  foreach (_fm_prog IN LISTS _fm_progs)
    if (NOT TARGET MARLEY::${_fm_prog})
      string(MAKE_C_IDENTIFIER "${_fm_prog}" _fm_prog_var)
      string(TOUPPER "${_fm_prog_var}_EXECUTABLE" _fm_prog_var)
      add_executable(MARLEY::${_fm_prog} IMPORTED)
      set_target_properties(MARLEY::${_fm_prog} PROPERTIES
        IMPORTED_LOCATION "${${_fm_prog_var}}"
      )
    endif()
  endforeach()
endif()

unset(_fm_lib)
unset(_fm_libs)
unset(_fm_prog)
unset(_fm_progs)
