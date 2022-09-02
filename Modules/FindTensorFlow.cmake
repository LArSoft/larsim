#[============================================================[.rst:
FindTensorFlow
============
#]============================================================]
if (UNIX)
  set(_ftf_m m)
  if (APPLE)
    set(_ftf_rt)
  else()
    set(_ftf_rt rt)
  endif()
else()
  set(_ftf_m)
  set(_ftf_rt)
endif()
set(_ftf_libs "_framework" "" "_cc")
set(_ftf_transitive_deps TensorFlow::framework)
set(_ftf_transitive_deps_cc TensorFlow::framework)
set(_ftf_transitive_deps_framework protobuf::libprotobuf Eigen3::Eigen Threads::Threads ZLIB::ZLIB ${_ftf_m} ${CMAKE_DL_LIBS})
set(_ftf_deps GIF gRPC JPEG PNG Protobuf SQLite3 Threads ZLIB)
unset(_ftf_m)
unset(_ftf_rt)

unset(_ftf_fphsa_extra_args)
if (NOT TensorFlow_FOUND)
  find_path(TensorFlow_INCLUDE_DIR tensorflow/core/public/version.h PATH_SUFFIXES include)
  mark_as_advanced(TensorFlow_INCLUDE_DIR)
  if (TensorFlow_INCLUDE_DIR)
    file(STRINGS "${TensorFlow_INCLUDE_DIR}/tensorflow/core/public/version.h"
      _ftf_version_bits REGEX "^[ \t]*#[ \t]*define[ \t]+TF_(MAJOR|MINOR|PATCH)_VERSION[ \t]+")
    list(TRANSFORM _ftf_version_bits REPLACE "^.*[ \t]+([^ \t]+)$" "\\1")
    list(JOIN _ftf_version_bits "." TensorFlow_VERSION)
    unset(_ftf_version_bits)
    foreach (_ftf_lib IN LISTS _ftf_libs)
      find_library(TensorFlow${_ftf_lib}_LIBRARY NAMES tensorflow${_ftf_lib})
      mark_as_advanced(TensorFlow${_ftf_lib}_LIBRARY)
      if (TensorFlow${_ftf_lib}_LIBRARY)
        if (_ftf_lib)
          string(REGEX REPLACE "^_" "TensorFlow::" _ftf_tgt "${_ftf_lib}")
        else()
          set(_ftf_tgt TensorFlow::TensorFlow)
        endif()
        list(APPEND TensorFlow_LIBRARIES ${_ftf_tgt})
      endif()
    endforeach()
  endif()
endif()

list(TRANSFORM _ftf_deps APPEND _FOUND
  OUTPUT_VARIABLE _ftf_fphsa_extra_required_vars)

if (TensorFlow_FOUND OR (TensorFlow_LIBRARIES AND TensorFlow_INCLUDE_DIR))
  unset(_ftf_missing_deps)
  foreach (_ftf_dep IN LISTS _ftf_deps)
    get_property(_ftf_${_ftf_dep}_alreadyTransitive GLOBAL PROPERTY
      _CMAKE_${_ftf_dep}_TRANSITIVE_DEPENDENCY)
    find_package(${_ftf_dep} ${_ftf_fp_${_ftf_dep}_args} QUIET)
    if (NOT DEFINED cet_${_ftf_dep}_alreadyTransitive OR cet_${_ftf_dep}_alreadyTransitive)
      set_property(GLOBAL PROPERTY _CMAKE_${_ftf_dep}_TRANSITIVE_DEPENDENCY TRUE)
    endif()
    unset(_ftf_${_ftf_dep}_alreadyTransitive)
    if (NOT ${_ftf_dep}_FOUND)
      list(APPEND _ftf_missing_deps ${_ftf_dep})
    endif()
  endforeach()
  if (NOT "${_ftf_missing_deps}" STREQUAL "")
    set(_ftf_fphsa_extra_args
      REASON_FAILURE_MESSAGE "missing dependencies: ${_ftf_missing_deps}"
    )
    unset(_ftf_missing_deps)
  endif()
endif()
unset(_ftf_dep)
unset(_ftf_deps)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(TensorFlow
  REQUIRED_VARS TensorFlow_INCLUDE_DIR TensorFlow_LIBRARIES
  ${_ftf_fphsa_extra_required_vars}
  VERSION_VAR TensorFlow_VERSION
  ${_ftf_fphsa_extra_args}
)
unset(_ftf_fphsa_extra_required_vars)
unset(_ftf_fphsa_extra_args)

if (TensorFlow_FOUND)
  foreach (_ftf_lib IN LISTS _ftf_libs)
    if (_ftf_lib)
      string(REGEX REPLACE "^_" "TensorFlow::" _ftf_tgt "${_ftf_lib}")
    else()
      set(_ftf_tgt TensorFlow::TensorFlow)
    endif()
    if (TensorFlow${_ftf_lib}_LIBRARY AND NOT TARGET ${_ftf_tgt})
      add_library(${_ftf_tgt} SHARED IMPORTED)
      set_target_properties(${_ftf_tgt} PROPERTIES
        IMPORTED_LOCATION "${TensorFlow${_ftf_lib}_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${TensorFlow_INCLUDE_DIR}")
      foreach (_ftf_tdep IN LISTS _ftf_transitive_deps${_ftf_lib})
        set_property(TARGET ${_ftf_tgt}
          APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${_ftf_tdep}")
      endforeach()
      unset(_ftf_tdep)
      unset(_ftf_transitive_deps_${_ftf_lib})
    endif()
  endforeach()
endif()

unset(_ftf_lib)
unset(_ftf_libs)
