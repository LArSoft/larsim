art_make(NO_PLUGINS
  EXCLUDE
    TFLoaderGAN_tool.cc
  LIB_LIBRARIES
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    cetlib
    cetlib_except
    fhiclcpp
    ${CLHEP}
    ${ART_UTILITIES}
    ${PROTOBUF}
    ${TENSORFLOW}
    ${FOR_TENSORFLOW}
  )

simple_plugin(TFLoaderGAN "tool"
              larsim_PhotonPropagation_TFLoaderTools
              )

install_headers()
install_fhicl()
install_source()

if( DEFINED ENV{TENSORFLOW_DIR} )
  if ( "${OSTYPE}" MATCHES "slf6" )
    set(FOR_TENSORFLOW "-lrt")
  endif ()
endif ()

include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories ( $ENV{PROTOBUF_INC} )
include_directories ( $ENV{TENSORFLOW_INC} )
include_directories ( $ENV{TENSORFLOW_INC}/eigen )
include_directories ( $ENV{TENSORFLOW_INC}/absl )
