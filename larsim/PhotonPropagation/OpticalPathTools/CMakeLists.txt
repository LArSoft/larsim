cet_make_library(LIBRARY_NAME OpticalPath INTERFACE
  SOURCE OpticalPath.h)

cet_write_plugin_builder(phot::OpticalPath art::tool Modules
  INSTALL_BUILDER)

include (phot::OpticalPath)

cet_build_plugin(StandardOpticalPath phot::OpticalPath
  LIBRARIES PRIVATE
  larsim::OpticalPath
  larcoreobj::geo_vectors
)

install_headers()
install_fhicl()
install_source()
