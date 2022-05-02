cet_make_library(LIBRARY_NAME ScintTimeTool INTERFACE
  SOURCE ScintTime.h)

cet_write_plugin_builder(lar::ScintTimeTool art::tool Modules
  INSTALL_BUILDER)

include (lar::ScintTimeTool)

cet_build_plugin(ScintTimeLAr lar::ScintTimeTool
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
  CLHEP::Random
)

install_headers()
install_fhicl()
install_source()
