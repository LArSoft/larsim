cet_make_library(LIBRARY_NAME TFLoaderTool PUBLIC
  SOURCE TFLoader.cc
  LIBRARIES PRIVATE
  TensorFlow::cc
  TensorFlow::framework
)

cet_write_plugin_builder(lar::TFLoaderTool art::tool Modules
  INSTALL_BUILDER
  larsim::TFLoaderTool)

include(lar::TFLoaderTool)
cet_build_plugin(TFLoaderMLP lar::TFLoaderTool
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
  TensorFlow::framework
)

install_headers()
install_fhicl()
install_source()
