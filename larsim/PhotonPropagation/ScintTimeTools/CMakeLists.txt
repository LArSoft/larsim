cet_build_plugin(ScintTimeLAr art::tool
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
  CLHEP::Random
)

install_headers()
install_fhicl()
install_source()
