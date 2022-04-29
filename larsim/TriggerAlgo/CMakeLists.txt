cet_make_library(SOURCE TriggerAlgoBase.cc
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
