cet_make_library(SOURCE GenericCRT.cxx
  LIBRARIES PUBLIC
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  art::Framework_Services_Registry
  PRIVATE
  messagefacility::MF_MessageLogger
  CLHEP::Evaluator
)

cet_build_plugin(GenericCRT art::EDProducer
  LIBRARIES PRIVATE
  larsim::GenericCRT
  lardataobj::Simulation
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()
