cet_build_plugin(CosmicsGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  larcoreobj::SummaryData
  nugen::EventGeneratorBase
  nutools::EventGeneratorBase_CRY
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Hist
)

install_headers()
install_fhicl()
install_source()
