cet_build_plugin(CosmicsGen art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  larcoreobj::SummaryData
  nutools::EventGeneratorBase_CRY
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Hist
)
         
# For evgenbase.h:
if (TARGET nugen::EventGeneratorBase)
  target_link_libraries(larsim_EventGenerator_CRY_CosmicsGen_module
    PRIVATE
    nugen::EventGeneratorBase
  )
else()
  target_link_libraries(larsim_EventGenerator_CRY_CosmicsGen_module
    PRIVATE
    ROOT::MathCore
  )
  target_include_directories(larsim_EventGenerator_CRY_CosmicsGen_module
    BEFORE PRIVATE ${nugen_INCLUDE_DIR}
  )
endif()

install_headers()
install_fhicl()
install_source()
