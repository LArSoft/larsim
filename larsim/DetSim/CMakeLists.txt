simple_plugin(SimWireAna "module"
              lardataobj_RawData
              art::Framework_Services_Registry
              art_root_io::tfile_support
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Hist)

simple_plugin(SimWire "module"
              lardataalg_DetectorInfo
              lardataobj_RawData
              lardataobj_Simulation
              larcorealg_Geometry
              lardata_Utilities_LArFFT_service
              nurandom_RandomUtils_NuRandomService_service
              art::Framework_Services_Registry
              art_root_io::TFileService_service
              art_root_io::tfile_support
              ${CLHEP_LIBRARIES}
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Hist
              ROOT::MathCore)

simple_plugin(WienerFilterAna "module"
              larcorealg_Geometry
              art::Framework_Services_Registry
              art_root_io::tfile_support
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Hist)

install_headers()
install_fhicl()
install_source()
