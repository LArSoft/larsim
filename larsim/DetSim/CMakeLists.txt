simple_plugin(SimWireAna "module"
              lardataobj_RawData
              art::Framework_Services_Registry
              art_root_io::tfile_support
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Hist) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(SimWire "module"
              lardataalg_DetectorInfo
              lardataobj_RawData
              lardataobj_Simulation
              larcorealg_Geometry
              lardata_Utilities_LArFFT_service
              nurandom::RandomUtils_NuRandomService_service
              art::Framework_Services_Registry
              art_root_io::TFileService_service
              art_root_io::tfile_support
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Hist
              ROOT::MathCore) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(WienerFilterAna "module"
              larcorealg_Geometry
              art::Framework_Services_Registry
              art_root_io::tfile_support
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Hist) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
