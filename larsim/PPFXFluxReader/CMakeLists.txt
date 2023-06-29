cet_build_plugin(PPFXFluxReader art::source
  LIBRARIES PRIVATE
    larsim::PPFXFluxReader
)

cet_make_library(LIBRARY_NAME FluxInterface USE_PROJECT_NAME
  SOURCE
    DK2NuInterface.cxx
    GSimpleInterface.cxx
  LIBRARIES PUBLIC
    GENIE::GTlFlx
    ROOT::Physics
  PRIVATE
    CLHEP::CLHEP
    GENIE::GFwAlg
    GENIE::GFwEG
    GENIE::GFwGHEP
    GENIE::GFwMsg
    ROOT::Tree
    dk2nu::Genie
    dk2nu::Tree
    fhiclcpp::fhiclcpp
    log4cpp::log4cpp
    messagefacility::MF_MessageLogger
    nusimdata::SimulationBase
)


cet_make_library(
  SOURCE PPFXFluxReader.cxx
  LIBRARIES PUBLIC
    larsim::FluxInterface
    nusimdata::SimulationBase
  PRIVATE
    ROOT::Hist
    art::Framework_Core
    art::Framework_IO_Sources
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    canvas::canvas
    dk2nu::Genie
    dk2nu::Tree
    fhiclcpp::fhiclcpp
    larcoreobj::SummaryData
    lardata::Utilities
    messagefacility::MF_MessageLogger
    nugen::EventGeneratorBase_GENIE
)

install_headers()
install_source()
install_fhicl()
