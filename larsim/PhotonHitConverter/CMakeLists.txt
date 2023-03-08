set(
    MODULE_LIBRARIES
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardata::Utilities
        lardataobj::RecoBase
        larsim::Simulation
        lardataobj::Simulation
        artg4tk::pluginDetectors_gdml
        nusimdata::SimulationBase
        canvas::canvas
        cetlib::cetlib
        cetlib_except::cetlib_except
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        ROOT::Core
        ROOT::Tree
)

cet_build_plugin(PhotonHitConverter art::module SOURCE PhotonHitConverter_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
