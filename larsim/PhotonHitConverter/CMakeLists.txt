cet_build_plugin(PhotonHitConverter art::module
      SOURCE PhotonHitConverter_module.cc
      LIBRARIES PRIVATE
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardataobj::Simulation
        artg4tk::pluginDetectors_gdml
)

install_source()
