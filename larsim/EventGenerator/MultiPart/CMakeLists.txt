cet_build_plugin( MultiPartRain art::module
        LIBRARIES PRIVATE
            larcorealg::Geometry
            lardata::ArtDataHelper
            larsim::Simulation nug4::ParticleNavigation
            lardataobj::Simulation
            lardata::Utilities
            lardataobj::MCBase
            larcoreobj::SummaryData
            nusimdata::SimulationBase
            nurandom::RandomUtils_NuRandomService_service
            CLHEP::CLHEP
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support ROOT::Core
            art_root_io::TFileService_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            ROOT::EG
        )

cet_build_plugin( MultiPartVertex art::module
        LIBRARIES PRIVATE
            larcorealg::Geometry
            lardata::ArtDataHelper
            larsim::Simulation nug4::ParticleNavigation
            lardataobj::Simulation
            lardata::Utilities
            lardataobj::MCBase
            larcoreobj::SummaryData
            nusimdata::SimulationBase
            nurandom::RandomUtils_NuRandomService_service
            CLHEP::CLHEP
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support ROOT::Core
            art_root_io::TFileService_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            ROOT::EG
        )

install_source()

