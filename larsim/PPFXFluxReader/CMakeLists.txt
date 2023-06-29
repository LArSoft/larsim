


cet_build_plugin( PPFXFluxReader art::source
  LIBRARIES PRIVATE
        larsim::PPFXFluxReader
  	nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_IO_Sources
                ROOT::Tree
             )

cet_make_library(LIBRARY_NAME FluxInterface
  SOURCE
                 DK2NuInterface.cxx
                 GSimpleInterface.cxx
  LIBRARIES PUBLIC
	       nusimdata::SimulationBase
             	log4cpp::log4cpp
               	GENIE::GFwEG
               	GENIE::GFwAlg
               	GENIE::GFwMsg
               	GENIE::GFwGHEP
               	GENIE::GTlFlx
               	dk2nu::Tree
               	dk2nu::Genie
                ROOT::Tree
  PRIVATE
                      cetlib::cetlib
                      cetlib_except::cetlib_except
		      fhiclcpp::fhiclcpp
                      messagefacility::MF_MessageLogger
		      CLHEP::CLHEP
)


cet_make_library(
   SOURCE PPFXFluxReader.cxx
  LIBRARIES PUBLIC
  	nusimdata::SimulationBase
	  		larsim::FluxInterface
        lardata::Utilities
			larcoreobj::SummaryData
                        nugen::EventGeneratorBase_GENIE
                        art_root_io::TFileService_service
                        canvas::canvas
			messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                 	dk2nu::Tree
                 	dk2nu::Genie
                        ROOT::Hist
  PRIVATE
                        art::Framework_Core
                        art::Framework_IO_Sources
                   	art::Utilities
                   	art::Persistency_Common
                  	art::Persistency_Provenance
                        art::Framework_Principal
                        art::Framework_Services_Registry
        )                                                           

install_headers()
install_source()
install_fhicl()
