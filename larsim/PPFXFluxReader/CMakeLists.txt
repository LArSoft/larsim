#cet_enable_asserts()

file(GLOB INTERFACEFILES *Interface*.cxx )
simple_plugin( PPFXFluxReader "source"
	       PPFXFluxReader
               BASENAME_ONLY
             )

art_make_library(LIBRARY_NAME FluxInterface 
                 SOURCE	      ${INTERFACEFILES}
                 LIBRARIES
		       nusimdata::SimulationBase
#GENIE STUFF########################################
                 	log4cpp::log4cpp
                 	GENIE::GFwEG
                 	GENIE::GFwAlg
                 	GENIE::GFwMsg
                 	GENIE::GFwGHEP
                 	GENIE::GTlFlx
                 	dk2nu::Tree
                 	dk2nu::Genie
####################################################
                      cetlib::cetlib
                      cetlib_except::cetlib_except
		      fhiclcpp::fhiclcpp
                      messagefacility::MF_MessageLogger
		      CLHEP::CLHEP
                      ROOT::Tree
                      ${Boost_SYSTEM_LIBRARY}
	        ) 

art_make( EXCLUDE PPFXFluxReader_source.cc ${INTERFACEFILES}
	  BASENAME_ONLY
          LIBRARY_NAME  PPFXFluxReader
          LIB_LIBRARIES	nusimdata::SimulationBase
	  		FluxInterface
			larcoreobj::SummaryData
                        nugen::EventGeneratorBase_GENIE
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                   	art::Persistency_Common
                  	art::Persistency_Provenance
                   	art::Utilities
                        art::Utilities
                        art_root_io::TFileService_service
                        canvas::canvas
			messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                 	dk2nu::Tree
                 	dk2nu::Genie
                        ROOT::Hist
        )                                                           

install_headers()
install_source()
install_fhicl()
