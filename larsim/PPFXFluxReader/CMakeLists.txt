#cet_enable_asserts()

#include_directories ( $ENV{DK2NUGENIE_INC}/GENIE )
#include_directories ( $ENV{DK2NUDATA_INC} )
#include_directories ( $ENV{NUSIMDATA_INC} )
#include_directories ( $ENV{ROOT_INC} )
#include_directories ( $ENV{ART_ROOT_IO_INC} )
#link_directories ( $ENV{DK2NUDATA_LIB} )

file(GLOB INTERFACEFILES *Interface* )
simple_plugin( PPFXFluxReader "source"
	       PPFXFluxReader
               # ${ART_FRAMEWORK_PERSISTENCY_COMMON}
               BASENAME_ONLY
             )

art_make_library(LIBRARY_NAME FluxInterface 
                 SOURCE	      ${INTERFACEFILES}
                 LIBRARIES
		       nusimdata::SimulationBase
#GENIE STUFF########################################
			gsl
			xml2
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
                      #${ROOT_BASIC_LIB_LIST}
                      ROOT::Tree
                      ${Boost_SYSTEM_LIBRARY}
	        ) 

art_make( EXCLUDE PPFXFluxReader_source.cc ${INTERFACEFILES}
	  BASENAME_ONLY
          LIBRARY_NAME  PPFXFluxReader
          LIB_LIBRARIES	nusimdata::SimulationBase
	  		FluxInterface
			larcoreobj::SummaryData
                        art::Framework_Core
                        art_Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        #${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        #${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                   	art::Persistency_Common
                  	art::Persistency_Provenance
                   	art_Utilities
                        art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
                        #${MF_UTILITIES}
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        #${ROOT_BASIC_LIB_LIST}

        )                                                           

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
