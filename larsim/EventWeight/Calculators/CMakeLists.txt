# Headers needed for the GenieWeightCalc class
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

# Define a macro that indicates whether a GENIE version before 3.0.0 is in use
parse_ups_version( ${GENIE_VERSION} )
if( ${VMAJ} LESS 3 )
  add_definitions(-DGENIE_PRE_R3)
endif()

# Builds of MicroBooNE-specific larsim releases should uncomment the
# "add_definitions" line below. It enables extra reweighting functionality
# which is present only in "MicroBooNE patch" releases of GENIE v3.0.4. These
# require building against a special ups product.
#add_definitions(-DGENIE_UB_PATCH)

art_make(LIB_LIBRARIES
         larsim_EventWeight_Base
         nurandom_RandomUtils_NuRandomService_service
         nugen_EventGeneratorBase_GENIE
         art::Framework_Principal
         art::Persistency_Provenance
         art::Framework_Services_Registry
         fhiclcpp::fhiclcpp fhiclcpp::types
         ${GENIE_LIB_LIST}
         ${LOG4CPP}
         messagefacility::MF_MessageLogger
         cetlib_except)

install_headers()
install_fhicl()
install_source()
