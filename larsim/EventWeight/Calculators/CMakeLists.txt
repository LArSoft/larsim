include_directories ( $ENV{PPFX_DIR}/include )

cet_make_library(SOURCE GenieWeightCalc.cxx
                        PPFXCVWeightCalc.cxx
                        PPFXMIPPKaonWeightCalc.cxx
                        PPFXMIPPPionWeightCalc.cxx
                        PPFXOtherWeightCalc.cxx
                        PPFXTargAttenWeightCalc.cxx
                        PPFXThinKaonWeightCalc.cxx
                        PPFXThinMesonWeightCalc.cxx
                        PPFXThinNeutronPionWeightCalc.cxx
                        PPFXThinNucAWeightCalc.cxx
                        PPFXThinNucWeightCalc.cxx
                        PPFXThinPionWeightCalc.cxx
                        PPFXTotAbsorpWeightCalc.cxx
                        PPFXWeightCalc.cxx
  LIBRARIES
  PRIVATE
  ppfx::ppfx
  larsim::EventWeight_Base
  nugen::EventGeneratorBase_GENIE
  nusimdata::SimulationBase
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ${GENIE_LIB_LIST}
  dk2nu::Tree
  dk2nu::Genie
  CLHEP::Random
  log4cpp::log4cpp # FIXME Should be transitive from GENIE target(s).
)

find_file(GENIE_HAS_UB_PATCH
  NAMES GENIE/RwCalculators/GReWrightDeltaradAngle.h
  HINTS ${GENIE_INCLUDE_DIRS} NO_DEFAULT_PATH)
mark_as_advanced(GENIE_HAS_UB_PATCH)
if (GENIE_HAS_UB_PATCH)
  # Enable extra reweighting functionality which is present only in
  # "MicroBooNE patch" releases of GENIE v3.0.4. These require building
  # against a special ups product.
  target_compile_definitions(larsim_EventWeight_Calculators
    PRIVATE GENIE_UB_PATCH)
endif()

install_headers()
install_fhicl()
install_source()
