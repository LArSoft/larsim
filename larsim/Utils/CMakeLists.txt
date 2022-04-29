cet_make_library(SOURCE
  SCEOffsetBounds.cxx
  TruthMatchUtils.cc
  LIBRARIES
  PUBLIC
  larcoreobj::geo_vectors
  canvas::canvas
  PRIVATE
  lardataobj::RecoBase
  larsim::MCCheater_BackTrackerService_service
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  canvas::canvas
)

install_headers()
install_source()
