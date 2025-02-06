cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH ${TRANSITIVE_PATHS_WITH_FHICL_DIR})

cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})

cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})

#
# Test GenerateEventTimestamp plug-in with SeedService's per-event policy,
# or test SeedService's per-event policy with GenerateEventTimestamp plug-in.
#
cet_build_plugin(TestGeneratedEventTimestamp art::EDAnalyzer NO_INSTALL
  LIBRARIES PRIVATE
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Framework_Services_Optional_RandomNumberGenerator_service
  art::Persistency_Provenance
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)

#
# GeneratedEventTimestamp tests
# ------------------------------
#
# Tests 1_1 and 1_2 are the same.
# They are expected to yield different results since the time stamps are
# different.
#
cet_test(GeneratedEventTimestamp_test1_1 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config GeneratedEventTimestamp_test1.fcl -n 10
  DATAFILES
    GeneratedEventTimestamp_test1.fcl
)

cet_test(GeneratedEventTimestamp_test1_2 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config GeneratedEventTimestamp_test1.fcl -n 10
  DATAFILES
    GeneratedEventTimestamp_test1.fcl
)

#
# Test 2 should fail since no time stamp is set.
#
cet_test(GeneratedEventTimestamp_test2 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config GeneratedEventTimestamp_test2.fcl -n 10
  TEST_PROPERTIES WILL_FAIL true
  DATAFILES
    GeneratedEventTimestamp_test2.fcl
)

#
# single particle production test
#
# kludge for now to only run with mrb
set( mrb_build_dir $ENV{MRB_BUILDDIR} )
if( mrb_build_dir )
cet_test(gensingle HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config gensingle_test.fcl -n 2
  DATAFILES gensingle_test.fcl
)
endif( mrb_build_dir )
