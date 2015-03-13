#
# Test GenerateEventTimestamp plug-in with SeedService's per-event policy,
# or test SeedService's per-event policy with GenerateEventTimestamp plug-in.
#
simple_plugin( TestGeneratedEventTimestamp "module"
  SeedService_service # artextensions
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
  ${ART_PERSISTENCY_PROVENANCE}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  ${CETLIB}
  ${CLHEP}
  NO_INSTALL
  BASENAME_ONLY
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

add_subdirectory(CRY)
# add_subdirectory(GENIE)
