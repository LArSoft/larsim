#
# Test GenerateEventTimestamp plug-in with SeedService's per-event policy,
# or test SeedService's per-event policy with GenerateEventTimestamp plug-in.
#
simple_plugin( TestGeneratedEventTimestamp "module"
  ${SEEDSERVICE_SERVICE} # artextensions
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
# Tests 2_1 and 2_2 are the same.
# They are expected to yield to the same results, since the event timestamps
# will be the same on every run.
#
cet_test(GeneratedEventTimestamp_test2_1 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config GeneratedEventTimestamp_test2.fcl -n 10
  DATAFILES
    GeneratedEventTimestamp_test2.fcl
)

cet_test(GeneratedEventTimestamp_test2_2 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config GeneratedEventTimestamp_test2.fcl -n 10
  DATAFILES
    GeneratedEventTimestamp_test2.fcl
)

