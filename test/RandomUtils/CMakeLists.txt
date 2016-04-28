add_subdirectory(Providers)

art_make(
  EXCLUDE GeneratedEventTimestampTestPlugin_plugin.cc
  LIB_LIBRARIES 
    larsim_RandomUtils_LArSeedService_service
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_UTILITIES}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${FHICLCPP}
    ${CETLIB}
    ${CLHEP}
  MODULE_LIBRARIES
    larsim_test_RandomUtils
    larsim_RandomUtils_LArSeedService_service
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${FHICLCPP}
    ${CETLIB}
    ${CLHEP}
  SERVICE_LIBRARIES
    larsim_RandomUtils_LArSeedService_service
    ${ART_UTILITIES}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${FHICLCPP}
#     ${CETLIB}
    ${ROOT_CORE}
    ${ROOT_MATHCORE}
  NO_INSTALL
  USE_PRODUCT_NAME
)


#
# The following tests verify that engine declaration may happen only on
# construction, and that seed query can happen at any time.
#
cet_test( SeedTestQuerySeeds_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c seedtest_query_seeds.fcl
  DATAFILES seedtest_query_seeds.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
)

cet_test( SeedTestRegisterOnProcessEvent_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c seedtest_register_on_process_event.fcl
  DATAFILES seedtest_register_on_process_event.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
  TEST_PROPERTIES WILL_FAIL true
)

cet_test( SeedTestRegisterOnSubRun_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c seedtest_register_on_subrun.fcl
  DATAFILES seedtest_register_on_subrun.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
  TEST_PROPERTIES WILL_FAIL true
)

cet_test( SeedTestRegisterOnRunRegistered_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c seedtest_register_on_run_registered.fcl
  DATAFILES seedtest_register_on_run_registered.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
  TEST_PROPERTIES WILL_FAIL false
)

cet_test( SeedTestRegisterOnRunUnregistered_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c seedtest_register_on_run_unregistered.fcl
  DATAFILES seedtest_register_on_run_unregistered.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
  TEST_PROPERTIES WILL_FAIL true
)

cet_test( GlobalSeedTestLinear_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c globalseedtest_linear.fcl
  DATAFILES globalseedtest_linear.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
)

cet_test( GlobalSeedTestPredefinedSeed_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c globalseedtest_predefined.fcl
  DATAFILES globalseedtest_predefined.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
)

cet_test( GlobalSeedTestPerEvent_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all -c globalseedtest_perevent.fcl
  DATAFILES globalseedtest_perevent.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
)

#
# Some tests are going to fail at configuration phase. Those are "failing".
# Some other tests will have failures when getting a seed:
#   the errors are expected and those tests are considered successful.
# Some more are art-specific and don't test anything in SeedMaster.
#


# The first batch of tests is for LArSeedService ("integration" tests shared with SeedMaster)
set( SuccessfulServiceSharedTests
  LinearMap01
  LinearMapDepr01
  PredefinedOfs01
  PredefinedOfs02
  PredefinedSeed01
  PredefinedSeed02
  PredefinedSeed03
  PredefinedSeed04
  Random01
  Random02
  )
foreach( ServiceTestName ${SuccessfulServiceSharedTests} )
  cet_test( test${ServiceTestName} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestName}.fcl
            DATAFILES Providers/test${ServiceTestName}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestName )

set( FailingServiceSharedTests
  InvalidPolicy
  LinearMapErr01
  PredefinedOfsErr01
  PredefinedOfsErr02
  PredefinedOfsErr04
  PredefinedOfsErr05
  PredefinedSeedErr01
  PredefinedSeedErr02
  PredefinedSeedErr03
  PredefinedOfsErr03
  )
foreach( ServiceTestName ${FailingServiceSharedTests} )
  cet_test( test${ServiceTestName} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestName}.fcl
            TEST_PROPERTIES WILL_FAIL true
            DATAFILES Providers/test${ServiceTestName}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestName )


# The second batch of tests is for LArSeedService ("integration" tests for art service only)
set( SuccessfulServiceOnlyTests
  PerEvent01
  )
foreach( ServiceTestName ${SuccessfulServiceOnlyTests} )
  cet_test( test${ServiceTestName} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestName}.fcl
            DATAFILES test${ServiceTestName}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestName )

set( FailingServiceOnlyTests
  PerEventErr01
  )
foreach( ServiceTestName ${FailingServiceOnlyTests} )
  cet_test( test${ServiceTestName} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestName}.fcl
            TEST_PROPERTIES WILL_FAIL true
            DATAFILES test${ServiceTestName}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestName )


# The third batch of tests is for LArSeedService engine manager features
set( ServiceManagingTests
  RandomManagerTest01
  RandomManagerTest03
  )
foreach( ServiceManagingTestName ${ServiceManagingTests} )
  cet_test( ${ServiceManagingTestName} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config ${ServiceManagingTestName}.fcl
            DATAFILES ${ServiceManagingTestName}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceManagingTestName )

set( FailingServiceManagingTests
  RandomManagerTest02
  )
foreach( ServiceManagingTestName ${FailingServiceManagingTests} )
  cet_test( ${ServiceManagingTestName} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config ${ServiceManagingTestName}.fcl
            TEST_PROPERTIES WILL_FAIL true
            DATAFILES ${ServiceManagingTestName}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceManagingTestName )

