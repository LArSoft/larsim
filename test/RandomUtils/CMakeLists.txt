# using art_make here would interferes with cet_test non-HANDBUILT targets
simple_plugin( SeedTest01 "module"
    larsim_RandomUtils_LArSeedService_service
    art_Framework_Core
    art_Framework_Principal
    art_Framework_Services_Registry
    art_Framework_Services_Optional_RandomNumberGenerator_service
    MF_MessageLogger
    MF_Utilities
    fhiclcpp
    cetlib
    ${CLHEP}
  NO_INSTALL
)

simple_plugin( RandomManagerTest "module"
    larsim_RandomUtils_LArSeedService_service
    art_Framework_Core
    art_Framework_Principal
    art_Framework_Services_Registry
    art_Framework_Services_Optional_RandomNumberGenerator_service
    MF_MessageLogger
    MF_Utilities
    fhiclcpp
    cetlib
    ${CLHEP}
  NO_INSTALL
)

simple_plugin( GeneratedEventTimestampTestPlugin "plugin"
    art_Framework_Core
    art_Persistency_Provenance
    MF_MessageLogger
    MF_Utilities
    cetlib
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
  BASENAME_ONLY
)


#
# Some tests are going to fail at configuration phase. Those are "failing".
# Some other tests will have failures when getting a seed:
#   the errors are expected and those tests are considered successful.
# Some more are art-specific and don't test anything in SeedMaster.
#

# This batch of tests is for LArSeedService ("integration" tests)
set( SuccessfulServiceTestIndices 01 02 03 04 05 10 11 12 13 14 15 16 17 18 19 21 22 23 24 )
foreach( ServiceTestIndex ${SuccessfulServiceTestIndices} )
  cet_test( test${ServiceTestIndex} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestIndex}.fcl
            DATAFILES Providers/test${ServiceTestIndex}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestIndex )

# This batch of tests is for LArSeedService ("integration" tests) that do not have a corresponding unit test
set( SuccessfulServiceTestIndices 06 07 08 25 )
foreach( ServiceTestIndex ${SuccessfulServiceTestIndices} )
  cet_test( test${ServiceTestIndex} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestIndex}.fcl
            DATAFILES test${ServiceTestIndex}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestIndex )


set( FailingServiceTestIndices 09 20 )
foreach( ServiceTestIndex ${FailingServiceTestIndices} )
  cet_test( test${ServiceTestIndex} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config test${ServiceTestIndex}.fcl
            TEST_PROPERTIES WILL_FAIL true
            DATAFILES Providers/test${ServiceTestIndex}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceTestIndex )


# This batch of tests is for LArSeedService engine manager features
set( ServiceManagingTestIndices 01 03 )
foreach( ServiceManagingTestIndex ${ServiceManagingTestIndices} )
  cet_test( RandomManagerTest${ServiceManagingTestIndex} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config RandomManagerTest${ServiceManagingTestIndex}.fcl
            DATAFILES RandomManagerTest${ServiceManagingTestIndex}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceManagingTestIndex )

set( FailingServiceManagingTestIndices 02 )
foreach( ServiceManagingTestIndex ${FailingServiceManagingTestIndices} )
  cet_test( RandomManagerTest${ServiceManagingTestIndex} HANDBUILT
            TEST_EXEC art
            TEST_ARGS --rethrow-all --config RandomManagerTest${ServiceManagingTestIndex}.fcl
            TEST_PROPERTIES WILL_FAIL true
            DATAFILES RandomManagerTest${ServiceManagingTestIndex}.fcl Providers/messageService.fcl Providers/standardMessageDestinations.fcl
          )
endforeach( ServiceManagingTestIndex )


add_subdirectory(Providers)
