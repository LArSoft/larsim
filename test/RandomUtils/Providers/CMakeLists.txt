# unit tests
cet_test( SeedMaster_test NO_AUTO
          LIBRARIES art_Utilities
                    MF_MessageLogger
                    MF_Utilities
                    fhiclcpp
                    cetlib
)

#
# Some tests are going to fail at configuration phase. Those are "failing".
# Some other tests will have failures when getting a seed:
#   the errors are expected and those tests are considered successful.
# Some more are art-specific and don't test anything in SeedMaster.
#

# This first batch of tests is for SeedMaster (no art involved, full unit test).
# skipped tests: 06 07 08 (art-specific) 25 (requires an event loop)
set( SuccessfulSeedMasterTestIndices 01 02 03 04 05 10 11 12 13 14 15 16 17 18 19 21 22 23 24 )
foreach( SeedMasterTestIndex ${SuccessfulSeedMasterTestIndices} )
  cet_test( SeedMaster_test${SeedMasterTestIndex} HANDBUILT
            TEST_EXEC SeedMaster_test
            TEST_ARGS test${SeedMasterTestIndex}.fcl "services.LArSeedService"
            DATAFILES test${SeedMasterTestIndex}.fcl messageService.fcl standardMessageDestinations.fcl
          )
endforeach( SeedMasterTestIndex )

set( FailingSeedMasterTestIndices 09 20 )
foreach( SeedMasterTestIndex ${FailingSeedMasterTestIndices} )
  cet_test( SeedMaster_test${SeedMasterTestIndex} HANDBUILT
            TEST_EXEC SeedMaster_test
            TEST_ARGS test${SeedMasterTestIndex}.fcl "services.LArSeedService"
            TEST_PROPERTIES WILL_FAIL true
            DATAFILES test${SeedMasterTestIndex}.fcl messageService.fcl standardMessageDestinations.fcl
          )
endforeach( SeedMasterTestIndex )

