# Build script for the TriggerAlgo package.

# For a standard analysis task, you probably won't need to change these lines:
PACKAGE     := TriggerAlgo
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

# The SRT/LArSoft standard extensions. 
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

# OpticalDetector specific data structure

override CPPFLAGS += -I$(LARSOFT)/include

# This line you may have to change: If you require a library in addition
# to those listed, you'd add "-l<library-name>" to the end of this line.
# For example, if you had to read in raw::RawDigit objects from an event
# record, you'd add -lRawData to the end of the line, since raw::RawDigit
# is defined in package RawData.
override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(NUTOOLS_LIB) -lSimulationBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lSimulation -lGeometry -lRecoBase -lPhotonPropagation 
