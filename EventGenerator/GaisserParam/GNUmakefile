#
# $Id: GNUmakefile,v 1.5 2009/09/29 21:10:31 brebel Exp $
#
# messier@indiana.edu
#
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_GaisserParam.mk

#VERBOSE     := true
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)GaisserParam
JOBFILES    := $(wildcard *.fcl)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(NUTOOLS_LIB) -lEventGeneratorBaseGaisserParam -lSimulationBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lSummaryData 
