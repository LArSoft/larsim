#
# $Id: GNUmakefile,v 1.5 2009/09/29 21:10:31 echurch Exp $
#
# messier@indiana.edu
#
include SoftRelTools/arch_spec_root.mk
#include SoftRelTools/arch_spec_nuwro.mk

#VERBOSE     := true
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)NuWro
LIBCXXFILES := $(wildcard *.cxx) 
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE) 

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override CXXFLAGS += -Wno-sign-compare -Wno-reorder -Wno-unused-variable
override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(NUTOOLS_LIB) -lEventGeneratorNuWro -lSimulationBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lSummaryData -L$(NUWROHOME)/bin -levent1
