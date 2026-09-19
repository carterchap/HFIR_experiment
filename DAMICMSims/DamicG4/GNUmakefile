# $Id: GNUmakefile 68058 2013-03-13 14:47:43Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := DAMICG4
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

############
# ROOT
############
G4NOHIST := true CPPFLAGS += $(shell root-config --cflags -ggdb) LDFLAGS += $(shell root-config --libs -ggdb)

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

ifdef ROOTSYS
  CPPFLAGS += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags -g)
  LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

   ifeq ($(G4SYSTEM),Linux-g++)
     LDLIBS += -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir)
   endif
endif


