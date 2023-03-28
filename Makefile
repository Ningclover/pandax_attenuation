
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lMinuit
GLIBS         = $(ROOTGLIBS)
LIBS         += -lstdc++ -lz

IDIR          = include
SRCDIR        = src
BUILDDIR      = build
TARGETDIR     = bin


_DEPS         = General_function.h
DEPS          = $(patsubst %,$(IDIR)/%,$(_DEPS))


TARGETS = $(TARGETDIR)/ADM_MC_back $(TARGETDIR)/ADM_MC_front

SRCEXT := C
SRCF :=$(SRCDIR)/ADM_MC_front.$(SRCEXT)
SRCB :=$(SRCDIR)/ADM_MC_back.$(SRCEXT)


OBJF := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SRCF:.$(SRCEXT)=.o))
OBJB := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SRCB:.$(SRCEXT)=.o))

.PHONY: all  clean

all: $(TARGETS)

$(TARGETDIR)/ADM_MC_back: $(OBJB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(GLIBS)

$(TARGETDIR)/ADM_MC_front: $(OBJF)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(GLIBS)

$(OBJB): $(SRCB) $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJF): $(SRCF) $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	@rm -f $(OBJF) $(OBJB) $(TARGETS)  *~  core

