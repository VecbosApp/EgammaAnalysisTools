ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./ 
INCLUDEDIRCOMMON = ../../../
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = ./lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

$(OUTLIB)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)EfficiencyEvaluator.o $<
$(OUTLIB)HZZEleIDSelector.o: $(INCLUDEDIRCOMMON)EgammaAnalysisTools/src/HZZEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)HZZEleIDSelector.o $<

fakerate: fakerate.cc estimateFakeRate.C \
	$(OUTLIB)EfficiencyEvaluator.o \
	$(OUTLIB)HZZEleIDSelector.o
	$(CXX) $(CXXFLAGS) -o fakerate $(OUTLIB)/*.o $(GLIBS) $ $<

fakeratehzz4l: fakerateHzz4l.cc estimateMuonFakeRateHzz4lTree.C \
	$(OUTLIB)EfficiencyEvaluator.o 
	$(CXX) $(CXXFLAGS) -o fakeratehzz4l $(OUTLIB)/*.o $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o
	rm -f fakerate
	rm -f fakeratehzz4l

all: fakerate

