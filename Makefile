#ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
#ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

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

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
CXX           = c++ -lm
#CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual $(ROOTCFLAGS) $(CLHEPCFLAGS)
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  
LD           = c++
LDFLAGS       = -O2 -g
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
ifeq ($(MACOSX_MINOR),4)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.4 c++
else
ifeq ($(MACOSX_MINOR),3)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.3 c++
else
UNDEFOPT      = suppress
LD            = c++
endif
endif

endif



NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
INCLUDEDIRCOMMON = ../
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)EgammaBase.o: $(INCLUDEDIR)/src/EgammaBase.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)EgammaBase.o $<
$(OUTLIB)Egamma.o: $(INCLUDEDIR)/src/Egamma.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Egamma.o $<
$(OUTLIBCOMMON)Conditions.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Conditions.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Conditions.o $<
$(OUTLIBCOMMON)Utils.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Utils.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIBCOMMON)Utils.o $<
$(OUTLIBCOMMON)Counters.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Counters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Counters.o $<
$(OUTLIBCOMMON)Selection.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Selection.o $<
$(OUTLIB)CutBasedEleIDSelector.o: $(INCLUDEDIR)/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutBasedEleIDSelector.o $<
$(OUTLIB)CiCBasedEleSelector.o: $(INCLUDEDIR)/src/CiCBasedEleSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CiCBasedEleSelector.o $<
$(OUTLIB)eIDCiChzzSelector.o: $(INCLUDEDIR)/src/eIDCiChzzSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)eIDCiChzzSelector.o $<
$(OUTLIB)EcalCleaner.o: $(INCLUDEDIR)/src/EcalCleaner.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)EcalCleaner.o $<
$(OUTLIBCOMMON)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EfficiencyEvaluator.o $<
$(OUTLIBCOMMON)FiguresOfMeritEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/FiguresOfMeritEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)FiguresOfMeritEvaluator.o $<
$(OUTLIBCOMMON)Monitor.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Monitor.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Monitor.o $<
$(OUTLIBCOMMON)SprDataFiller.o: $(INCLUDEDIRCOMMON)/CommonTools/src/SprDataFiller.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)SprDataFiller.o $<
$(OUTLIBCOMMON)TriggerMask.o: $(INCLUDEDIRCOMMON)/CommonTools/src/TriggerMask.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)TriggerMask.o $<
$(OUTLIB)ElectronIDMVA.o: $(INCLUDEDIR)/src/ElectronIDMVA.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronIDMVA.o $<
$(OUTLIB)ElectronIDMVAHZZ.o: $(INCLUDEDIR)/src/ElectronIDMVAHZZ.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronIDMVAHZZ.o $<
$(OUTLIB)HZZEleIDSelector.o: $(INCLUDEDIR)/src/HZZEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HZZEleIDSelector.o $<
$(OUTLIB)eIDSimpleCutsSelector.o: $(INCLUDEDIR)/src/eIDSimpleCutsSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)eIDSimpleCutsSelector.o $<
$(OUTLIB)EGammaMvaEleEstimator.o: $(INCLUDEDIR)/src/EGammaMvaEleEstimator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)EGammaMvaEleEstimator.o $<
$(OUTLIB)SimpleCutsIDSelector.o: $(INCLUDEDIR)/src/SimpleCutsIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleCutsIDSelector.o $<
$(OUTLIB)RedEleIDTree.o: $(INCLUDEDIR)/src/RedEleIDTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDTree.o $<
$(OUTLIB)FakeTree.o: $(INCLUDEDIR)/src/FakeTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FakeTree.o $<
$(OUTLIB)LikelihoodAnalysis.o: $(INCLUDEDIR)/src/LikelihoodAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LikelihoodAnalysis.o $<
$(OUTLIB)FakeElectronSelector.o: $(INCLUDEDIR)/src/FakeElectronSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FakeElectronSelector.o $<
$(OUTLIB)FakeElectronSelectorWenuPlusOneJet.o: $(INCLUDEDIR)/src/FakeElectronSelectorWenuPlusOneJet.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FakeElectronSelectorWenuPlusOneJet.o $<
$(OUTLIB)FakeElectronSelectorWmunuPlusOneJet.o: $(INCLUDEDIR)/src/FakeElectronSelectorWmunuPlusOneJet.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FakeElectronSelectorWmunuPlusOneJet.o $<
$(OUTLIB)FakeElectronSelectorZllPlusOneFake.o: $(INCLUDEDIR)/src/FakeElectronSelectorZllPlusOneFake.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FakeElectronSelectorZllPlusOneFake.o $<
$(OUTLIB)ZeeTagAndProbe.o: $(INCLUDEDIR)/src/ZeeTagAndProbe.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ZeeTagAndProbe.o $<
$(OUTLIB)HZZ4LElectronSelector.o: $(INCLUDEDIR)/src/HZZ4LElectronSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HZZ4LElectronSelector.o $<
$(OUTLIB)LHPdfsProducer.o: $(INCLUDEDIR)/src/LHPdfsProducer.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LHPdfsProducer.o $<
$(OUTLIB)sPlotsPdfsComparison.o: $(INCLUDEDIR)/src/sPlotsPdfsComparison.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)sPlotsPdfsComparison.o $<
$(OUTLIB)IsolationPdfsProducer.o: $(INCLUDEDIR)/src/IsolationPdfsProducer.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)IsolationPdfsProducer.o $<
#$(OUTLIB)SuperClusterWSelection.o: $(INCLUDEDIR)/src/SuperClusterWSelection.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SuperClusterWSelection.o $<
$(OUTLIB)McTruthEvent.o: ../src/McTruthEvent.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)McTruthEvent.o $<
$(OUTLIB)PFElectronSeedingEfficiency.o: $(INCLUDEDIR)/src/PFElectronSeedingEfficiency.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)PFElectronSeedingEfficiency.o $<
$(OUTLIB)PFElectronSeedingDistributions.o: $(INCLUDEDIR)/src/PFElectronSeedingDistributions.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)PFElectronSeedingDistributions.o $<
$(OUTLIB)LikelihoodPdf.o: $(INCLUDEDIR)/src/LikelihoodPdf.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LikelihoodPdf.o $<
$(OUTLIB)LikelihoodSpecies.o: $(INCLUDEDIR)/src/LikelihoodSpecies.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LikelihoodSpecies.o $<
$(OUTLIB)LikelihoodPdfProduct.o: $(INCLUDEDIR)/src/LikelihoodPdfProduct.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LikelihoodPdfProduct.o $<
$(OUTLIB)ElectronLikelihood.o: $(INCLUDEDIR)/src/ElectronLikelihood.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronLikelihood.o $<
$(OUTLIB)TestTurnOnCurveGamma.o: $(INCLUDEDIR)/src/TestTurnOnCurveGamma.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TestTurnOnCurveGamma.o $<
$(OUTLIB)TurnOnTreeGamma.o: $(INCLUDEDIR)/src/TurnOnTreeGamma.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TurnOnTreeGamma.o $<

#----------------------------------------------------#

# ==================== EgammaAnalysis =============================================
EgammaAnalysis:  $(INCLUDEDIR)/src/EgammaAnalysis.C \
	$(OUTLIB)EgammaBase.o \
	$(OUTLIB)Egamma.o \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIB)CiCBasedEleSelector.o \
	$(OUTLIB)eIDCiChzzSelector.o \
	$(OUTLIB)CutBasedEleIDSelector.o \
	$(OUTLIB)EcalCleaner.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIBCOMMON)FiguresOfMeritEvaluator.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)Monitor.o \
	$(OUTLIBCOMMON)SprDataFiller.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIB)ElectronIDMVA.o \
	$(OUTLIB)ElectronIDMVAHZZ.o \
	$(OUTLIB)HZZEleIDSelector.o \
	$(OUTLIB)eIDSimpleCutsSelector.o \
	$(OUTLIB)EGammaMvaEleEstimator.o \
	$(OUTLIB)SimpleCutsIDSelector.o \
	$(OUTLIB)RedEleIDTree.o \
	$(OUTLIB)FakeTree.o \
	$(OUTLIB)sPlotsPdfsComparison.o \
	$(OUTLIB)LHPdfsProducer.o \
	$(OUTLIB)IsolationPdfsProducer.o \
	$(OUTLIB)LikelihoodPdf.o \
	$(OUTLIB)LikelihoodSpecies.o \
	$(OUTLIB)LikelihoodPdfProduct.o \
	$(OUTLIB)ElectronLikelihood.o \
	$(OUTLIB)LikelihoodAnalysis.o \
	$(OUTLIB)FakeElectronSelector.o \
	$(OUTLIB)FakeElectronSelectorWenuPlusOneJet.o \
	$(OUTLIB)FakeElectronSelectorWmunuPlusOneJet.o \
	$(OUTLIB)FakeElectronSelectorZllPlusOneFake.o \
	$(OUTLIB)ZeeTagAndProbe.o \
	$(OUTLIB)TestTurnOnCurveGamma.o \
	$(OUTLIB)TurnOnTreeGamma.o \
	$(OUTLIB)HZZ4LElectronSelector.o
	$(CXX) $(CXXFLAGS) -o EgammaAnalysis $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

# ======= likelihood libs
Likelihood: $(OUTLIB)LikelihoodPdf.o \
	$(OUTLIB)LikelihoodSpecies.o \
	$(OUTLIB)LikelihoodPdfProduct.o \
	$(OUTLIB)ElectronLikelihood.o

CompareEff: $(INCLUDEDIR)/src/CompareEff.C \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o
	$(CXX) $(CXXFLAGS) -o CompareEff $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

CompareMisId: $(INCLUDEDIR)/src/CompareMisId.C \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o
	$(CXX) $(CXXFLAGS) -o CompareMisId $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

CompareMisId_analysisPlots: $(INCLUDEDIR)/src/CompareMisId_analysisPlots.C \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o
	$(CXX) $(CXXFLAGS) -o CompareMisId_analysisPlots $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

produceMCSubtractedFake: $(INCLUDEDIR)/src/produceMCSubtractedFake.cxx \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o
	$(CXX) $(CXXFLAGS) -o produceMCSubtractedFake $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

makeFriendTrees: $(INCLUDEDIR)/src/MakeFriendTrees.cc \
	$(OUTLIB)HZZEleIDSelector.o \
	$(OUTLIB)SimpleCutsIDSelector.o \
	$(OUTLIB)eIDSimpleCutsSelector.o \
	$(OUTLIB)EGammaMvaEleEstimator.o \
	$(OUTLIB)eIDCiChzzSelector.o
	$(CXX) $(CXXFLAGS) -o makeFriendTrees $(OUTLIB)/*o $(GLIBS) $ $<

CompareClasses: $(INCLUDEDIR)/src/CompareClasses.C \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o
	$(CXX) $(CXXFLAGS) -o CompareClasses $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

MakeNotePdfPlots: $(INCLUDEDIR)/src/MakeNotePdfPlots.C
	$(CXX) $(CXXFLAGS) -o MakeNotePdfPlots $(GLIBS) $ $<

MakeNoteBkgPdfPlots: $(INCLUDEDIR)/src/MakeNoteBkgPdfPlots.C
	$(CXX) $(CXXFLAGS) -o MakeNoteBkgPdfPlots $(GLIBS) $ $<

MakeNoteIsolationPlots: $(INCLUDEDIR)/src/MakeNoteIsolationPlots.C
	$(CXX) $(CXXFLAGS) -o MakeNoteIsolationPlots $(GLIBS) $ $<

CompareSignalPdfs: $(INCLUDEDIR)/src/CompareSignalPdfs.C
	$(CXX) $(CXXFLAGS) -o CompareSignalPdfs $(GLIBS) $ $<

makeSPlotsPdfs: $(INCLUDEDIR)/src/ComparesPlotsPdfs.C
	$(CXX) $(CXXFLAGS) -o makeSPlotsPdfs $(GLIBS) $ $<

makeSPlotsClassesPdfs: $(INCLUDEDIR)/src/ComparesPlotsClassesPdfs.C
	$(CXX) $(CXXFLAGS) -o makeSPlotsClassesPdfs $(GLIBS) $ $<

CompareSPlots: $(INCLUDEDIR)/src/ComparesPlotsPdfs.C 
	$(CXX) $(CXXFLAGS) -o CompareSPlots $(GLIBS) $ $<

ComparesPlotsClassesPdfs: $(INCLUDEDIR)/src/ComparesPlotsClassesPdfs.C 
	$(CXX) $(CXXFLAGS) -o ComparesPlotsClassesPdfs $(GLIBS) $ $<

CompareSignalIsolation: $(INCLUDEDIR)/src/CompareSignalIsolation.C
	$(CXX) $(CXXFLAGS) -o CompareSignalIsolation $(GLIBS) $ $<

MakeROC: $(INCLUDEDIR)/src/MakeROC.C
	$(CXX) $(CXXFLAGS) -o MakeROC $(OUTLIBCOMMON)/*o $(GLIBS) $ $<

CompareROC: $(INCLUDEDIR)/src/CompareROC.C \
	$(OUTLIBCOMMON)FiguresOfMeritEvaluator.o 
	$(CXX) $(CXXFLAGS) -o CompareROC $(OUTLIBCOMMON)/*o $(GLIBS) $ $<


EgammaAnalysis.clean:
	rm -f EgammaAnalysis

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f EgammaAnalysis
	rm -f CompareEff
	rm -f CompareMisId
	rm -f CompareClasses
	rm -f produceMCSubtractedFake
	rm -f MakeNotePdfPlots
	rm -f MakeNoteBkgPdfPlots
	rm -f CompareSignalPdfs
	rm -f CompareSignalIsolation
	rm -f ComparesPlotsPdfs
	rm -f makeSPlotsPdfs
	rm -f CompareSPlots
	rm -f ComparesPlotsClassesPdfs
	rm -f MakeNoteIsolationPlots

all:  EgammaAnalysis \
	CompareEff \
	CompareMisId \
	CompareClasses \
	produceMCSubtractedFake \
	MakeNotePdfPlots \
	MakeNoteBkgPdfPlots \
	CompareSignalPdfs \
	CompareSignalIsolation \
	ComparesPlotsClassesPdfs \
	MakeNoteIsolationPlots

