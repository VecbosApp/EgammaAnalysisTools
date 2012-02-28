import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TagProbeFitBase = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                         # IO parameters:
                                         InputFileNames = cms.vstring("../results_data/tptree.root"),
                                         InputDirectoryName = cms.string("eleIDdir"),
                                         InputTreeName = cms.string("T1"),
                                         OutputFileName = cms.string("effData.root"),
                                         #numbrer of CPUs to use for fitting
                                         NumCPU = cms.uint32(4),
                                         # specifies wether to save the RooWorkspace containing the data for each bin and
                                         # the pdf object with the initial and final state snapshots
                                         SaveWorkspace = cms.bool(True),
                                         floatShapeParameters = cms.bool(True),
                                         fixVars = cms.vstring(""),
                                         
                                         # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                         Variables = cms.PSet( mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                               pt = cms.vstring("Probe p_{T}", "10", "200", "GeV/c"),
                                                               abseta = cms.vstring("Probe #eta", "0.0", "2.5", ""),
                                                               vertices = cms.vstring("vertices","1","20", "")
                                                               ),
                                         # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                         Categories = cms.PSet( bdthzz = cms.vstring("bdthzz", "dummy[pass=1,fail=0]")
                                                                ),
                                         # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
                                         # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
                                         PDFs = cms.PSet( gaussPlusLinear = cms.vstring("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                                                                                        
                                                                                        #"RooVoigtian::signal(mass, mean[91.2, 80.0, 100.0],gamma[10,0,100], sigma[2.3, 0.5, 10.0])",
                                                                                        "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                                                                                        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                                                                                        ##             "RooPolynomial::backgroundPass(mass, {cPass1[0,-100,100],cPass2[0,-100,100]},2)",
                                                                                        ##             "RooPolynomial::backgroundFail(mass, {cFail1[0,-100,100],cFail2[0,-100,100]},2)",
                                                                                        "efficiency[0.5,0,1]",
                                                                                        "signalFractionInPassing[0.9]"
                                                                                        ),
                                                          ),
                                         
                                         binnedFit = cms.bool(True),
                                         binsForFit = cms.uint32(50),
                                         
                                         # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                         # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                         Efficiencies = cms.PSet(etapt = cms.PSet( EfficiencyCategoryAndState = cms.vstring("bdthzz","pass"),
                                                                                   UnbinnedVariables = cms.vstring("mass"),
                                                                                   BinnedVariables = cms.PSet(abseta = cms.vdouble(0.,1.0,1.4442,1.566,2.0,2.5),
                                                                                                              pt = cms.vdouble(10,15,20,25,50,200)
                                                                                                              ),
                                                                                   BinToPDFmap = cms.vstring("gaussPlusLinear")
                                                                                   ),
                                                                 )
                                         )


# the dedicated fits to MC to determine the signal shape. Mostly will be determined from data itself, but the low mass
# tail has to be fixed from MC. Tail (alpha_L) is different for pt < 20 GeV / pt > 20 GeV, but the approximately the same in EB/EE (from MC fits to signal only).
# done for EB/EE pt </> 20 GeV

# pt<20 GeV
process.TagProbeFitPt10To20 = process.TagProbeFitBase.clone()
process.TagProbeFitPt10To20.Efficiencies.etapt.BinnedVariables.pt = (10,15,20)
process.TagProbeFitPt10To20.Efficiencies.etapt.BinnedVariables.abseta = (0.,1.0,1.4442,1.566,2.0,2.5)
process.TagProbeFitPt10To20.OutputFileName = "efficiency_Data2011_Pt10To20.root"
process.TagProbeFitPt10To20.PDFs.gaussPlusLinear = cms.vstring("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.40],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                                                               "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                                                               "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                                                               "efficiency[0.5,0,1]",
                                                               "signalFractionInPassing[0.9]"
                                                               )

# pt>20 GeV
process.TagProbeFitPt20To200 = process.TagProbeFitBase.clone()
process.TagProbeFitPt20To200.Efficiencies.etapt.BinnedVariables.pt = (20,25,50,200)
process.TagProbeFitPt20To200.Efficiencies.etapt.BinnedVariables.abseta = (0.,1.0,1.4442,1.566,2.0,2.5)
process.TagProbeFitPt20To200.OutputFileName = "efficiency_Data2011_Pt20To200.root"
process.TagProbeFitPt20To200.PDFs.gaussPlusLinear = cms.vstring("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                                                               "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                                                               "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                                                               "efficiency[0.5,0,1]",
                                                               "signalFractionInPassing[0.9]"
                                                               )

process.fit = cms.Path(process.TagProbeFitPt10To20 *
                       process.TagProbeFitPt20To200
                       )

