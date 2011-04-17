import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TagProbeFitTreeAnalyzerTrack = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                                      # IO parameters:
                                                      InputFileNames = cms.vstring("results_data/merged.root"),
                                                      InputDirectoryName = cms.string("eleIDdir"),
                                                      InputTreeName = cms.string("T1"),
                                                      OutputFileName = cms.string("effData_Eta_WP95.root"),
                                                      #numbrer of CPUs to use for fitting
                                                      NumCPU = cms.uint32(4),
                                                      # specifies wether to save the RooWorkspace containing the data for each bin and
                                                      # the pdf object with the initial and final state snapshots
                                                      SaveWorkspace = cms.bool(True),
                                                      floatShapeParameters = cms.bool(True),
                                                      fixVars = cms.vstring(""),
                                                      
                                                      # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                                      Variables = cms.PSet( mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                                            pt = cms.vstring("Probe p_{T}", "10", "100", "GeV/c"),
                                                                            eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
                                                                            vertices = cms.vstring("vertices","1","20", "")
                                                                            ),
                                                      # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                                      Categories = cms.PSet( WP95 = cms.vstring("WP95", "dummy[pass=1,fail=0]")
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

                                                      # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                                      # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                                      Efficiencies = cms.PSet(eta = cms.PSet( EfficiencyCategoryAndState = cms.vstring("WP95","pass"),
                                                                                              UnbinnedVariables = cms.vstring("mass"),
                                                                                              BinnedVariables = cms.PSet(eta = cms.vdouble(-2.5,-2.3,-1.9,-1.566,-1.4442,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.4442,1.566,1.9,2.3,2.5)
                                                                                                                         ),
                                                                                              BinToPDFmap = cms.vstring("gaussPlusLinear")
                                                                                              ),
                                                                              )
                                                      )



process.fit = cms.Path( process.TagProbeFitTreeAnalyzerTrack)
