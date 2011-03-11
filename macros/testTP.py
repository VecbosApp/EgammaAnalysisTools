import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.TagProbeFitTreeAnalyzerTrack = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                                      # IO parameters:
                                                      InputFileNames = cms.vstring("results_data/merged.root"),
                                                      InputDirectoryName = cms.string("eleIDdir"),
                                                      InputTreeName = cms.string("T1"),
                                                      OutputFileName = cms.string("effData_WP80.root"),
                                                      #numbrer of CPUs to use for fitting
                                                      NumCPU = cms.uint32(6),
                                                      # specifies wether to save the RooWorkspace containing the data for each bin and
                                                      # the pdf object with the initial and final state snapshots
                                                      SaveWorkspace = cms.bool(True),
                                                      floatShapeParameters = cms.bool(True),
                                                      #fixVars = cms.vstring(""),
                                                      
                                                      # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                                      Variables = cms.PSet( mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                                            pt = cms.vstring("Probe p_{T}", "10", "100", "GeV/c"),
                                                                            eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
                                                                            vertices = cms.vstring("n_{PV}", "0.5", "24.5", "")
                                                                            ),
                                                      # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                                      Categories = cms.PSet( WP80 = cms.vstring("WP80", "dummy[pass=1,fail=0]")
                                                                             ),
                                                      # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
                                                      # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values
                                                      PDFs = cms.PSet( gaussPlusLinear = cms.vstring(
                                                                               "RooExponential::backgroundPass(mass, cPass[-0.055,-3,0.1])",
                                                                               "RooExponential::backgroundFail(mass, cFail[-0.012,-3,0.1])",
                                                                               "EXPR::signalPass('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',pass_mean[91.2, 80.0, 100.0],mass, pass_sigmaL[2.3, 0.5, 10.0],pass_alphaL[0.23],pass_sigmaR[2.3, 0.5, 10.0],pass_alphaR[0.2,0,3])",
                                                                               "EXPR::signalFail('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',pass_mean,mass, fail_sigmaL[2.3, 0.5, 10.0],fail_alphaL[0.23],fail_sigmaR[2.3, 0.5, 10.0],pass_alphaR)",
                                                                                                     "efficiency[0.5,0,1]",
                                                                                                     "signalFractionInPassing[0.9]"
                                                                                                     ),
                                                                       ),

                                                      # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                                      # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                                      Efficiencies = cms.PSet(pt = cms.PSet(EfficiencyCategoryAndState = cms.vstring("WP80","pass"),
                                                                                            UnbinnedVariables = cms.vstring("mass"),
                                                                                            BinnedVariables = cms.PSet(pt = cms.vdouble(10,20,30,40,50,100),
                                                                                                                       ),
                                                                                            BinToPDFmap = cms.vstring("gaussPlusLinear")
                                                                                            ),
 
                                                                               vertices = cms.PSet(EfficiencyCategoryAndState = cms.vstring("WP80","pass"),
                                                                                                UnbinnedVariables = cms.vstring("mass"),
                                                                                                BinnedVariables = cms.PSet(vertices = cms.vdouble(0.5,1.5,2.5,3.5,5.5,7.5,11.5,15.5,24.5),
                                                                                                                            ),
                                                                                                BinToPDFmap = cms.vstring("gaussPlusLinear")
                                                                                                ),

                                                                               eta = cms.PSet( EfficiencyCategoryAndState = cms.vstring("WP80","pass"),
                                                                                               UnbinnedVariables = cms.vstring("mass"),
                                                                                               BinnedVariables = cms.PSet(eta = cms.vdouble(0,0.8,1.4442,1.566,2.0,2.5)
                                                                                                                          ),
                                                                                               BinToPDFmap = cms.vstring("gaussPlusLinear")
                                                                                              ),        
                                                                              )
                                                      )






process.fit = cms.Path( process.TagProbeFitTreeAnalyzerTrack)
