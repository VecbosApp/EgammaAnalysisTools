import FWCore.ParameterSet.Config as cms
import sys

WPTOTEST = "newhzzWP"
DOKINE = 1
YEAR = 2011
ISMC = 1

if YEAR == 2011:
    PREFIX = "eff_Data7TeV" if (ISMC==0) else "eff_MC7TeV"
elif YEAR == 2012:
    PREFIX = "eff_Data8TeV" if (ISMC==0) else "eff_MC8TeV"
else:
    PREFIX = "eff_Data"
IDTOTEST=WPTOTEST+"idonly"
ISOTOTEST=WPTOTEST+"isoonly"

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TagProbeFitBase = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                         # IO parameters:
                                         InputFileNames = cms.vstring("../results_data/electrons_hzzidbitsFriend.root") if (ISMC==0) else cms.vstring("../results_data/electrons_zeemc_hzzidbitsFriend.root"),
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
                                                               pt = cms.vstring("Probe p_{T}", "5", "1000", "GeV/c"),
                                                               abseta = cms.vstring("Probe #eta", "0.0", "2.5", ""),
                                                               vertices = cms.vstring("vertices","1","35", ""),
                                                               puW = cms.vstring("puW","0","10","")
                                                               ),
                                         # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                         Categories = cms.PSet( denom = cms.vstring("denom", "dummy[pass=1,fail=0]"),
                                                                wp70trg = cms.vstring("wp70trg", "dummy[pass=1,fail=0]"),
                                                                wp80trg = cms.vstring("wp80trg", "dummy[pass=1,fail=0]"),
                                                                hzzWP = cms.vstring("hzzWP", "dummy[pass=1,fail=0]"),
                                                                hzzWPidonly = cms.vstring("hzzWPidonly", "dummy[pass=1,fail=0]"),
                                                                hzzWPisoonly = cms.vstring("hzzWPisoonly", "dummy[pass=1,fail=0]"),
                                                                newhzzWP = cms.vstring("newhzzWP", "dummy[pass=1,fail=0]"),
                                                                newhzzWPidonly = cms.vstring("newhzzWPidonly", "dummy[pass=1,fail=0]"),
                                                                newhzzWPisoonly = cms.vstring("newhzzWPisoonly", "dummy[pass=1,fail=0]")
                                                                ),
                                         # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
                                         # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
                                         PDFs = cms.PSet( cruijffPlusExpo = cms.vstring("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                                                                                        #"RooVoigtian::signal(mass, mean[91.2, 80.0, 100.0],gamma[10,0,100], sigma[2.3, 0.5, 10.0])",
                                                                                        "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                                                                                        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                                                                                        ##             "RooPolynomial::backgroundPass(mass, {cPass1[0,-100,100],cPass2[0,-100,100]},2)",
                                                                                        ##             "RooPolynomial::backgroundFail(mass, {cFail1[0,-100,100],cFail2[0,-100,100]},2)",
                                                                                        "efficiency[0.5,0,1]",
                                                                                        "signalFractionInPassing[0.9]"
                                                                                        ),
                                                          CBVoigtianPlusChebychevBackground = cms.vstring("mean[90,80,100]",
                                                                                                          "sigma[3,1,20]",
                                                                                                          "CBShape::cbs(mass, mean, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
                                                                                                          "Voigtian::vs(mass, mean, width[2.495], sigma)",
                                                                                                          "SUM::signal(f1[0.6, 0.4, 0.8]*cbs,f2[0.3, 0.2, 0.8]*vs)",
                                                                                                          "Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
                                                                                                          "Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
                                                                                                          "efficiency[0.9,0,1]",
                                                                                                          "signalFractionInPassing[0.9]"
                                                                                                          )
                                                          ),
                                         
                                         binnedFit = cms.bool(True),
                                         binsForFit = cms.uint32(30),

                                         WeightVariable = cms.string("puW"),
                                         
                                         # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                         # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                         Efficiencies = cms.PSet(etapt = cms.PSet( EfficiencyCategoryAndState = cms.vstring(WPTOTEST,"pass",IDTOTEST,"pass",ISOTOTEST,"pass"),
                                                                                   UnbinnedVariables = cms.vstring("mass","puW"),
                                                                                   BinnedVariables = cms.PSet(abseta = cms.vdouble(0.,0.8,1.4442,1.566,2.0,2.5),
                                                                                                              pt = cms.vdouble(5,10,15,20,30,40,50,1000),
                                                                                                              vertices = cms.vdouble(1,35)
                                                                                                              ),
                                                                                   BinToPDFmap = cms.vstring("cruijffPlusExpo")
                                                                                   ),
                                                                 )
                                         )


# the dedicated fits to MC to determine the signal shape. Mostly will be determined from data itself, but the low mass
# tail has to be fixed from MC. Tail (alpha_L) is different for pt < 20 GeV / pt > 20 GeV, but the approximately the same in EB/EE (from MC fits to signal only).
# done for EB/EE pt </> 20 GeV

verticesbinshighpt = cms.vdouble(1,35)
verticesbinslowpt = cms.vdouble(1,35)
if YEAR == 2011:
    verticesbinshighpt = cms.vdouble(1,3,5,7,10,15,20,35)
    verticesbinslowpt = cms.vdouble(1,3,5,7,10,15,20,35)
elif YEAR == 2012:
    verticesbinshighpt = cms.vdouble(1,5,8,12,14,16,18,20,25,30,35)
    verticesbinslowpt = cms.vdouble(1,7,12,17,22,35)
else:
    verticesbinshighpt = cms.vdouble(1,35)
    verticesbinslowpt = cms.vdouble(1,35)

# pt<20 GeV
process.TagProbeFitPt7To10 = process.TagProbeFitBase.clone()
process.TagProbeFitPt7To10.Efficiencies.etapt.BinnedVariables.pt = (7,10)
process.TagProbeFitPt7To10.Efficiencies.etapt.BinnedVariables.abseta = (0.,0.8,1.4442,1.566,2.0,2.5)
process.TagProbeFitPt7To10.OutputFileName = WPTOTEST+"_"+PREFIX+"_Kine_Pt7To10.root"
process.TagProbeFitPt7To10.PDFs.cruijffPlusExpo = cms.vstring("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.40],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                                                               "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                                                               "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                                                               "efficiency[0.5,0,1]",
                                                               "signalFractionInPassing[0.9]"
                                                               )

# pt<20 GeV, barrel: as a function of vertices
process.TagProbeFitVerticesBarrelPt7To10 = process.TagProbeFitPt7To10.clone()
process.TagProbeFitVerticesBarrelPt7To10.Efficiencies.etapt.BinnedVariables.pt = (7,10)
process.TagProbeFitVerticesBarrelPt7To10.Efficiencies.etapt.BinnedVariables.abseta = (0.,1.4442)
process.TagProbeFitVerticesBarrelPt7To10.Efficiencies.etapt.BinnedVariables.vertices = verticesbinslowpt
process.TagProbeFitVerticesBarrelPt7To10.OutputFileName = WPTOTEST+"_"+PREFIX+"_Vertices_barrel_Pt7To10.root"

# pt<20 GeV, endcap: as a function of vertices
process.TagProbeFitVerticesEndcapPt7To10 = process.TagProbeFitPt7To10.clone()
process.TagProbeFitVerticesEndcapPt7To10.Efficiencies.etapt.BinnedVariables.pt = (7,10)
process.TagProbeFitVerticesEndcapPt7To10.Efficiencies.etapt.BinnedVariables.abseta = (1.566,2.5)
process.TagProbeFitVerticesEndcapPt7To10.Efficiencies.etapt.BinnedVariables.vertices = verticesbinslowpt
process.TagProbeFitVerticesEndcapPt7To10.OutputFileName = WPTOTEST+"_"+PREFIX+"_Vertices_endcap_Pt7To10.root"


# pt>20 GeV
process.TagProbeFitPt10To1000 = process.TagProbeFitBase.clone()
process.TagProbeFitPt10To1000.Efficiencies.etapt.BinnedVariables.pt = (10,15,20,30,40,50,1000)
process.TagProbeFitPt10To1000.Efficiencies.etapt.BinnedVariables.abseta = (0.,0.8,1.4442,1.566,2.0,2.5)
process.TagProbeFitPt10To1000.OutputFileName = WPTOTEST+"_"+PREFIX+"_Kine_Pt10To1000.root"
process.TagProbeFitPt10To1000.PDFs.cruijffPlusExpo = cms.vstring("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                                                               "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                                                               "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                                                               "efficiency[0.5,0,1]",
                                                               "signalFractionInPassing[0.9]"
                                                               )

# pt>20 GeV, barrel: as a function of vertices
process.TagProbeFitVerticesBarrelPt10To1000 = process.TagProbeFitPt10To1000.clone()
process.TagProbeFitVerticesBarrelPt10To1000.Efficiencies.etapt.BinnedVariables.pt = (10,1000)
process.TagProbeFitVerticesBarrelPt10To1000.Efficiencies.etapt.BinnedVariables.abseta = (0.,1.4442)
process.TagProbeFitVerticesBarrelPt10To1000.Efficiencies.etapt.BinnedVariables.vertices = verticesbinshighpt
process.TagProbeFitVerticesBarrelPt10To1000.OutputFileName = WPTOTEST+"_"+PREFIX+"_Vertices_barrel_Pt10To1000.root"

# pt>20 GeV, endcap: as a function of vertices
process.TagProbeFitVerticesEndcapPt10To1000 = process.TagProbeFitPt10To1000.clone()
process.TagProbeFitVerticesEndcapPt10To1000.Efficiencies.etapt.BinnedVariables.pt = (10,1000)
process.TagProbeFitVerticesEndcapPt10To1000.Efficiencies.etapt.BinnedVariables.abseta = (1.566,2.5)
process.TagProbeFitVerticesEndcapPt10To1000.Efficiencies.etapt.BinnedVariables.vertices = verticesbinshighpt
process.TagProbeFitVerticesEndcapPt10To1000.OutputFileName = WPTOTEST+"_"+PREFIX+"_Vertices_endcap_Pt10To1000.root"

if DOKINE:
    process.fit = cms.Path(process.TagProbeFitPt7To10 *
                           process.TagProbeFitPt10To1000 )
else:
    process.fit = cms.Path(process.TagProbeFitVerticesBarrelPt7To10 *
                           process.TagProbeFitVerticesEndcapPt7To10 *
                           process.TagProbeFitVerticesBarrelPt10To1000 *
                           process.TagProbeFitVerticesEndcapPt10To1000 )
