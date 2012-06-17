
import FWCore.ParameterSet.Config as cms

import sys
TYPE = sys.argv[2]
FUNC = sys.argv[3]
VARIABLE = sys.argv[4]

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

if TYPE == 'MCSummer12'      : INPUTFILE = "../egammaresults/EleTrees2012_May25/electrons_hlt_summer12mc.root"
elif TYPE == 'DATA_2012'     : INPUTFILE = "../egammaresults/EleTrees2012_May25/electrons_hlt.root"

process.TnP_EleID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(4),
    SaveWorkspace = cms.bool(True),

    InputFileNames = cms.vstring(INPUTFILE),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElEl"),
    OutputFileName = cms.string("hlteff.root"),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60", "120", "GeV/c^{2}"),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV"),        
        eta    = cms.vstring("Probe #eta", "-3.0", "3.0", ""),
        phi    = cms.vstring("Probe #phi", "-3.0", "3.0", "rad"),
        nVtx   = cms.vstring("Number PV", "0", "60", ""),
        #weight = cms.vstring("weight","0","10",""),
        pair_nJet  = cms.vstring("nJets30","0","100",""),
    ),

    Categories = cms.PSet(
        passID  = cms.vstring("passID","dummy[pass=1,fail=0]"),
        passIP  = cms.vstring("passIP","dummy[pass=1,fail=0]"),
        passIso = cms.vstring("passIso","dummy[pass=1,fail=0]"),
        passConvR = cms.vstring("passConvR","dummy[pass=1,fail=0]"),
        HLT_Ele20_Ele4_TnP_Ele4Leg = cms.vstring("HLT_Ele20_Ele4_TnP_Ele4Leg","dummy[pass=1,fail=0]"), # this probe for low pT
        HLT_Ele32_SC17_TnP_SC17Leg = cms.vstring("HLT_Ele32_SC17_TnP_SC17Leg","dummy[pass=1,fail=0]"),
        HLT_Ele17_Ele8_Ele17Leg = cms.vstring("HLT_Ele17_Ele8_Ele17Leg","dummy[pass=1,fail=0]"),
        HLT_Ele17_Ele8_Ele8Leg = cms.vstring("HLT_Ele17_Ele8_Ele8Leg","dummy[pass=1,fail=0]"),
        tag_HLT_Ele20_Ele4_TnP_Ele20Leg = cms.vstring("tag_HLT_Ele20_Ele4_TnP_Ele20Leg","dummy[pass=1,fail=0]"), # this tag for low pT
        tag_HLT_Ele32_SC17_TnP_Ele32Leg = cms.vstring("tag_HLT_Ele32_SC17_TnP_Ele32Leg","dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
         
        CBVoigtianPlusChebychevBackground = cms.vstring(
            "mean[90,80,100]",
            "sigma[3,1,20]",
            "CBShape::cbs(mass, mean, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
            "Voigtian::vs(mass, mean, width[2.495], sigma)",
            "SUM::signal(f1[0.6, 0.4, 0.8]*cbs,f2[0.3, 0.2, 0.8]*vs)",
            "Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
            "Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),                        

        CBVoigtianPlusExponentialBackground = cms.vstring(
                "mean[90,80,100]",
                "sigma[3,1,20]",
                "CBShape::cbs(mass, mean, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
                "Voigtian::vs(mass, mean, width[2.495], sigma)",
                "SUM::signal(f1[0.6, 0.4, 0.8]*cbs,f2[0.3, 0.2, 0.8]*vs)",
                "Exponential::backgroundPass(mass, lp[0,-5,5])",
                "Exponential::backgroundFail(mass, lf[0,-5,5])",
                "efficiency[0.9,0,1]",
                "signalFractionInPassing[0.9]"
        ),
        

        gaussPlusLinear = cms.vstring(

                #"EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[88.5, 80.0, 91.25],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.0, 0.5, 3.0],alphaR[0.1,0,3])",

                "EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[84.5, 83.0, 92.25],mass, sigmaL[5.0, 1.0, 12.0],alphaL[0.2,0.15,0.40],sigmaR[3.0, 2.0, 7.0],alphaR[0.1,0.1,2.0])",

                
                #"RooVoigtian::signal(mass, mean[91.2, 80.0, 100.0],gamma[10,0,100], sigma[2.3, 0.5, 10.0])",

                "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
                "RooExponential::backgroundFail(mass, cFail[0,-10,10])",

                ## "RooPolynomial::backgroundPass(mass, {cPass1[0,-100,100],cPass2[0,-100,100]},2)",
                ## "RooPolynomial::backgroundFail(mass, {cFail1[0,-100,100],cFail2[0,-100,100]},2)",
                "efficiency[0.5,0,1]",
                "signalFractionInPassing[0.9]"
                ),

        #gaussPlusLinear = cms.vstring(
        #        "EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])",
                
                #"RooVoigtian::signal(mass, mean[91.2, 80.0, 100.0],gamma[10,0,100], sigma[2.3, 0.5, 10.0])",
       #         "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
       #         "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
                ## "RooPolynomial::backgroundPass(mass, {cPass1[0,-100,100],cPass2[0,-100,100]},2)",
                ##             "RooPolynomial::backgroundFail(mass, {cFail1[0,-100,100],cFail2[0,-100,100]},2)",
       #         "efficiency[0.5,0,1]",
       #         "signalFractionInPassing[0.9]"
       #         ),

        vpvPlusExpo = cms.vstring(
                "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
                "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
                "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
                "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
                "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
                "efficiency[0.9,0,1]",
                "signalFractionInPassing[0.9]"
                ),

        gpOne = cms.vstring(
                    "Voigtian::signal1Pass(mass, mean1Pass[90,80,100], width[2.495], sigma1Pass[2,1,3])",
                    "Voigtian::signal2Pass(mass, mean2Pass[90,80,100], width,        sigma2Pass[4,2,10])",
                    "SUM::signalPass(vFracPass[0.8,0,1]*signal1Pass, signal2Pass)",
                    "Voigtian::signal1Fail(mass, mean1Fail[90,80,100], width[2.495], sigma1Fail[2,1,3])",
                    "Voigtian::signal2Fail(mass, mean2Fail[90,80,100], width, sigma2Fail[4,2,10])",
                    "SUM::signalFail(vFracFail[0.8,0,1]*signal1Fail, signal2Fail)",
                    "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
                    "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
                    "efficiency[0.9,0,1]",
                    "signalFractionInPassing[0.9]"
                    ),

        gpTwo = cms.vstring(
                    "Voigtian::signal1Pass(mass, mean1[90,80,100],     width[2.495], sigma1[2,1,3])",
                    "Voigtian::signal2Pass(mass, mean2Pass[90,80,100], width,        sigma2Pass[4,2,10])",
                    "SUM::signalPass(vFracPass[0.8,0,1]*signal1Pass, signal2Pass)",
                    "Voigtian::signal2Fail(mass, mean2Fail[90,80,100], width,        sigma2Fail[4,2,10])",
                    "SUM::signalFail(vFracFail[0.8,0,1]*signal1Pass, signal2Fail)",
                    "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
                    "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
                    "efficiency[0.9,0,1]",
                    "signalFractionInPassing[0.9]"
                    ),
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(21),
    saveDistributionsPlot = cms.bool(False),

    #WeightVariable = cms.string("weight"),                                

    Efficiencies = cms.PSet(), # will be filled later
)

PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble( 7, 10, 15, 20, 30, 40, 50, 100 ),
    abseta = cms.vdouble( 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5 ),
    nVtx = cms.vdouble( 0,40 )
)

PT_ETA_BINS_ELE17 = cms.PSet(
    pt     = cms.vdouble(  10, 15, 16, 12.25, 16.5, 16.75, 17, 17.25, 17.5, 18.75, 18.0, 19, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.479, 2.5 )
)
PT_ETA_BINS_ELE8  = cms.PSet(
    pt     = cms.vdouble(  3, 5, 6, 7, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 10, 12, 15, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.479, 2.5 ),
    nVtx = cms.vdouble( 0,40 )
)

VTX_BINS  = cms.PSet(
    pt     = cms.vdouble(  20, 120 ),
    abseta = cms.vdouble(  0.0, 2.4),
    #nVtx = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5)
    nVtx = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,16.5,20.5,25,30,40)
)

OUTPUTFILE = "Electrons_%s_TnP_Z_%s_%s.root" % (FUNC,TYPE,VARIABLE)
process.TnP_EleID.OutputFileName = cms.string(OUTPUTFILE)

numString = cms.vstring()
if "selection" in VARIABLE:
    numString += cms.vstring("passConvR", "pass", "passID", "pass", "passIso", "pass", "passIP", "pass")
    DEN = PT_ETA_BINS.clone(
        #HLT_Ele20_Ele4_TnP_Ele4Leg = cms.vstring("pass"),
        tag_HLT_Ele20_Ele4_TnP_Ele20Leg = cms.vstring("pass")
        )
if "HLT8" in VARIABLE:
    numString += cms.vstring("HLT_Ele17_Ele8_Ele8Leg","pass")
if "HLT17" in VARIABLE:
    numString += cms.vstring("HLT_Ele17_Ele8_Ele17Leg","pass")
if "HLT" in VARIABLE:
    DEN = PT_ETA_BINS.clone(
        passID = cms.vstring("pass"),
        passIP = cms.vstring("pass"),
        passIso = cms.vstring("pass"),
        passConvR = cms.vstring("pass"),
    )

# instead of reweighting, take the core around data
if "MC" in TYPE:
    DEN.nVtx = cms.vdouble(8,17)

process.TnP_EleID.Efficiencies.PASSING_all = cms.PSet(
    EfficiencyCategoryAndState = numString,
    #UnbinnedVariables = cms.vstring("mass","weight"),
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = DEN,
    BinToPDFmap = cms.vstring(FUNC),
    #BinToPDFmap = cms.vstring("gpOne","abseta_bin1__pt_bin0","gaussPlusLinear","*abseta_bin2__pt_bin0*","gaussPlusLinear","*abseta_bin0__pt_bin1*","gaussPlusLinear","*abseta_bin2__pt_bin2*","gaussPlusLinear"),
    )


process.p = cms.Path(process.TnP_EleID)
