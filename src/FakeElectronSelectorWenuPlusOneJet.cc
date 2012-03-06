
#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/FakeElectronSelectorWenuPlusOneJet.hh"
#include "CommonTools/include/Counters.hh"
#include "EgammaAnalysisTools/include/SimpleCutsIDSelector.hh"

using namespace bits;
using namespace std;

FakeElectronSelectorWenuPlusOneJet::FakeElectronSelectorWenuPlusOneJet(TTree *tree)
  : Egamma(tree) {
  
  _isData = false;       // chiara
  
  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1gt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;

  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useHoverE = false;        
  defaultSwitches.m_useOneOverEMinusOneOverP = true;

  LH = new ElectronLikelihood(&(*EB0lt15dir), &(*EB1lt15dir), &(*EElt15dir), &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);

  // configuring the electron BDT
  fMVAHWW = new ElectronIDMVA();
  fMVAHWWNoIP = new ElectronIDMVA();
  fMVAHWW->Initialize("BDTG method",
                      "elebdtweights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml" ,                
                      ElectronIDMVA::kWithIPInfo);

  fMVAHWWNoIP->Initialize("BDTG method",
                          "elebdtweights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml",
                          "elebdtweights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml",
                          "elebdtweights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml",
                          "elebdtweights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml",
                          "elebdtweights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml",
                          "elebdtweights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml" ,                
                      ElectronIDMVA::kNoIPInfo);
  
  // configuring the electron BDT for H->ZZ
  fMVAHZZMC = new ElectronIDMVAHZZ();
  fMVAHZZ = new ElectronIDMVAHZZ();
  fMVAHZZNoIP = new ElectronIDMVAHZZ();
  // Default H->ZZ MC training
  fMVAHZZMC->Initialize("BDTSimpleCat",
                        "elebdtweights/HZZBDT_BDTSimpleCat.weights.xml",
                        ElectronIDMVAHZZ::kBDTSimpleCat);

  // New H->ZZ DATA training, with IP
  fMVAHZZ->Initialize("BDTSimpleCat",
                      "elebdtweights/HZZBDT_BDTSimpleCat_EBSplit_Data.weights.xml",
                      ElectronIDMVAHZZ::kBDTSimpleCatData);

  // New H->ZZ DATA training, no IP
  fMVAHZZNoIP->Initialize("BDTSimpleCat",
                          "elebdtweights/HZZBDT_BDTSimpleCatNoIP_EBSplit_Data.weights.xml",
                          ElectronIDMVAHZZ::kBDTSimpleCatNoIPData);

  // chiara
  // to read good run list
  if (_isData) {
    std::string goodRunGiasoneFile = "config/json/goodCollisions2011.json";
    setJsonGoodRunList(goodRunGiasoneFile); 
    fillRunLSMap();
  }
  
  // counter initialize
  myCounter.SetTitle("EVENT_COUNTER");
  myCounter.AddVar("event");
  myCounter.AddVar("trigger");
  myCounter.AddVar("welectron");
  myCounter.AddVar("hltmatch");
  myCounter.AddVar("zveto");
  myCounter.AddVar("wmt");
  myCounter.AddVar("onejet");
  myCounter.AddVar("bveto");
  myCounter.AddVar("denom");
}

FakeElectronSelectorWenuPlusOneJet::~FakeElectronSelectorWenuPlusOneJet() { }


void FakeElectronSelectorWenuPlusOneJet::Loop(const char *outname) {

  // study vs eta
  float minEta = -2.5;
  float maxEta =  2.5;

  // ---------------------------------------------------------------------
  // to study the event selection
  char filename[200];
  sprintf(filename,"%s_FakeKineTree.root",outname);
  myOutKineTree = new FakeTree(filename);
  sprintf(filename,"%s_FakeIDTree.root",outname);
  myOutIDTree = new RedEleIDTree(filename);
  myOutIDTree->addDenominatorFakeBits();
  myOutIDTree->addIsolations();
  myOutIDTree->addRunInfos();
  myOutIDTree->addMore();

  // trigger: electrons - chiara
  cout << "using electrons triggers" << endl;
  requiredTriggers.push_back("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v");
  requiredTriggers.push_back("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v");
  requiredTriggers.push_back("HLT_Ele25_WP80_PFMT40_v");
  requiredTriggers.push_back("HLT_Ele27_WP80_PFMT50_v");
  requiredTriggers.push_back("HLT_Ele25_WP70_PFMT50_v");
  requiredTriggers.push_back("HLT_Ele32_WP70_PFMT50_v");

  // loop on events
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // good runs selection 
    if (_isData && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_isData && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    // skipping even numbered events (2011A, used for eleID BDT training) - chiara
    // if (runNumber<=173692 && eventNumber%2==0) continue;
    
    // all events
    myCounter.IncrVar("event",1);

    // event selection: trigger
    reloadTriggerMask(true);           
    bool passedHLT = hasPassedHLT();   
    if ( _isData && !passedHLT ) continue;
    myCounter.IncrVar("trigger",1);

    // electrons passing the denominator selection (removed to get it unbiased. Set a bit at the end)
    SimpleCutsIDSelector aCutSel;
    Utils anaUtils;
    std::vector<int> tightElectrons, looseElectrons; 
    for(int iele=0; iele<nEle; iele++) {

      float eta = etaEle[iele];
      float pt = GetPt(pxEle[iele],pyEle[iele]);
      float deta = deltaEtaAtVtxEle[iele];
      float dphi = deltaPhiAtVtxEle[iele];
      float hoe = hOverEEle[iele];
      float see = SigmaiEiE(iele);
      float tkiso = dr03TkSumPtEle[iele]/pt;
      float ecalIsolAbs = 0.0;
      if ( fabs(eta)<1.479 ) ecalIsolAbs = max(0.0,dr03EcalRecHitSumEtEle[iele]-1.0);
      else ecalIsolAbs = dr03EcalRecHitSumEtEle[iele];
      float ecaliso = ecalIsolAbs/pt;
      float hcaliso = dr03HcalTowerSumEtEle[iele]/pt;
      
      bool isGoodWEle = pt>25. && 
	aCutSel.outputid(eta,see,dphi,deta,hoe,SimpleCutsIDSelector::kWP70) &&
	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP70);
      
      bool isLooseEle = aCutSel.outputid(eta,see,dphi,deta,hoe,SimpleCutsIDSelector::kWP95) &&
	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP95);

      if (isGoodWEle) tightElectrons.push_back(iele);
      if (isLooseEle) looseElectrons.push_back(iele);
    }

    // event selection: at least one candidate for W->enu
    if (tightElectrons.size()==0) continue;
    myCounter.IncrVar("welectron",1);    

    // best electron as W candidate
    std::pair<int,int> bestpair = getBestGoodElePair(tightElectrons);     
    int theWEle(bestpair.first);
    TLorentzVector wEleP4;
    wEleP4.SetXYZT(pxEle[theWEle],pyEle[theWEle],pzEle[theWEle],energyEle[theWEle]);

    // match with the HLT firing candidates
    if(_isData) {
      bool HLTmatch = triggerMatch(wEleP4.Eta(),wEleP4.Phi(),0.2);
      if (!HLTmatch) continue;
    }
    myCounter.IncrVar("hltmatch",1);    

    // veto the Z on WP70 x WP95
    bool zmass = false;
    for(int iloose=0;iloose<(int)looseElectrons.size() && iloose!=theWEle && !zmass; ++iloose) {
      TLorentzVector looseP4;
      looseP4.SetXYZT(pxEle[iloose],pyEle[iloose],pzEle[iloose],energyEle[iloose]);
      float mass = (wEleP4+looseP4).M();
      if (mass>60. && mass<120.) zmass = true;  
    }
    
    if(zmass) continue;
    myCounter.IncrVar("zveto",1);    

    if (theWEle==-1) cout << "sanity check: impossibile!" << endl;

    // MT cut
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    float WmT = sqrt(2*wEleP4.Pt()*p3Met.Pt()*(1-cos(wEleP4.Vect().Angle(p3Met))) );      
    if(WmT<50) continue;
    myCounter.IncrVar("wmt",1);  


    float maxTCHE = -1000.;
    std::vector<int> jets30;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.DeltaPhi(wEleP4.Vect())) < 0.3 ) continue;
      if ( p3Jet.Pt()>30 && fabs(p3Jet.Eta())<3.0 ) jets30.push_back(jet); 
      if ( p3Jet.Pt()>15 && fabs(p3Jet.Eta())<3.0 && trackCountingHighEffBJetTagsAK5Jet[jet]>maxTCHE ) maxTCHE= trackCountingHighEffBJetTagsAK5Jet[jet];
    }

    // require <=1 jet with pt > 30 GeV
    if(jets30.size()>1) continue;
    myCounter.IncrVar("onejet",1);

    // veto on max btag in jets with pt > 15 GeV
    if(maxTCHE>2.1) continue;
    myCounter.IncrVar("bveto",1);  

    // look for probes among the reco electrons different from non W
    std::vector<int> probes;
    for(int iele=0; iele<nEle && iele!=theWEle; iele++) {
      TVector3 probeP3;
      probeP3.SetXYZ(pxEle[iele],pyEle[iele],pzEle[iele]);
      // exclude possible brem candidates in a cone 0.3 from W ele
      if(fabs(probeP3.Eta())<2.5 && probeP3.Pt()>5. && fabs(probeP3.DeltaPhi(wEleP4.Vect()))>0.3) probes.push_back(iele);
    }

    // take the highest pt probe
    std::pair<int,int> bestprobes = getBestGoodElePair(probes);     
    int probe(bestprobes.first);
    TLorentzVector probeP4;
    probeP4.SetXYZT(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);

    // at least one probe candidate
    if(probe<0) continue;
    myCounter.IncrVar("denom",1);
    
    // to store the inv mass probe - tag electron
    float mass = (probeP4+wEleP4).M();

    // fill the kine tree - after HLT and denominator
    myOutKineTree -> fill( p3Met.Pt(), WmT, mass, wEleP4.Pt(), probeP4.Pt(), fabs(wEleP4.Vect().DeltaPhi(probeP4.Vect())) );
    myOutKineTree -> store();
    
    // fill the denominator: take only the highest pT denominator candidate
    float etaFake = fabs(probeP4.Eta()); 
    float etFake  = probeP4.Pt();

    // some eleID variables
    float HoE, s1s9, s9s25, phiwidth, etawidth, deta, dphi, fbrem, see, spp, eleopout, eopout, eop, nbrems, recoFlag, EleSCEta;
    float oneoveremoneoverp, eledeta, d0, ip3d, ip3ds, kfhits, kfchi2, e1x5e5x5, dcot, dist;
    float detacalo, dphicalo, sep, dz, gsfchi2, emaxovere, etopovere, ebottomovere, eleftovere, erightovere,
      e2ndovere, e2x5rightovere, e2x5leftovere, e2x5topovere, e2x5bottomovere, 
      e2x5maxovere, e1x5overe, e2x2overe, e3x3overe, e5x5overe, r9,
      EleSCPhi, scenergy, scrawenergy, scesenergy;

    int gsfTrack = gsfTrackIndexEle[probe];
    int kfTrack = trackIndexEle[probe];
    double gsfsign   = (-eleDxyPV(probe,0) >=0 ) ? 1. : -1.;
    int matchConv = (hasMatchedConversionEle[probe]) ? 1 : 0;

    d0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
    dz = eleDzPV(probe,0);
    ip3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
    ip3ds = ip3d/impactPar3DErrorGsfTrack[gsfTrack];
    kfchi2 = (kfTrack>-1) ? trackNormalizedChi2Track[kfTrack] : 0.0;
    kfhits = (kfTrack>-1) ? trackerLayersWithMeasurementTrack[kfTrack] : -1.0;
    gsfchi2 = trackNormalizedChi2GsfTrack[gsfTrack];
    int misshits = expInnerLayersGsfTrack[gsfTrack];
    dcot = convDistEle[probe];
    dist = convDcotEle[probe];
    bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[probe], isEcalDriven);
    int ecalseed;
    HoE = hOverEEle[probe];
    eledeta = deltaEtaEleClusterTrackAtCaloEle[probe];
    deta = deltaEtaAtVtxEle[probe];
    dphi = deltaPhiAtVtxEle[probe];
    detacalo = deltaEtaAtCaloEle[probe];
    dphicalo = deltaPhiAtCaloEle[probe];
    fbrem = fbremEle[probe];
    nbrems = nbremsEle[probe];
    eleopout = eEleClusterOverPoutEle[probe];
    eopout = eSeedOverPoutEle[probe];
    eop = eSuperClusterOverPEle[probe];
    if(ecaldriven) {
      ecalseed = 1;
      int sc = superClusterIndexEle[probe];
      float seedEnergy = seedEnergySC[sc];
      s1s9 = eMaxSC[sc]/eMaxSC[sc];
      s9s25 = e3x3SC[sc]/e5x5SC[sc];
      e1x5e5x5 = (e5x5SC[sc] - e1x5SC[sc])/e5x5SC[sc];
      phiwidth = phiWidthSC[sc];
      etawidth = etaWidthSC[sc];
      see = sqrt(covIEtaIEtaSC[sc]);
      sep = sqrt(covIEtaIPhiSC[sc]);
      spp = sqrt(covIPhiIPhiSC[sc]);
      oneoveremoneoverp = 1./energySC[sc]  - 1./probeP4.Vect().Mag();
      emaxovere = eMaxSC[sc]/seedEnergy;
      etopovere = eTopSC[sc]/seedEnergy;
      ebottomovere = eBottomSC[sc]/seedEnergy;
      eleftovere = eLeftSC[sc]/seedEnergy;
      erightovere = eRightSC[sc]/seedEnergy;
      e2ndovere = e2ndSC[sc]/seedEnergy;
      e2x5rightovere = e2x5RightSC[sc]/seedEnergy;
      e2x5leftovere = e2x5LeftSC[sc]/seedEnergy;
      e2x5topovere = e2x5TopSC[sc]/seedEnergy;
      e2x5bottomovere = e2x5BottomSC[sc]/seedEnergy;
      e2x5maxovere = e2x5MaxSC[sc]/seedEnergy;
      e1x5overe = e1x5SC[sc]/seedEnergy;
      e2x2overe = e2x2SC[sc]/seedEnergy;
      e3x3overe = e3x3SC[sc]/seedEnergy;
      e5x5overe = e5x5SC[sc]/seedEnergy;
      r9 = e3x3SC[sc]/rawEnergySC[sc];
      recoFlag = recoFlagSC[sc];
      EleSCEta = etaSC[sc];
      EleSCPhi = phiSC[sc];            
      scenergy = energySC[sc];
      scrawenergy = rawEnergySC[sc];
      scesenergy = esEnergySC[sc];
    } else {
      ecalseed = 0;
      int sc = PFsuperClusterIndexEle[probe];
      if(sc>-1) {
        float seedEnergy = seedEnergyPFSC[sc];
        s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
        s1s9 = eMaxPFSC[sc]/eMaxPFSC[sc];
        e1x5e5x5 = (e5x5PFSC[sc] - e1x5PFSC[sc])/e5x5PFSC[sc];
        phiwidth = phiWidthPFSC[sc];
        etawidth = etaWidthPFSC[sc];
        see = sqrt(covIEtaIEtaPFSC[sc]);
        sep = sqrt(covIEtaIPhiPFSC[sc]);
        spp = sqrt(covIPhiIPhiPFSC[sc]);
        oneoveremoneoverp = 1./energyPFSC[sc]  - 1./probeP4.Vect().Mag();
        emaxovere = eMaxPFSC[sc]/seedEnergy;
        etopovere = eTopPFSC[sc]/seedEnergy;
        ebottomovere = eBottomPFSC[sc]/seedEnergy;
        eleftovere = eLeftPFSC[sc]/seedEnergy;
        erightovere = eRightPFSC[sc]/seedEnergy;
        e2ndovere = e2ndPFSC[sc]/seedEnergy;
        e2x5rightovere = e2x5RightPFSC[sc]/seedEnergy;
        e2x5leftovere = e2x5LeftPFSC[sc]/seedEnergy;
        e2x5topovere = e2x5TopPFSC[sc]/seedEnergy;
        e2x5bottomovere = e2x5BottomPFSC[sc]/seedEnergy;
        e2x5maxovere = e2x5MaxPFSC[sc]/seedEnergy;
        e1x5overe = e1x5PFSC[sc]/seedEnergy;
        e2x2overe = e2x2PFSC[sc]/seedEnergy;
        e3x3overe = e3x3PFSC[sc]/seedEnergy;
        e5x5overe = e5x5PFSC[sc]/seedEnergy;
        r9 = e3x3PFSC[sc]/rawEnergyPFSC[sc];
        recoFlag = recoFlagPFSC[sc];
        EleSCEta = etaPFSC[sc];
        EleSCPhi = phiPFSC[sc];     
        scenergy = energyPFSC[sc];
        scrawenergy = rawEnergyPFSC[sc];
        scesenergy = esEnergyPFSC[sc];
      } else {
        s9s25 = 999.;
        see = 999.;
        spp = 999.;
      }
    }

    // some MVAs...
    float pfmva = pflowMVAEle[probe];
    float lh=likelihoodRatio(probe,*LH);
    float bdthww = eleBDT(fMVAHWW,probe);
    float bdthwwnoip = eleBDT(fMVAHWWNoIP,probe);
    float bdthzz = eleBDT(fMVAHZZ,probe);
    float bdthzznoip = eleBDT(fMVAHZZNoIP,probe);
    float bdthzzmc = eleBDT(fMVAHZZMC,probe);

    // fill the reduced tree
    float pt = probeP4.Pt();
    float eta = probeP4.Eta();
    int charge = chargeEle[probe];
    myOutIDTree->fillVariables(eleopout,eopout,eop,HoE,deta,dphi,s9s25,s1s9,see,spp,fbrem,
                              nbrems,misshits,dcot,dist,pt,eta,charge,phiwidth,etawidth,
                              oneoveremoneoverp,eledeta,d0,ip3d,ip3ds,kfhits,kfchi2,e1x5e5x5,ecalseed,matchConv);
    myOutIDTree->fillVariables2(detacalo, dphicalo, sep, dz, gsfchi2, emaxovere, etopovere, ebottomovere, eleftovere, erightovere,
                               e2ndovere, e2x5rightovere, e2x5leftovere, e2x5topovere, e2x5bottomovere, 
                               e2x5maxovere, e1x5overe, e2x2overe, e3x3overe, e5x5overe, r9,
                               EleSCPhi, scenergy, scrawenergy, scesenergy);    
    myOutIDTree->fillIsolations(dr03TkSumPtEle[probe] - rhoFastjet*TMath::Pi()*0.3*0.3,
                               dr03EcalRecHitSumEtEle[probe] - rhoFastjet*TMath::Pi()*0.3*0.3,
                               dr03HcalTowerSumEtFullConeEle[probe] - rhoFastjet*TMath::Pi()*0.3*0.3,
                               pfCombinedIsoEle[probe],
                               pfCandChargedIsoEle[probe],pfCandNeutralIsoEle[probe],pfCandPhotonIsoEle[probe]);
    myOutIDTree->fillFakeRateDenomBits(isDenomFake_HwwEgamma(probe),isDenomFake_smurfs(probe));
    myOutIDTree->fillMore(nPV,rhoFastjet,bdthww,bdthzz);
    myOutIDTree->fillRunInfos(runNumber, lumiBlock, eventNumber, nPU, -1);
    myOutIDTree->store();
    
  } // loop events

  // saving the counters
  sprintf(filename,"%sCounters.root",outname);
  myCounter.Save(filename,"recreate");
  myCounter.Draw();
  
  // saving the output tree
  myOutKineTree -> save();
  myOutIDTree   -> save();

}

