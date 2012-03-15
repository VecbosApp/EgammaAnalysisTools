
#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/FakeElectronSelectorZllPlusOneFake.hh"
#include "CommonTools/include/Counters.hh"
#include "EgammaAnalysisTools/include/SimpleCutsIDSelector.hh"

using namespace bits;
using namespace std;

FakeElectronSelectorZllPlusOneFake::FakeElectronSelectorZllPlusOneFake(TTree *tree)
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
  myCounter.AddVar("twoleptons");
  myCounter.AddVar("zmass");
  myCounter.AddVar("bveto");
  myCounter.AddVar("denom");
}

FakeElectronSelectorZllPlusOneFake::~FakeElectronSelectorZllPlusOneFake() { }


void FakeElectronSelectorZllPlusOneFake::Loop(const char *outname) {

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
  myOutIDTree->addAttributesSignal();
  myOutIDTree->addIsolations();
  myOutIDTree->addRunInfos();
  myOutIDTree->addMore();

  // trigger: electrons - chiara
  cout << "using electrons triggers" << endl;
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  requiredTriggers.push_back("HLT_DoubleMu7_v");
  requiredTriggers.push_back("HLT_Mu13_Mu8_v");
  requiredTriggers.push_back("HLT_Mu17_Mu8_v");
  requiredTriggers.push_back("HLT_Mu17_TkMu8_v");

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
    //if ( _isData && !passedHLT ) continue; // not used: just running on DoubleLepton datasets
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
      
      bool isGoodWEle = pt>20. && 
	aCutSel.outputid(eta,see,dphi,deta,hoe,SimpleCutsIDSelector::kWP90) &&
	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP90);
      
      bool isLooseEle = pt>10. && 
      	aCutSel.outputid(eta,see,dphi,deta,hoe,SimpleCutsIDSelector::kWP95) &&
      	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP95);

      if (isGoodWEle) tightElectrons.push_back(iele);
      if (isLooseEle) looseElectrons.push_back(iele);
    }

    std::vector<int> tightMuons, looseMuons; 
    for(int imu=0; imu<nMuon; imu++) {

      float eta = etaMuon[imu];
      float pt = GetPt(pxMuon[imu],pyMuon[imu]);
      float tkiso = sumPt03Muon[imu]/pt;
      float ecalIsolAbs = 0.0;
      if ( fabs(eta)<1.479 ) ecalIsolAbs = max(0.0,emEt03Muon[imu]-1.0);
      else ecalIsolAbs = emEt03Muon[imu];
      float ecaliso = ecalIsolAbs/pt;
      float hcaliso = hadEt03Muon[imu]/pt;

      // use the same iolation as WP70 for electrons
      bool isGoodWMu = pt>20. && 
	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP90);
      
      bool isLooseMu = pt>10. && 
	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP95);

      if (isGoodWMu) tightMuons.push_back(imu);
      if (isLooseMu) looseMuons.push_back(imu);
    }

    // event selection: at least one lepton pair (1tight x 1loose)
    int zdecay=0; // 0 = mumu
    std::vector<int> zmuons, zelectrons;
    if(tightMuons.size()>0 && looseMuons.size()>1) {
      std::pair<int,int> besttightpair = getBestGoodMuonPair(tightMuons);
      zmuons.push_back(besttightpair.first);
      if(besttightpair.second>-1) zmuons.push_back(besttightpair.second);
      else {
	 std::pair<int,int> bestloosepair = getBestGoodMuonPair(looseMuons);
	 if(bestloosepair.first != besttightpair.first) zmuons.push_back(bestloosepair.first);
	 else zmuons.push_back(bestloosepair.second);
      }
    } else {
      zdecay=1; // 1 = ee
      // exclusive: just to run once and get Zee and Zmm
      if(tightElectrons.size()>0 && looseElectrons.size()>1) {
	std::pair<int,int> besttightpair = getBestGoodElePair(tightElectrons);
	zelectrons.push_back(besttightpair.first);
	if(besttightpair.second>-1) zelectrons.push_back(besttightpair.second);
	else {
	  std::pair<int,int> bestloosepair = getBestGoodElePair(looseElectrons);
	  if(bestloosepair.first != besttightpair.first) zelectrons.push_back(bestloosepair.first);
	  else zelectrons.push_back(bestloosepair.second);
	}
      }
    }


    if(zmuons.size()<2 && zelectrons.size()<2) continue; 
    myCounter.IncrVar("twoleptons",1);    

    std::vector<TLorentzVector> leptons;
    if(zmuons.size()>1) {
      leptons.push_back(TLorentzVector(pxMuon[zmuons[0]],pyMuon[zmuons[0]],pzMuon[zmuons[0]],energyMuon[zmuons[0]]));
      leptons.push_back(TLorentzVector(pxMuon[zmuons[1]],pyMuon[zmuons[1]],pzMuon[zmuons[1]],energyMuon[zmuons[1]]));
    } else if(zelectrons.size()>1) {
      leptons.push_back(TLorentzVector(pxEle[zelectrons[0]],pyEle[zelectrons[0]],pzEle[zelectrons[0]],energyEle[zelectrons[0]]));
      leptons.push_back(TLorentzVector(pxEle[zelectrons[1]],pyEle[zelectrons[1]],pzEle[zelectrons[1]],energyEle[zelectrons[1]]));      
    } else cout << "sanity check: impossibile!" << endl;

    // make the Z
    float mass = (leptons[0]+leptons[1]).M();
    if (mass<75. || mass>106.) continue;
    myCounter.IncrVar("zmass",1);    

    float maxTCHE = -1000.;
    std::vector<int> jets30;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.DeltaPhi(leptons[0].Vect())) < 0.3 || fabs(p3Jet.DeltaPhi(leptons[1].Vect())) < 0.3 ) continue;
      if ( p3Jet.Pt()>30 && fabs(p3Jet.Eta())<3.0 ) jets30.push_back(jet); 
      if ( p3Jet.Pt()>15 && fabs(p3Jet.Eta())<3.0 && trackCountingHighEffBJetTagsAK5Jet[jet]>maxTCHE ) maxTCHE= trackCountingHighEffBJetTagsAK5Jet[jet];
    }

    // veto on max btag in jets with pt > 15 GeV
    if(maxTCHE>2.1) continue;
    myCounter.IncrVar("bveto",1);  

    // look for probes among the reco electrons different from non Z
    std::vector<int> probes;
    for(int iele=0; iele<nEle; iele++) {
      // if Z->ee, probe should not overlap with the Zee electrons
      if(zelectrons.size()>1 && (iele==zelectrons[0] || iele==zelectrons[1])) continue;
      TVector3 probeP3;
      probeP3.SetXYZ(pxEle[iele],pyEle[iele],pzEle[iele]);
      // exclude possible brem candidates in a cone 0.3 from Z leptons
      if(fabs(probeP3.Eta())<2.5 && probeP3.Pt()>5. 
	 && fabs(probeP3.DeltaPhi(leptons[0].Vect()))>0.3
	 && fabs(probeP3.DeltaPhi(leptons[1].Vect()))>0.3 ) probes.push_back(iele);
    }

    // take the highest pt probe
    std::pair<int,int> bestprobes = getBestGoodElePair(probes);     
    int probe(bestprobes.first);
    TLorentzVector probeP4;
    probeP4.SetXYZT(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);

    // at least one probe candidate
    if(probe<0) continue;
    myCounter.IncrVar("denom",1);
    
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    // fill the kine tree - after HLT and denominator
    myOutKineTree -> fill( p3Met.Pt(), -1, mass, leptons[0].Pt(), probeP4.Pt(), 
			   min(fabs(leptons[0].Vect().DeltaPhi(probeP4.Vect())), fabs(leptons[1].Vect().DeltaPhi(probeP4.Vect()))) );
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
      sep = covIEtaIPhiSC[sc]/(sqrt(covIEtaIEtaSC[sc])*sqrt(covIPhiIPhiSC[sc]));
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
        sep = covIEtaIPhiPFSC[sc]/(sqrt(covIEtaIEtaPFSC[sc])*sqrt(covIPhiIPhiPFSC[sc]));
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
    myOutIDTree->fillAttributesSignal(mass,zdecay);
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

// two highest pT muons
std::pair<int,int> FakeElectronSelectorZllPlusOneFake::getBestGoodMuonPair(std::vector<int> goodMuons) {
  
  int theMuon1=-1;
  int theMuon2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;
  for(int iMuon=0;iMuon<goodMuons.size();iMuon++) {
    int muonIndex = goodMuons[iMuon];
    TVector3 pMuon(pxMuon[muonIndex],pyMuon[muonIndex],pzMuon[muonIndex]);
    float thisPt=pMuon.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theMuon2 = theMuon1; theMuon1 = muonIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theMuon2 = muonIndex; }
  }
  return std::make_pair(theMuon1,theMuon2);
}
