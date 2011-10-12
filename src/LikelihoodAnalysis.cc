
#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"
#include "CommonTools/include/Counters.hh"

using namespace bits;
using namespace std;

LikelihoodAnalysis::LikelihoodAnalysis(TTree *tree)
  : Egamma(tree) {

  _isData = true;       // chiara
  
  EgammaCutBasedIDWPs.push_back("WP95"); // [0]
  EgammaCutBasedIDWPs.push_back("WP90"); // [1]
  EgammaCutBasedIDWPs.push_back("WP85"); // [2]
  EgammaCutBasedIDWPs.push_back("WP80"); // [3]
  EgammaCutBasedIDWPs.push_back("WP70"); // [4]
  EgammaCutBasedIDWPs.push_back("WP80Smurf"); // [5]
  EgammaCutBasedIDWPs.push_back("WP70Smurf"); // [6]

  EgammaCiCBasedIDWPs.push_back("CiCVeryLoose");
  EgammaCiCBasedIDWPs.push_back("CiCLoose");
  EgammaCiCBasedIDWPs.push_back("CiCMedium");
  EgammaCiCBasedIDWPs.push_back("CiCTight");
  EgammaCiCBasedIDWPs.push_back("CiCSuperTight");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight2");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight3");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight4");

  EgammaLHBasedIDWPs.push_back("LHVeryLoose");  // [0]
  EgammaLHBasedIDWPs.push_back("LHLoose");      // [1] 
  EgammaLHBasedIDWPs.push_back("LHMedium");     // [2] 
  EgammaLHBasedIDWPs.push_back("LHTight");      // [3]
  EgammaLHBasedIDWPs.push_back("LHHyperTight"); // [4]
  EgammaLHBasedIDWPs.push_back("LHTightPFIso"); // [5]
  EgammaLHBasedIDWPs.push_back("LHHyperTightPFIso"); // [6]

  
  // single electron efficiency, simple cuts 
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;
      
      char configDir[50];
      sprintf(configDir,"config/%s",EgammaCutBasedIDWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaCutBasedIDWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
      //  aSelector.ConfigureEcalCleaner("config/");
      EgammaCutBasedID.push_back(aSelector);
    }  

  // single electron efficiency, likelihood
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;

      char configDir[50];
      sprintf(configDir,"config/%s",EgammaLHBasedIDWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaLHBasedIDWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
      //  aSelector.ConfigureEcalCleaner("config/");
      EgammaLHBasedID.push_back(aSelector);
    }  

  // single electron efficiency, cics
  for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
    {
      CiCBasedEleSelector aSelector;

      char configDir[50];
      std::cout << "===== Configuring " <<  EgammaCiCBasedIDWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.Configure(EgammaCiCBasedIDWPs[i],1,1,2);
      //  aSelector.ConfigureEcalCleaner("config/");
      EgammaCiCBasedID.push_back(aSelector);
    }  

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

  // chiara
  // to read good run list
  if (_isData) {
    // std::string goodRunGiasoneFile = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/goodCollisions2011.json";  
    // std::string goodRunGiasoneFile = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/goodCollisions2011_preapproval_184pb.json";
    // std::string goodRunGiasoneFile = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/goodCollisions2011_V7.json";
    // std::string goodRunGiasoneFile = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/goodCollisions2011_815pb.json";
    //// std::string goodRunGiasoneFile    = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/Cert_EPSFINAL_May10ReReco_v2_PromptReco_160404_167913_JSON.txt";
    std::string goodRunGiasoneFile = "config/json/goodCollisions2011.json";
    // std::string goodRunGiasoneFile = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/HWW.conservativeCertificationLP11.json";

    setJsonGoodRunList(goodRunGiasoneFile); 
    fillRunLSMap();
  }

  // counter initialize
  myCounter.SetTitle("EVENT_COUNTER");
  myCounter.AddVar("event");
  myCounter.AddVar("trigger");
  myCounter.AddVar("denom");
  myCounter.AddVar("met");
  myCounter.AddVar("trasvMass");
  myCounter.AddVar("Zmass");
  myCounter.AddVar("leadingExist");
  myCounter.AddVar("deltaPhi");
  myCounter.AddVar("leadingPT");
}



LikelihoodAnalysis::~LikelihoodAnalysis() {
}



float LikelihoodAnalysis::findEquivalentLHCut(float wantEfficiency) {

  int nbins = 500;
  TH1F *LHBinnedHisto = new TH1F("LHBinnedHisto", "LHBinnedHisto", nbins, 0.0, 1.0);


  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    for(int iele=0; iele<nEle; iele++) {
      TVector3 pEle(pxEle[iele],pyEle[iele],pzEle[iele]);
      if(pEle.Pt()>15.0) LHBinnedHisto->Fill(likelihoodRatio(iele,*LH));
    }

  }


  // scan the likelihood to find the cut
  float nEntries = LHBinnedHisto->GetEntries();
  float efficiency = 0.0;
  int bin = nbins+1;

  while (efficiency < wantEfficiency) {

    float integral = LHBinnedHisto->Integral(bin,nbins+1);
    efficiency = integral / nEntries;
    std::cout << "integral = " << integral 
	      << "\tefficiency = " << efficiency
	      << "bin = " << bin << std::endl;
    
    bin--;

  }

  float equivalentLHCut = LHBinnedHisto->GetBinLowEdge(bin);

  std::cout << "Equivalent cut on LH is LH > " << equivalentLHCut << std::endl
	    << "efficiency is: " << efficiency << std::endl
	    << "while wanted is: " << wantEfficiency << std::endl;

  TFile *outfile = TFile::Open("likHisto.root","recreate");
  LHBinnedHisto->Write();
  outfile->Close();

}




void LikelihoodAnalysis::reproduceEgammaCutID() {

  int nEvtGoodElectron = 0;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    int iSelected = -1;
    for(int iele=0; iele<nEle; iele++) {

      for (int icut=0;icut<EgammaCutBasedID.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCutBasedID[icut],iele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if( isEleIDCutBased ) iSelected=iele;
	}
    }

    nEvtGoodElectron++;
    
  }
  
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      EgammaCutBasedID[icut].displayEfficiencies();
    }
    
}



void LikelihoodAnalysis::estimateIDEfficiency(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *GenEta   = new TH1F( "GenEta",  "generated #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaHighPt  = new TH1F( "RecoEtaHighPt", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaLowPt  = new TH1F( "RecoEtaLowPt", "reconstructed #eta", nbinsEta, minEta, maxEta );

  std::vector<TH1F*> CutIdEta;
  std::vector<TH1F*> CutIdOnlyIDEta;
  std::vector<TH1F*> CutIdOnlyIsoEta;
  std::vector<TH1F*> CutIdOnlyConvEta;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEta;
  std::vector<TH1F*> LHIdOnlyIDEta;
  std::vector<TH1F*> LHIdOnlyIsoEta;
  std::vector<TH1F*> LHIdOnlyConvEta;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEta;
  std::vector<TH1F*> CiCIdOnlyIDEta;
  std::vector<TH1F*> CiCIdOnlyIsoEta;
  std::vector<TH1F*> CiCIdOnlyConvEta;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CutIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIDEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> LHIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CutIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEtaLowPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIDEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> LHIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEtaLowPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEtaLowPt.push_back(aHisto);
    }

//   TH1F *CutIdTightCICEta = new TH1F( "CutIdTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *CutIdSuperTightCICEta = new TH1F( "CutIdSuperTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH Loose ID #eta",        nbinsEta, minEta, maxEta );
//   TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH Tight ID #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 18;
  float minPt = 10.0;
  float maxPt = 100.;
  
  TH1F *GenPt   = new TH1F( "GenPt",  "generated p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *GenPtBarrel   = new TH1F( "GenPtBarrel",  "generated p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *GenPtEndcap   = new TH1F( "GenPtEndcap",  "generated p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtBarrel  = new TH1F( "RecoPtBarrel", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtEndcap  = new TH1F( "RecoPtEndcap", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );

  std::vector<TH1F*> CutIdPt;
  std::vector<TH1F*> CutIdOnlyIDPt;
  std::vector<TH1F*> CutIdOnlyIsoPt;
  std::vector<TH1F*> CutIdOnlyConvPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPt;
  std::vector<TH1F*> LHIdOnlyIDPt;
  std::vector<TH1F*> LHIdOnlyIsoPt;
  std::vector<TH1F*> LHIdOnlyConvPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPt;
  std::vector<TH1F*> CiCIdOnlyIDPt;
  std::vector<TH1F*> CiCIdOnlyIsoPt;
  std::vector<TH1F*> CiCIdOnlyConvPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPtBarrel;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPtBarrel;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIDPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CiCIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPtEndcap;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPtEndcap.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPtEndcap;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPtEndcap.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIDPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CiCIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPtEndcap.push_back(aHisto);
    }

//   TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH Loose ID p_{T} (GeV)",        nbinsPt, minPt, maxPt );
//   TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH Tight ID p_{T} (GeV)",        nbinsPt, minPt, maxPt );


  int nbinsPU = 26;
  float minPU = 0;
  float maxPU = 25;
    
  TH1F *GenPU   = new TH1F( "GenPU",  "generated nPU",     nbinsPU, minPU, maxPU );
  TH1F *RecoPU  = new TH1F( "RecoPU", "reconstructed nPU", nbinsPU, minPU, maxPU );
  TH1F *RecoPUHighPt  = new TH1F( "RecoPUHighPt", "reconstructed nPU", nbinsPU, minPU, maxPU );
  TH1F *RecoPULowPt  = new TH1F( "RecoPULowPt", "reconstructed nPU", nbinsPU, minPU, maxPU );

  std::vector<TH1F*> CutIdPU;
  std::vector<TH1F*> CutIdOnlyIDPU;
  std::vector<TH1F*> CutIdOnlyIsoPU;
  std::vector<TH1F*> CutIdOnlyConvPU;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdPU.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIDPU.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIsoPU.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyConvPU.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPU;
  std::vector<TH1F*> LHIdOnlyIDPU;
  std::vector<TH1F*> LHIdOnlyIsoPU;
  std::vector<TH1F*> LHIdOnlyConvPU;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdPU.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIDPU.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIsoPU.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyConvPU.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPU;
  std::vector<TH1F*> CiCIdOnlyIDPU;
  std::vector<TH1F*> CiCIdOnlyIsoPU;
  std::vector<TH1F*> CiCIdOnlyConvPU;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdPU.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIDPU.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIsoPU.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyConvPU.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPUHighPt;
  std::vector<TH1F*> CutIdOnlyIDPUHighPt;
  std::vector<TH1F*> CutIdOnlyIsoPUHighPt;
  std::vector<TH1F*> CutIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIDPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIsoPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyConvPUHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPUHighPt;
  std::vector<TH1F*> LHIdOnlyIDPUHighPt;
  std::vector<TH1F*> LHIdOnlyIsoPUHighPt;
  std::vector<TH1F*> LHIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIDPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIsoPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyConvPUHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPUHighPt;
  std::vector<TH1F*> CiCIdOnlyIDPUHighPt;
  std::vector<TH1F*> CiCIdOnlyIsoPUHighPt;
  std::vector<TH1F*> CiCIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIDPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIsoPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyConvPUHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPULowPt;
  std::vector<TH1F*> CutIdOnlyIDPULowPt;
  std::vector<TH1F*> CutIdOnlyIsoPULowPt;
  std::vector<TH1F*> CutIdOnlyConvPULowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIDPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIsoPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyConvPULowPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPULowPt;
  std::vector<TH1F*> LHIdOnlyIDPULowPt;
  std::vector<TH1F*> LHIdOnlyIsoPULowPt;
  std::vector<TH1F*> LHIdOnlyConvPULowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdPULowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIDPULowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIsoPULowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyConvPULowPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPULowPt;
  std::vector<TH1F*> CiCIdOnlyIDPULowPt;
  std::vector<TH1F*> CiCIdOnlyIsoPULowPt;
  std::vector<TH1F*> CiCIdOnlyConvPULowPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIDPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIsoPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyConvPULowPt.push_back(aHisto);
    }


  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  //  Long64_t nentries = 10000;
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // PU interactions cut
    //    if(nPU>0) continue;

    int mceleindex = -1;
    for(int imc=0; imc<50; imc++) {
      if ( fabs(idMc[imc])==11 && fabs(idMc[mothMc[imc]])==24 ) mceleindex = imc;
    }
    
    if(mceleindex==-1) continue;

    TVector3 mcEle(pMc[mceleindex]*cos(phiMc[mceleindex])*sin(thetaMc[mceleindex]), 
                   pMc[mceleindex]*sin(phiMc[mceleindex])*sin(thetaMc[mceleindex]),
                   pMc[mceleindex]*cos(thetaMc[mceleindex]));

    float mcEta=etaMc[mceleindex];
    float mcPt=pMc[mceleindex] * fabs(sin(thetaMc[mceleindex]));

    // exclude forward electrons
    if(mcPt < minPt || fabs(mcEta) > maxEta) continue;

    GenEta->Fill(mcEta);
    GenPt->Fill(mcPt);

    // loop over ALL reconstructed electrons and find the closest one to the generated one
    float deltaR_min=0.3;
    int matchedRecoEle=-1;
    for(int iele=0; iele<nEle; iele++) {

      TVector3 eleP3(pxEle[iele],pyEle[iele],pzEle[iele]);
      float deltaR = eleP3.DeltaR(mcEle);
      if(deltaR < deltaR_min) {
        matchedRecoEle=iele;
        deltaR_min=deltaR;
      }

    }

    if(matchedRecoEle > -1) {
      Utils anaUtils;
      if (anaUtils.fiducialFlagECAL(fiducialFlagsEle[matchedRecoEle], isEBEEGap))
	continue;

      bool isInEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[matchedRecoEle], isEB);
      bool isInEE = anaUtils.fiducialFlagECAL(fiducialFlagsEle[matchedRecoEle], isEE);
      TVector3 eleP3(pxEle[matchedRecoEle],pyEle[matchedRecoEle],pzEle[matchedRecoEle]);
      bool highPt= (eleP3.Perp()>20.);
      bool lowPt= (eleP3.Perp()<=20.);
      RecoEta->Fill(mcEta);
      RecoPt->Fill(mcPt);
      RecoPU->Fill(nPU);
      if (highPt) {
        RecoEtaHighPt->Fill(mcEta);
        RecoPUHighPt->Fill(nPU);
      }
      if (lowPt) {
        RecoEtaLowPt->Fill(mcEta);
        RecoPULowPt->Fill(nPU);
      }
      if (isInEB) RecoPtBarrel->Fill(mcPt);
      if (isInEE) RecoPtEndcap->Fill(mcPt);

      for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCutBasedID[icut],matchedRecoEle,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

          // for the smurf selection, add further cut for pT<20 GeV. --> move to WP70
          if(TString(EgammaCutBasedIDWPs[icut].c_str()).Contains("Smurf") && eleP3.Perp()<20.) {
            // apply the VBTF70 smurfs
            isEleID(&EgammaCutBasedID[6],matchedRecoEle,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
            isEleIDCutBased = isEleIDCutBased && (fbremEle[matchedRecoEle]>0.15 || (fabs(etaEle[matchedRecoEle])<1.0 && eSuperClusterOverPEle[matchedRecoEle]>0.95));
          }

	  if ( isEleIDCutBased ) {
	    CutIdOnlyIDEta[icut]->Fill(mcEta);
	    CutIdOnlyIDPt[icut]->Fill(mcPt);
            CutIdOnlyIDPU[icut]->Fill(nPU);
	    if (highPt) {
              CutIdOnlyIDEtaHighPt[icut]->Fill(mcEta);
              CutIdOnlyIDPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CutIdOnlyIDEtaLowPt[icut]->Fill(mcEta);
              CutIdOnlyIDPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CutIdOnlyIDPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdOnlyIDPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isIsolCutBased ) {
	    CutIdOnlyIsoEta[icut]->Fill(mcEta);
	    CutIdOnlyIsoPt[icut]->Fill(mcPt);
            CutIdOnlyIsoPU[icut]->Fill(nPU);
	    if (highPt) {
              CutIdOnlyIsoEtaHighPt[icut]->Fill(mcEta);
              CutIdOnlyIsoPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CutIdOnlyIsoEtaLowPt[icut]->Fill(mcEta);
              CutIdOnlyIsoPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CutIdOnlyIsoPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdOnlyIsoPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isConvRejCutBased ) {
	    CutIdOnlyConvEta[icut]->Fill(mcEta);
	    CutIdOnlyConvPt[icut]->Fill(mcPt);
            CutIdOnlyConvPU[icut]->Fill(nPU);
	    if (highPt) {
              CutIdOnlyConvEtaHighPt[icut]->Fill(mcEta);
              CutIdOnlyConvPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CutIdOnlyConvEtaLowPt[icut]->Fill(mcEta);
              CutIdOnlyConvPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CutIdOnlyConvPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdOnlyConvPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CutIdEta[icut]->Fill(mcEta);
	    CutIdPt[icut]->Fill(mcPt);
            CutIdPU[icut]->Fill(nPU);
	    if (highPt) {
              CutIdEtaHighPt[icut]->Fill(mcEta);
              CutIdPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CutIdEtaLowPt[icut]->Fill(mcEta);
              CutIdPULowPt[icut]->Fill(nPU);
            }
            if (isInEB) CutIdPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdPtEndcap[icut]->Fill(mcPt);
	  }
	}
      
//       

//       if ( anaUtils.electronIdVal(eleIdCutsEle[matchedRecoEle], eleIdTightCIC) ) {

//         CutIdTightCICEta->Fill(mcEta);
//         CutIdTightCICPt->Fill(mcPt);
        
//       }

//       if ( anaUtils.electronIdVal(eleIdCutsEle[matchedRecoEle], eleIdSuperTightCIC) ) {

//         CutIdSuperTightCICEta->Fill(mcEta);
//         CutIdSuperTightCICPt->Fill(mcPt);
        
//       }
      
      for (int icut=0;icut<EgammaLHBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaLHBasedID[icut],matchedRecoEle,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    LHIdOnlyIDEta[icut]->Fill(mcEta);
	    LHIdOnlyIDPt[icut]->Fill(mcPt);
            LHIdOnlyIDPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdOnlyIDEtaHighPt[icut]->Fill(mcEta);
              LHIdOnlyIDPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdOnlyIDEtaLowPt[icut]->Fill(mcEta);
              LHIdOnlyIDPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdOnlyIDPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdOnlyIDPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isIsolCutBased ) {
	    LHIdOnlyIsoEta[icut]->Fill(mcEta);
	    LHIdOnlyIsoPt[icut]->Fill(mcPt);
            LHIdOnlyIsoPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdOnlyIsoEtaHighPt[icut]->Fill(mcEta);
              LHIdOnlyIsoPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdOnlyIsoEtaLowPt[icut]->Fill(mcEta);
              LHIdOnlyIsoPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdOnlyIsoPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdOnlyIsoPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isConvRejCutBased ) {
	    LHIdOnlyConvEta[icut]->Fill(mcEta);
	    LHIdOnlyConvPt[icut]->Fill(mcPt);
            LHIdOnlyConvPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdOnlyConvEtaHighPt[icut]->Fill(mcEta);
              LHIdOnlyConvPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdOnlyConvEtaLowPt[icut]->Fill(mcEta);
              LHIdOnlyConvPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdOnlyConvPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdOnlyConvPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    LHIdEta[icut]->Fill(mcEta);
	    LHIdPt[icut]->Fill(mcPt);
            LHIdPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdEtaHighPt[icut]->Fill(mcEta);
              LHIdPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdEtaLowPt[icut]->Fill(mcEta);
              LHIdPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdPtEndcap[icut]->Fill(mcPt);
	  }	  

	}

      for (int icut=0;icut<EgammaCiCBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCiCBasedID[icut],matchedRecoEle,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    CiCIdOnlyIDEta[icut]->Fill(mcEta);
	    CiCIdOnlyIDPt[icut]->Fill(mcPt);
            CiCIdOnlyIDPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdOnlyIDEtaHighPt[icut]->Fill(mcEta);
              CiCIdOnlyIDPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdOnlyIDEtaLowPt[icut]->Fill(mcEta);
              CiCIdOnlyIDPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdOnlyIDPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdOnlyIDPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isIsolCutBased ) {
	    CiCIdOnlyIsoEta[icut]->Fill(mcEta);
	    CiCIdOnlyIsoPt[icut]->Fill(mcPt);
            CiCIdOnlyIsoPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdOnlyIsoEtaHighPt[icut]->Fill(mcEta);
              CiCIdOnlyIsoPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdOnlyIsoEtaLowPt[icut]->Fill(mcEta);
              CiCIdOnlyIsoPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdOnlyIsoPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdOnlyIsoPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isConvRejCutBased ) {
	    CiCIdOnlyConvEta[icut]->Fill(mcEta);
	    CiCIdOnlyConvPt[icut]->Fill(mcPt);
            CiCIdOnlyConvPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdOnlyConvEtaHighPt[icut]->Fill(mcEta);
              CiCIdOnlyConvPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdOnlyConvEtaLowPt[icut]->Fill(mcEta);
              CiCIdOnlyConvPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdOnlyConvPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdOnlyConvPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CiCIdEta[icut]->Fill(mcEta);
	    CiCIdPt[icut]->Fill(mcPt);
            CiCIdPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdEtaHighPt[icut]->Fill(mcEta);
              CiCIdPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdEtaLowPt[icut]->Fill(mcEta);
              CiCIdPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdPtEndcap[icut]->Fill(mcPt);
	  }	  	  

	}
//       

//       if ( lik > minLikelihoodLoose && !eleInGap) {
        
//         LHIdLooseEta->Fill(mcEta);
//         LHIdLoosePt->Fill(mcPt);
        
//       }
      
//       if ( lik > minLikelihoodTight && !eleInGap) {
        
//         LHIdTightEta->Fill(mcEta);
//         LHIdTightPt->Fill(mcPt);
        
//       }
      
    }

  } // loop events

  char filename[200];
  sprintf(filename,"%s-EleEfficiencyEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  //  ElectronEffEta.AddNumerator(GenEta);
  ElectronEffEta.AddNumerator(RecoEta);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(CutIdEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyConvEta[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(LHIdEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyConvEta[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(CiCIdEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyConvEta[icut]);
    }
  ElectronEffEta.SetDenominator(RecoEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("electron efficiency vs #eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("efficiency");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleEfficiencyEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  //  ElectronEffEtaHighPt.AddNumerator(GenEtaHighPt);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(CutIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyConvEtaHighPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(LHIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyConvEtaHighPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(CiCIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyConvEtaHighPt[icut]);
    }
  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("electron efficiency vs #eta");
  ElectronEffEtaHighPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaHighPt.SetYaxisTitle("efficiency");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleEfficiencyEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  //  ElectronEffEtaLowPt.AddNumerator(GenEtaLowPt);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(CutIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyConvEtaLowPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(LHIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyConvEtaLowPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(CiCIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyConvEtaLowPt[icut]);
    }
  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("electron efficiency vs #eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("efficiency");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  sprintf(filename,"%s-EleEfficiencyPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  //  ElectronEffPt.AddNumerator(GenPt);
  ElectronEffPt.AddNumerator(RecoPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(CutIdPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyConvPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(LHIdPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyConvPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(CiCIdPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyConvPt[icut]);
    }

  ElectronEffPt.SetDenominator(RecoPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("electron efficiency vs p_{T}");
  ElectronEffPt.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPt.SetYaxisTitle("efficiency");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleEfficiencyPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  //  ElectronEffPtBarrel.AddNumerator(GenPtBarrel);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(CutIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyConvPtBarrel[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(LHIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyConvPtBarrel[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(CiCIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyConvPtBarrel[icut]);
    }

  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("electron efficiency vs p_{T}");
  ElectronEffPtBarrel.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel.SetYaxisTitle("efficiency");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-EleEfficiencyPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  //  ElectronEffPtEndcap.AddNumerator(GenPtEndcap);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(CutIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyConvPtEndcap[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(LHIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyConvPtEndcap[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(CiCIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyConvPtEndcap[icut]);
    }

  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("electron efficiency vs p_{T}");
  ElectronEffPtEndcap.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap.SetYaxisTitle("efficiency");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();


  sprintf(filename,"%s-EleEfficiencyPU.root",outname);
  EfficiencyEvaluator ElectronEffPU(filename);
  //  ElectronEffPU.AddNumerator(GenPU);
  ElectronEffPU.AddNumerator(RecoPU);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPU.AddNumerator(CutIdPU[icut]);
      ElectronEffPU.AddNumerator(CutIdOnlyIDPU[icut]);
      ElectronEffPU.AddNumerator(CutIdOnlyIsoPU[icut]);
      ElectronEffPU.AddNumerator(CutIdOnlyConvPU[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPU.AddNumerator(LHIdPU[icut]);
      ElectronEffPU.AddNumerator(LHIdOnlyIDPU[icut]);
      ElectronEffPU.AddNumerator(LHIdOnlyIsoPU[icut]);
      ElectronEffPU.AddNumerator(LHIdOnlyConvPU[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPU.AddNumerator(CiCIdPU[icut]);
      ElectronEffPU.AddNumerator(CiCIdOnlyIDPU[icut]);
      ElectronEffPU.AddNumerator(CiCIdOnlyIsoPU[icut]);
      ElectronEffPU.AddNumerator(CiCIdOnlyConvPU[icut]);
    }
  ElectronEffPU.SetDenominator(RecoPU);
  ElectronEffPU.ComputeEfficiencies();
  ElectronEffPU.SetTitle("electron efficiency vs n PU");
  ElectronEffPU.SetXaxisTitle("n PU");
  ElectronEffPU.SetYaxisTitle("efficiency");
  ElectronEffPU.SetYaxisMin(0.0);
  ElectronEffPU.Write();

  sprintf(filename,"%s-EleEfficiencyPUHighPt.root",outname);
  EfficiencyEvaluator ElectronEffPUHighPt(filename);
  //  ElectronEffPUHighPt.AddNumerator(GenPUHighPt);
  ElectronEffPUHighPt.AddNumerator(RecoPUHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPUHighPt.AddNumerator(CutIdPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CutIdOnlyIDPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CutIdOnlyIsoPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CutIdOnlyConvPUHighPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPUHighPt.AddNumerator(LHIdPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(LHIdOnlyIDPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(LHIdOnlyIsoPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(LHIdOnlyConvPUHighPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPUHighPt.AddNumerator(CiCIdPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CiCIdOnlyIDPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CiCIdOnlyIsoPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CiCIdOnlyConvPUHighPt[icut]);
    }
  ElectronEffPUHighPt.SetDenominator(RecoPUHighPt);
  ElectronEffPUHighPt.ComputeEfficiencies();
  ElectronEffPUHighPt.SetTitle("electron efficiency vs n PU");
  ElectronEffPUHighPt.SetXaxisTitle("electron n PU");
  ElectronEffPUHighPt.SetYaxisTitle("efficiency");
  ElectronEffPUHighPt.SetYaxisMin(0.0);
  ElectronEffPUHighPt.Write();

  sprintf(filename,"%s-EleEfficiencyPULowPt.root",outname);
  EfficiencyEvaluator ElectronEffPULowPt(filename);
  //  ElectronEffPULowPt.AddNumerator(GenPULowPt);
  ElectronEffPULowPt.AddNumerator(RecoPULowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPULowPt.AddNumerator(CutIdPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CutIdOnlyIDPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CutIdOnlyIsoPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CutIdOnlyConvPULowPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPULowPt.AddNumerator(LHIdPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(LHIdOnlyIDPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(LHIdOnlyIsoPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(LHIdOnlyConvPULowPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPULowPt.AddNumerator(CiCIdPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CiCIdOnlyIDPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CiCIdOnlyIsoPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CiCIdOnlyConvPULowPt[icut]);
    }
  ElectronEffPULowPt.SetDenominator(RecoPULowPt);
  ElectronEffPULowPt.ComputeEfficiencies();
  ElectronEffPULowPt.SetTitle("electron efficiency vs n PU");
  ElectronEffPULowPt.SetXaxisTitle("electron n PU");
  ElectronEffPULowPt.SetYaxisTitle("efficiency");
  ElectronEffPULowPt.SetYaxisMin(0.0);
  ElectronEffPULowPt.Write();

}



// this is for W+jets
void LikelihoodAnalysis::estimateFakeRate(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;

  // to study which ET cut we want to apply on the jet (for Wjets HWW studies)
  TH1F *FakeJetForThresholdPT = new TH1F("FakeJetForThresholdPT", "fakeable jets p_{T}", 20, 0., 100.);
    
  TH1F *FakeableJetsEta   = new TH1F( "FakeableJetsEta",  "fakeable jets #eta",nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaHighPt  = new TH1F( "RecoEtaHighPt", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaLowPt  = new TH1F( "RecoEtaLowPt", "reconstructed #eta", nbinsEta, minEta, maxEta );
  std::vector<TH1F*> CutIdEta;
  std::vector<TH1F*> CutIdOnlyIDEta;
  std::vector<TH1F*> CutIdOnlyIsoEta;
  std::vector<TH1F*> CutIdOnlyConvEta;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEta;
  std::vector<TH1F*> LHIdOnlyIDEta;
  std::vector<TH1F*> LHIdOnlyIsoEta;
  std::vector<TH1F*> LHIdOnlyConvEta;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEta;
  std::vector<TH1F*> CiCIdOnlyIDEta;
  std::vector<TH1F*> CiCIdOnlyIsoEta;
  std::vector<TH1F*> CiCIdOnlyConvEta;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CutIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIDEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> LHIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CutIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEtaLowPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIDEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> LHIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEtaLowPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEtaLowPt.push_back(aHisto);
    }

//   TH1F *CutIdWP80Eta = new TH1F( "CutIdWP80Eta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *CutIdTightCICEta = new TH1F( "CutIdTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *CutIdSuperTightCICEta = new TH1F( "CutIdSuperTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH ID tight #eta",        nbinsEta, minEta, maxEta );
//   TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH ID loose #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 18;
  float minPt = 10.0;
  float maxPt = 100.;
  
  TH1F *FakeableJetsPt   = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtBarrel  = new TH1F( "RecoPtBarrel", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtEndcap  = new TH1F( "RecoPtEndcap", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  std::vector<TH1F*> CutIdPt;
  std::vector<TH1F*> CutIdOnlyIDPt;
  std::vector<TH1F*> CutIdOnlyIsoPt;
  std::vector<TH1F*> CutIdOnlyConvPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPt;
  std::vector<TH1F*> LHIdOnlyIDPt;
  std::vector<TH1F*> LHIdOnlyIsoPt;
  std::vector<TH1F*> LHIdOnlyConvPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPt;
  std::vector<TH1F*> CiCIdOnlyIDPt;
  std::vector<TH1F*> CiCIdOnlyIsoPt;
  std::vector<TH1F*> CiCIdOnlyConvPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPtBarrel;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPtBarrel;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIDPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CiCIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPtEndcap;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPtEndcap.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPtEndcap;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPtEndcap.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIDPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CiCIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPtEndcap.push_back(aHisto);
    }

//   TH1F *CutIdWP80Pt = new TH1F( "CutIdWP80Pt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH ID tight p_{T} (GeV)",        nbinsPt, minPt, maxPt );
//   TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH ID loose p_{T} (GeV)",        nbinsPt, minPt, maxPt );


  int nbinsPU = 26;
  float minPU = 0;
  float maxPU = 25;
    
  TH1F *GenPU   = new TH1F( "GenPU",  "generated nPU",     nbinsPU, minPU, maxPU );
  TH1F *RecoPU  = new TH1F( "RecoPU", "reconstructed nPU", nbinsPU, minPU, maxPU );
  TH1F *RecoPUHighPt  = new TH1F( "RecoPUHighPt", "reconstructed nPU", nbinsPU, minPU, maxPU );
  TH1F *RecoPULowPt  = new TH1F( "RecoPULowPt", "reconstructed nPU", nbinsPU, minPU, maxPU );

  std::vector<TH1F*> CutIdPU;
  std::vector<TH1F*> CutIdOnlyIDPU;
  std::vector<TH1F*> CutIdOnlyIsoPU;
  std::vector<TH1F*> CutIdOnlyConvPU;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdPU.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIDPU.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIsoPU.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyConvPU.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPU;
  std::vector<TH1F*> LHIdOnlyIDPU;
  std::vector<TH1F*> LHIdOnlyIsoPU;
  std::vector<TH1F*> LHIdOnlyConvPU;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdPU.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIDPU.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIsoPU.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyConvPU.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPU;
  std::vector<TH1F*> CiCIdOnlyIDPU;
  std::vector<TH1F*> CiCIdOnlyIsoPU;
  std::vector<TH1F*> CiCIdOnlyConvPU;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdPU.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIDPU.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIsoPU.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyConvPU.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPUHighPt;
  std::vector<TH1F*> CutIdOnlyIDPUHighPt;
  std::vector<TH1F*> CutIdOnlyIsoPUHighPt;
  std::vector<TH1F*> CutIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIDPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIsoPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyConvPUHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPUHighPt;
  std::vector<TH1F*> LHIdOnlyIDPUHighPt;
  std::vector<TH1F*> LHIdOnlyIsoPUHighPt;
  std::vector<TH1F*> LHIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIDPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIsoPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyConvPUHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPUHighPt;
  std::vector<TH1F*> CiCIdOnlyIDPUHighPt;
  std::vector<TH1F*> CiCIdOnlyIsoPUHighPt;
  std::vector<TH1F*> CiCIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIDPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIsoPUHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyConvPUHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPULowPt;
  std::vector<TH1F*> CutIdOnlyIDPULowPt;
  std::vector<TH1F*> CutIdOnlyIsoPULowPt;
  std::vector<TH1F*> CutIdOnlyConvPULowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIDPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyIsoPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CutIdOnlyConvPULowPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPULowPt;
  std::vector<TH1F*> LHIdOnlyIDPULowPt;
  std::vector<TH1F*> LHIdOnlyIsoPULowPt;
  std::vector<TH1F*> LHIdOnlyConvPULowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdPULowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIDPULowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyIsoPULowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      LHIdOnlyConvPULowPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPULowPt;
  std::vector<TH1F*> CiCIdOnlyIDPULowPt;
  std::vector<TH1F*> CiCIdOnlyIsoPULowPt;
  std::vector<TH1F*> CiCIdOnlyConvPULowPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIDPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyIsoPULowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
      CiCIdOnlyConvPULowPt.push_back(aHisto);
    }

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // PU interactions cut
    //    if(nPU>0) continue;

    bool tauPresence=false;

    int mceleindex = -1;
    for(int imc=0; imc<50; imc++) {
      // not only ele from W->enu: there is enu emission (V_ud?) in madgraph
      if ( (fabs(idMc[imc])==11) ) {
        mceleindex=imc;
        break;
      }
      // since the stable particle list is truncated, if there is a tau 
      // not possible to say what happens...
      if ( (fabs(idMc[imc])==15) ) {
        tauPresence=true;
        break;
      }
    }

    if(tauPresence) continue;

    //debug
//     if(mceleindex==-1) {
//       for(int imc=0; imc<20; imc++) {
//         cout << "imc = " << imc << "\tidMc = " << idMc[imc] << "\tmothMc = " << mothMc[imc] << endl;
//       }
//     }
    //enddebug

    TVector3 mcEle(0,0,0);
    if(mceleindex>-1) mcEle = TVector3(pMc[mceleindex]*cos(phiMc[mceleindex])*sin(thetaMc[mceleindex]), 
                                       pMc[mceleindex]*sin(phiMc[mceleindex])*sin(thetaMc[mceleindex]),
                                       pMc[mceleindex]*cos(thetaMc[mceleindex]));
    
//     int nJetsReco=0;
//     int nEleReco=0;

    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 10.0 ) {
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Jet.DeltaR(mcEle);
        // remove from denominator the electron reconstructed as jet
        if( deltaR>0.5 ) { 
          FakeableJetsEta->Fill( p3Jet.Eta() );
          FakeableJetsPt->Fill( p3Jet.Pt() );
          //          nJetsReco++;
        }
      }
    }

    //  consider the electrons (numerator) 
    for ( int ele=0; ele<nEle; ele++ ) {
      TVector3 p3Ele(pxEle[ele], pyEle[ele], pzEle[ele]);
      if ( fabs(etaEle[ele]) < 2.5 && p3Ele.Pt() > 10.0 ) {
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Ele.DeltaR(mcEle);
        if(deltaR<0.3) continue;   

        float dREleJet_min = 1000;
        int closestJet=-1;

        for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
          TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
          if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 10.0 ) {          
            float dREleJet = p3Jet.DeltaR(p3Ele);
            if(dREleJet<dREleJet_min) {
              closestJet=jet;
              dREleJet_min=dREleJet;
            }
          }
        }

        //        nEleReco++;

        if(closestJet > -1) {

          TVector3 p3ClosestJet(pxAK5PFPUcorrJet[closestJet],pyAK5PFPUcorrJet[closestJet],pzAK5PFPUcorrJet[closestJet]);

	  //!!//
//           float etFake=p3ClosestJet.Pt();
//           float etaFake=p3ClosestJet.Eta();

          float etFake=p3Ele.Pt();
          float etaFake=etaEle[ele];

	  // to study the dependence on the jet threshold for Wjets in HWW
          if ( etFake > 10.0 ) FakeJetForThresholdPT ->Fill( p3ClosestJet.Pt() ); // hardcoded 
	  
          Utils anaUtils;
	  bool isInEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEB);
	  bool isInEE = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEE);
	  bool highPt= (etFake>20.);
	  bool lowPt= (etFake<=20.);
	  RecoEta->Fill(etaFake);
	  RecoPt->Fill(etFake);
          RecoPU->Fill(nPU);
	  if (highPt) {
            RecoEtaHighPt->Fill(etaFake);
            RecoPUHighPt->Fill(nPU);
          }
	  if (lowPt) {
            RecoEtaLowPt->Fill(etaFake);
            RecoPULowPt->Fill(nPU);
          }
	  if (isInEB) RecoPtBarrel->Fill(etFake);
	  if (isInEE) RecoPtEndcap->Fill(etFake);

          for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
            {
              bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
              isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
              isEleID(&EgammaCutBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

              // for the smurf selection, add further cut for pT<20 GeV. --> move to WP70
              if(TString(EgammaCutBasedIDWPs[icut].c_str()).Contains("Smurf") && etFake<20.) {
                // apply the VBTF70 smurfs
                isEleID(&EgammaCutBasedID[6],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
                isEleIDCutBased = isEleIDCutBased && (fbremEle[ele]>0.15 || (fabs(etaEle[ele])<1.0 && eSuperClusterOverPEle[ele]>0.95));
              }

              if ( isEleIDCutBased ) {
                CutIdOnlyIDEta[icut]->Fill(etaFake);
                CutIdOnlyIDPt[icut]->Fill(etFake);
                CutIdOnlyIDPU[icut]->Fill(nPU);
                if (highPt) {
                  CutIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
                  CutIdOnlyIDPUHighPt[icut]->Fill(nPU);
                }
                if (lowPt) {
                  CutIdOnlyIDEtaLowPt[icut]->Fill(etaFake);
                  CutIdOnlyIDPULowPt[icut]->Fill(nPU);
                }
                if (isInEB) CutIdOnlyIDPtBarrel[icut]->Fill(etFake);
                if (isInEE) CutIdOnlyIDPtEndcap[icut]->Fill(etFake);
              }
              if ( isIsolCutBased ) {
                CutIdOnlyIsoEta[icut]->Fill(etaFake);
                CutIdOnlyIsoPt[icut]->Fill(etFake);
                CutIdOnlyIsoPU[icut]->Fill(nPU);
                if (highPt) {
                  CutIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
                  CutIdOnlyIsoPUHighPt[icut]->Fill(nPU);
                }
                if (lowPt) {
                  CutIdOnlyIsoEtaLowPt[icut]->Fill(etaFake);
                  CutIdOnlyIsoPULowPt[icut]->Fill(nPU);
                }
                if (isInEB) CutIdOnlyIsoPtBarrel[icut]->Fill(etFake);
                if (isInEE) CutIdOnlyIsoPtEndcap[icut]->Fill(etFake);
              }
              if ( isConvRejCutBased ) {
                CutIdOnlyConvEta[icut]->Fill(etaFake);
                CutIdOnlyConvPt[icut]->Fill(etFake);
                CutIdOnlyConvPU[icut]->Fill(nPU);
                if (highPt) {
                  CutIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
                  CutIdOnlyConvPUHighPt[icut]->Fill(nPU);
                }
                if (lowPt) {
                  CutIdOnlyConvEtaLowPt[icut]->Fill(etaFake);
                  CutIdOnlyConvPULowPt[icut]->Fill(nPU);
                }
                if (isInEB) CutIdOnlyConvPtBarrel[icut]->Fill(etFake);
                if (isInEE) CutIdOnlyConvPtEndcap[icut]->Fill(etFake);
              }
              if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
                CutIdEta[icut]->Fill(etaFake);
                CutIdPt[icut]->Fill(etFake);
                CutIdPU[icut]->Fill(nPU);
                if (highPt) {
                  CutIdEtaHighPt[icut]->Fill(etaFake);
                  CutIdPUHighPt[icut]->Fill(nPU);
                }
                if (lowPt) {
                  CutIdEtaLowPt[icut]->Fill(etaFake);
                  CutIdPULowPt[icut]->Fill(nPU);
                }
                if (isInEB) CutIdPtBarrel[icut]->Fill(etFake);
                if (isInEE) CutIdPtEndcap[icut]->Fill(etFake);
              }
            }

      for (int icut=0;icut<EgammaLHBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaLHBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    LHIdOnlyIDEta[icut]->Fill(etaFake);
	    LHIdOnlyIDPt[icut]->Fill(etFake);
            LHIdOnlyIDPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
              LHIdOnlyIDPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdOnlyIDEtaLowPt[icut]->Fill(etaFake);
              LHIdOnlyIDPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdOnlyIDPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdOnlyIDPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isIsolCutBased ) {
	    LHIdOnlyIsoEta[icut]->Fill(etaFake);
	    LHIdOnlyIsoPt[icut]->Fill(etFake);
            LHIdOnlyIsoPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
              LHIdOnlyIsoPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdOnlyIsoEtaLowPt[icut]->Fill(etaFake);
              LHIdOnlyIsoPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdOnlyIsoPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdOnlyIsoPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isConvRejCutBased ) {
	    LHIdOnlyConvEta[icut]->Fill(etaFake);
	    LHIdOnlyConvPt[icut]->Fill(etFake);
            LHIdOnlyConvPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
              LHIdOnlyConvPUHighPt[icut]->Fill(nPU);
            }
            if (lowPt) {
              LHIdOnlyConvEtaLowPt[icut]->Fill(etaFake);
              LHIdOnlyConvPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdOnlyConvPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdOnlyConvPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    LHIdEta[icut]->Fill(etaFake);
	    LHIdPt[icut]->Fill(etFake);
            LHIdPU[icut]->Fill(nPU);
	    if (highPt) {
              LHIdEtaHighPt[icut]->Fill(etaFake);
              LHIdPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              LHIdEtaLowPt[icut]->Fill(etaFake);
              LHIdPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) LHIdPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdPtEndcap[icut]->Fill(etFake);
	  }	  

	}

      for (int icut=0;icut<EgammaCiCBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCiCBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    CiCIdOnlyIDEta[icut]->Fill(etaFake);
	    CiCIdOnlyIDPt[icut]->Fill(etFake);
            CiCIdOnlyIDPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
              CiCIdOnlyIDPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdOnlyIDEtaLowPt[icut]->Fill(etaFake);
              CiCIdOnlyIDPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdOnlyIDPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdOnlyIDPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isIsolCutBased ) {
	    CiCIdOnlyIsoEta[icut]->Fill(etaFake);
	    CiCIdOnlyIsoPt[icut]->Fill(etFake);
            CiCIdOnlyIsoPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
              CiCIdOnlyIsoPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdOnlyIsoEtaLowPt[icut]->Fill(etaFake);
              CiCIdOnlyIsoPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdOnlyIsoPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdOnlyIsoPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isConvRejCutBased ) {
	    CiCIdOnlyConvEta[icut]->Fill(etaFake);
	    CiCIdOnlyConvPt[icut]->Fill(etFake);
            CiCIdOnlyConvPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
              CiCIdOnlyConvPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdOnlyConvEtaLowPt[icut]->Fill(etaFake);
              CiCIdOnlyConvPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdOnlyConvPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdOnlyConvPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CiCIdEta[icut]->Fill(etaFake);
	    CiCIdPt[icut]->Fill(etFake);
            CiCIdPU[icut]->Fill(nPU);
	    if (highPt) {
              CiCIdEtaHighPt[icut]->Fill(etaFake);
              CiCIdPUHighPt[icut]->Fill(nPU);
            }
	    if (lowPt) {
              CiCIdEtaLowPt[icut]->Fill(etaFake);
              CiCIdPULowPt[icut]->Fill(nPU);
            }
	    if (isInEB) CiCIdPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdPtEndcap[icut]->Fill(etFake);
	  }	  	  

	}
	
          
//           Utils anaUtils;
          
//           if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdTightCIC) ) {
            
//             CutIdTightCICEta->Fill(etaFake);
//             CutIdTightCICPt->Fill(etFake);
            
//           }

//           if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdSuperTightCIC) ) {
            
//             CutIdSuperTightCICEta->Fill(etaFake);
//             CutIdSuperTightCICPt->Fill(etFake);
            
//           }

//           float lik = 
//           bool eleInGap = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEBEEGap);

//           if ( lik > minLikelihoodLoose && !eleInGap) {

//             LHIdLooseEta->Fill(etaFake);
//             LHIdLoosePt->Fill(etFake);

//           }

//           if ( lik > minLikelihoodTight && !eleInGap) {

//             LHIdTightEta->Fill(etaFake);
//             LHIdTightPt->Fill(etFake);
	  
// 	}
	
	}	
      } // electron acceptance & pt cut
      
    } // loop ele
  
    //    cout << "nEleReco = " << nEleReco << " nJetsReco = " << nJetsReco << endl;
  
  }// loop events
  
  
  char filename[200];
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  //  ElectronEffEta.AddNumerator(GenEta);
  ElectronEffEta.AddNumerator(RecoEta);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(CutIdEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyConvEta[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(LHIdEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyConvEta[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(CiCIdEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyConvEta[icut]);
    }
  ElectronEffEta.SetDenominator(RecoEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("fake rate vs #eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("Fake rate");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleMisidEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  //  ElectronEffEtaHighPt.AddNumerator(GenEtaHighPt);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(CutIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyConvEtaHighPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(LHIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyConvEtaHighPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(CiCIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyConvEtaHighPt[icut]);
    }
  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaHighPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  //  ElectronEffEtaLowPt.AddNumerator(GenEtaLowPt);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(CutIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyConvEtaLowPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(LHIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyConvEtaLowPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(CiCIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyConvEtaLowPt[icut]);
    }
  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  //  ElectronEffPt.AddNumerator(GenPt);
  ElectronEffPt.AddNumerator(RecoPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(CutIdPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyConvPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(LHIdPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyConvPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(CiCIdPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyConvPt[icut]);
    }

  ElectronEffPt.SetDenominator(RecoPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("fake rate vs p_{T}");
  ElectronEffPt.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPt.SetYaxisTitle("Fake rate");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleMisidPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  //  ElectronEffPtBarrel.AddNumerator(GenPtBarrel);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(CutIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyConvPtBarrel[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(LHIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyConvPtBarrel[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(CiCIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyConvPtBarrel[icut]);
    }

  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("fake rate vs p_{T}");
  ElectronEffPtBarrel.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-EleMisidPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  //  ElectronEffPtEndcap.AddNumerator(GenPtEndcap);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(CutIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyConvPtEndcap[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(LHIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyConvPtEndcap[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(CiCIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyConvPtEndcap[icut]);
    }

  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("fake rate vs p_{T}");
  ElectronEffPtEndcap.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();


  sprintf(filename,"%s-EleMisidPU.root",outname);
  EfficiencyEvaluator ElectronEffPU(filename);
  //  ElectronEffPU.AddNumerator(GenPU);
  ElectronEffPU.AddNumerator(RecoPU);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPU.AddNumerator(CutIdPU[icut]);
      ElectronEffPU.AddNumerator(CutIdOnlyIDPU[icut]);
      ElectronEffPU.AddNumerator(CutIdOnlyIsoPU[icut]);
      ElectronEffPU.AddNumerator(CutIdOnlyConvPU[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPU.AddNumerator(LHIdPU[icut]);
      ElectronEffPU.AddNumerator(LHIdOnlyIDPU[icut]);
      ElectronEffPU.AddNumerator(LHIdOnlyIsoPU[icut]);
      ElectronEffPU.AddNumerator(LHIdOnlyConvPU[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPU.AddNumerator(CiCIdPU[icut]);
      ElectronEffPU.AddNumerator(CiCIdOnlyIDPU[icut]);
      ElectronEffPU.AddNumerator(CiCIdOnlyIsoPU[icut]);
      ElectronEffPU.AddNumerator(CiCIdOnlyConvPU[icut]);
    }
  ElectronEffPU.SetDenominator(RecoPU);
  ElectronEffPU.ComputeEfficiencies();
  ElectronEffPU.SetTitle("fake rate vs n PU");
  ElectronEffPU.SetXaxisTitle("n PU");
  ElectronEffPU.SetYaxisTitle("Fake rate");
  ElectronEffPU.SetYaxisMin(0.0);
  ElectronEffPU.Write();

  sprintf(filename,"%s-EleMisidPUHighPt.root",outname);
  EfficiencyEvaluator ElectronEffPUHighPt(filename);
  //  ElectronEffPUHighPt.AddNumerator(GenPUHighPt);
  ElectronEffPUHighPt.AddNumerator(RecoPUHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPUHighPt.AddNumerator(CutIdPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CutIdOnlyIDPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CutIdOnlyIsoPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CutIdOnlyConvPUHighPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPUHighPt.AddNumerator(LHIdPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(LHIdOnlyIDPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(LHIdOnlyIsoPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(LHIdOnlyConvPUHighPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPUHighPt.AddNumerator(CiCIdPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CiCIdOnlyIDPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CiCIdOnlyIsoPUHighPt[icut]);
      ElectronEffPUHighPt.AddNumerator(CiCIdOnlyConvPUHighPt[icut]);
    }
  ElectronEffPUHighPt.SetDenominator(RecoPUHighPt);
  ElectronEffPUHighPt.ComputeEfficiencies();
  ElectronEffPUHighPt.SetTitle("fake rate vs n PU");
  ElectronEffPUHighPt.SetXaxisTitle("electron n PU");
  ElectronEffPUHighPt.SetYaxisTitle("Fake rate");
  ElectronEffPUHighPt.SetYaxisMin(0.0);
  ElectronEffPUHighPt.Write();

  sprintf(filename,"%s-EleMisidPULowPt.root",outname);
  EfficiencyEvaluator ElectronEffPULowPt(filename);
  //  ElectronEffPULowPt.AddNumerator(GenPULowPt);
  ElectronEffPULowPt.AddNumerator(RecoPULowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPULowPt.AddNumerator(CutIdPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CutIdOnlyIDPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CutIdOnlyIsoPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CutIdOnlyConvPULowPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPULowPt.AddNumerator(LHIdPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(LHIdOnlyIDPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(LHIdOnlyIsoPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(LHIdOnlyConvPULowPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPULowPt.AddNumerator(CiCIdPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CiCIdOnlyIDPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CiCIdOnlyIsoPULowPt[icut]);
      ElectronEffPULowPt.AddNumerator(CiCIdOnlyConvPULowPt[icut]);
    }
  ElectronEffPULowPt.SetDenominator(RecoPULowPt);
  ElectronEffPULowPt.ComputeEfficiencies();
  ElectronEffPULowPt.SetTitle("fake rate vs n PU");
  ElectronEffPULowPt.SetXaxisTitle("nPU");
  ElectronEffPULowPt.SetYaxisTitle("Fake rate");
  ElectronEffPULowPt.SetYaxisMin(0.0);
  ElectronEffPULowPt.Write();

  sprintf(filename,"%s-EleMisidPt-Threshold.root",outname);
  TFile fileThreshold(filename,"RECREATE");
  fileThreshold.cd();
  FakeJetForThresholdPT->Write();
}

void LikelihoodAnalysis::estimateFakeRateQCD(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *FakeableJetsEta   = new TH1F( "FakeableJetsEta",  "fakeable jets #eta",     nbinsEta, minEta, maxEta );

  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaHighPt  = new TH1F( "RecoEtaHighPt", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaLowPt  = new TH1F( "RecoEtaLowPt", "reconstructed #eta", nbinsEta, minEta, maxEta );

  std::vector<TH1F*> CutIdEta;
  std::vector<TH1F*> CutIdOnlyIDEta;
  std::vector<TH1F*> CutIdOnlyIsoEta;
  std::vector<TH1F*> CutIdOnlyConvEta;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEta;
  std::vector<TH1F*> LHIdOnlyIDEta;
  std::vector<TH1F*> LHIdOnlyIsoEta;
  std::vector<TH1F*> LHIdOnlyConvEta;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEta;
  std::vector<TH1F*> CiCIdOnlyIDEta;
  std::vector<TH1F*> CiCIdOnlyIsoEta;
  std::vector<TH1F*> CiCIdOnlyConvEta;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEta.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEta.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CutIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIDEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> LHIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEtaHighPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEtaHighPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CutIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CutIdOnlyConvEtaLowPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIDEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> LHIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      LHIdOnlyConvEtaLowPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIDEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyIsoEtaLowPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
      CiCIdOnlyConvEtaLowPt.push_back(aHisto);
    }

//   TH1F *CutIdWP80Eta = new TH1F( "CutIdWP80Eta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *CutIdTightCICEta = new TH1F( "CutIdTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *CutIdSuperTightCICEta = new TH1F( "CutIdSuperTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
//   TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH ID tight #eta",        nbinsEta, minEta, maxEta );
//   TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH ID loose #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 18;
  float minPt = 10.0;
  float maxPt = 100.;
  
  TH1F *FakeableJetsPt   = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtBarrel  = new TH1F( "RecoPtBarrel", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtEndcap  = new TH1F( "RecoPtEndcap", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );

  std::vector<TH1F*> CutIdPt;
  std::vector<TH1F*> CutIdOnlyIDPt;
  std::vector<TH1F*> CutIdOnlyIsoPt;
  std::vector<TH1F*> CutIdOnlyConvPt;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPt;
  std::vector<TH1F*> LHIdOnlyIDPt;
  std::vector<TH1F*> LHIdOnlyIsoPt;
  std::vector<TH1F*> LHIdOnlyConvPt;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPt;
  std::vector<TH1F*> CiCIdOnlyIDPt;
  std::vector<TH1F*> CiCIdOnlyIsoPt;
  std::vector<TH1F*> CiCIdOnlyConvPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPt.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPt.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPtBarrel;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPtBarrel;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIDPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CiCIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPtBarrel.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPtBarrel.push_back(aHisto);
    }

  std::vector<TH1F*> CutIdPtEndcap;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCutBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CutIdOnlyConvPtEndcap.push_back(aHisto);
    }

  std::vector<TH1F*> LHIdPtEndcap;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaLHBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      LHIdOnlyConvPtEndcap.push_back(aHisto);
    }

  std::vector<TH1F*> CiCIdPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIDPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CiCIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCiCBasedID.size();++i)
    {
      TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIDPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyIsoPtEndcap.push_back(aHisto);
      aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
      CiCIdOnlyConvPtEndcap.push_back(aHisto);
    }

//   TH1F *CutIdWP80Pt = new TH1F( "CutIdWP80Pt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH ID tight p_{T} (GeV)",        nbinsPt, minPt, maxPt );
//   TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH ID loose p_{T} (GeV)",        nbinsPt, minPt, maxPt );

  // json 
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;

  // QCD trigger 
  requiredTriggers.push_back("HLT_Jet15U");
  //  requiredTriggers.push_back("HLT_Jet30U");
  //requiredTriggers.push_back("HLT_Jet50U");

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask
    bool newTriggerMask = false;
    if(_isData) newTriggerMask = true;
    reloadTriggerMask(newTriggerMask);
    
    // Good Run selection 
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

    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( !passedHLT ) continue;   


    TVector3 p3PFMET(pxPFMet[0],pyPFMet[0],pzPFMet[0]);
    //MET Cut to reduce W contamination
    if (p3PFMET.Pt()>30.)
      continue;
    // look for the leading jet (not considered as fakeable object to remove trigger bias)
    float maxEt = -1;
    int leadingJet = -1;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 25.0 && p3Jet.Pt() > maxEt) {
        maxEt = p3Jet.Pt();
        leadingJet = jet;
      }
    }

    // consider the other as fakes (denominator)
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 10.0 && jet!=leadingJet ) {
        FakeableJetsEta->Fill( p3Jet.Eta() );
        FakeableJetsPt->Fill( p3Jet.Pt() );
      }
    }

    // consider the electrons near a jet (numerator)
    for ( int ele=0; ele<nEle; ele++ ) {
      TVector3 p3Ele(pxEle[ele], pyEle[ele], pzEle[ele]);
      if ( fabs(etaEle[ele]) < 2.5 && p3Ele.Pt() > 10.0 ) {

        float dREleJet_min = 0.7;
        int closestJet=-1;

        for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
          TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
          if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 10.0 && jet!=leadingJet ) {          
            float dREleJet = p3Jet.DeltaR(p3Ele);
            if(dREleJet<dREleJet_min) {
              closestJet=jet;
              dREleJet_min=dREleJet;
            }
          }
        }

        //        nEleReco++;

        if(closestJet > -1) {

          TVector3 p3ClosestJet(pxAK5PFPUcorrJet[closestJet],pyAK5PFPUcorrJet[closestJet],pzAK5PFPUcorrJet[closestJet]);

//           float etFake=p3ClosestJet.Pt();
//           float etaFake=p3ClosestJet.Eta();

          float etFake=p3Ele.Pt();
          float etaFake=etaEle[ele];

	  bool isInEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEB);
	  bool isInEE = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEE);
	  bool highPt= (etFake>20.);
	  bool lowPt= (etFake<=20.);
	  RecoEta->Fill(etaFake);
	  RecoPt->Fill(etFake);
	  if (highPt) RecoEtaHighPt->Fill(etaFake);
	  if (lowPt) RecoEtaLowPt->Fill(etaFake);
	  if (isInEB) RecoPtBarrel->Fill(etFake);
	  if (isInEE) RecoPtEndcap->Fill(etFake);

      for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCutBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

          // for the smurf selection, add further cut for pT<20 GeV. --> move to WP70
          if(TString(EgammaCutBasedIDWPs[icut].c_str()).Contains("Smurf") && etFake<20.) {
            // apply the VBTF70 smurfs
            isEleID(&EgammaCutBasedID[6],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
            isEleIDCutBased = isEleIDCutBased && (fbremEle[ele]>0.20 || (fabs(etaEle[ele])<1.0 && eSuperClusterOverPEle[ele]>0.95));
          }

	  if ( isEleIDCutBased ) {
	    CutIdOnlyIDEta[icut]->Fill(etaFake);
	    CutIdOnlyIDPt[icut]->Fill(etFake);
	    if (highPt) CutIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CutIdOnlyIDEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CutIdOnlyIDPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CutIdOnlyIDPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isIsolCutBased ) {
	    CutIdOnlyIsoEta[icut]->Fill(etaFake);
	    CutIdOnlyIsoPt[icut]->Fill(etFake);
	    if (highPt) CutIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CutIdOnlyIsoEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CutIdOnlyIsoPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CutIdOnlyIsoPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isConvRejCutBased ) {
	    CutIdOnlyConvEta[icut]->Fill(etaFake);
	    CutIdOnlyConvPt[icut]->Fill(etFake);
	    if (highPt) CutIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CutIdOnlyConvEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CutIdOnlyConvPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CutIdOnlyConvPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CutIdEta[icut]->Fill(etaFake);
	    CutIdPt[icut]->Fill(etFake);
	    if (highPt) CutIdEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CutIdEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CutIdPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CutIdPtEndcap[icut]->Fill(etFake);
	  }
	}
      
//       

//       if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdTightCIC) ) {

//         CutIdTightCICEta->Fill(etaFake);
//         CutIdTightCICPt->Fill(etFake);
        
//       }

//       if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdSuperTightCIC) ) {

//         CutIdSuperTightCICEta->Fill(etaFake);
//         CutIdSuperTightCICPt->Fill(etFake);
        
//       }
      
      for (int icut=0;icut<EgammaLHBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaLHBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    LHIdOnlyIDEta[icut]->Fill(etaFake);
	    LHIdOnlyIDPt[icut]->Fill(etFake);
	    if (highPt) LHIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) LHIdOnlyIDEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) LHIdOnlyIDPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdOnlyIDPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isIsolCutBased ) {
	    LHIdOnlyIsoEta[icut]->Fill(etaFake);
	    LHIdOnlyIsoPt[icut]->Fill(etFake);
	    if (highPt) LHIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) LHIdOnlyIsoEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) LHIdOnlyIsoPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdOnlyIsoPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isConvRejCutBased ) {
	    LHIdOnlyConvEta[icut]->Fill(etaFake);
	    LHIdOnlyConvPt[icut]->Fill(etFake);
	    if (highPt) LHIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) LHIdOnlyConvEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) LHIdOnlyConvPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdOnlyConvPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    LHIdEta[icut]->Fill(etaFake);
	    LHIdPt[icut]->Fill(etFake);
	    if (highPt) LHIdEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) LHIdEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) LHIdPtBarrel[icut]->Fill(etFake);
	    if (isInEE) LHIdPtEndcap[icut]->Fill(etFake);
	  }	  

	}

      for (int icut=0;icut<EgammaCiCBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCiCBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    CiCIdOnlyIDEta[icut]->Fill(etaFake);
	    CiCIdOnlyIDPt[icut]->Fill(etFake);
	    if (highPt) CiCIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CiCIdOnlyIDEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CiCIdOnlyIDPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdOnlyIDPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isIsolCutBased ) {
	    CiCIdOnlyIsoEta[icut]->Fill(etaFake);
	    CiCIdOnlyIsoPt[icut]->Fill(etFake);
	    if (highPt) CiCIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CiCIdOnlyIsoEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CiCIdOnlyIsoPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdOnlyIsoPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isConvRejCutBased ) {
	    CiCIdOnlyConvEta[icut]->Fill(etaFake);
	    CiCIdOnlyConvPt[icut]->Fill(etFake);
	    if (highPt) CiCIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CiCIdOnlyConvEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CiCIdOnlyConvPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdOnlyConvPtEndcap[icut]->Fill(etFake);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CiCIdEta[icut]->Fill(etaFake);
	    CiCIdPt[icut]->Fill(etFake);
	    if (highPt) CiCIdEtaHighPt[icut]->Fill(etaFake);
	    if (lowPt) CiCIdEtaLowPt[icut]->Fill(etaFake);
	    if (isInEB) CiCIdPtBarrel[icut]->Fill(etFake);
	    if (isInEE) CiCIdPtEndcap[icut]->Fill(etFake);
	  }	  	  

	}
	
          
//           Utils anaUtils;
          
//           if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdTightCIC) ) {
            
//             CutIdTightCICEta->Fill(etaFake);
//             CutIdTightCICPt->Fill(etFake);
            
//           }

//           if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdSuperTightCIC) ) {
            
//             CutIdSuperTightCICEta->Fill(etaFake);
//             CutIdSuperTightCICPt->Fill(etFake);
            
//           }

//           float lik = likelihoodRatio(ele,*LH);
//           bool eleInGap = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEBEEGap);

//           if ( lik > minLikelihoodLoose && !eleInGap) {

//             LHIdLooseEta->Fill(etaFake);
//             LHIdLoosePt->Fill(etFake);

//           }

//           if ( lik > minLikelihoodTight && !eleInGap) {

//             LHIdTightEta->Fill(etaFake);
//             LHIdTightPt->Fill(etFake);

//           }

        }

      } // electron acceptance & pt cut
      
    } // loop ele

    //    cout << "nEleReco = " << nEleReco << " nJetsReco = " << nJetsReco << endl;
    
  } // loop events
  

  char filename[200];
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  //  ElectronEffEta.AddNumerator(GenEta);
  ElectronEffEta.AddNumerator(RecoEta);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(CutIdEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(CutIdOnlyConvEta[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(LHIdEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(LHIdOnlyConvEta[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEta.AddNumerator(CiCIdEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyIDEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyIsoEta[icut]);
      ElectronEffEta.AddNumerator(CiCIdOnlyConvEta[icut]);
    }
  ElectronEffEta.SetDenominator(RecoEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("fake rate vs #eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("Fake rate");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleMisidEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  //  ElectronEffEtaHighPt.AddNumerator(GenEtaHighPt);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(CutIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CutIdOnlyConvEtaHighPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(LHIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(LHIdOnlyConvEtaHighPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEtaHighPt.AddNumerator(CiCIdEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIDEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIsoEtaHighPt[icut]);
      ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyConvEtaHighPt[icut]);
    }
  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaHighPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  //  ElectronEffEtaLowPt.AddNumerator(GenEtaLowPt);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(CutIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CutIdOnlyConvEtaLowPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(LHIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(LHIdOnlyConvEtaLowPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffEtaLowPt.AddNumerator(CiCIdEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIDEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIsoEtaLowPt[icut]);
      ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyConvEtaLowPt[icut]);
    }
  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  //  ElectronEffPt.AddNumerator(GenPt);
  ElectronEffPt.AddNumerator(RecoPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(CutIdPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(CutIdOnlyConvPt[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(LHIdPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(LHIdOnlyConvPt[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPt.AddNumerator(CiCIdPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyIDPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyIsoPt[icut]);
      ElectronEffPt.AddNumerator(CiCIdOnlyConvPt[icut]);
    }

  ElectronEffPt.SetDenominator(RecoPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("fake rate vs p_{T}");
  ElectronEffPt.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPt.SetYaxisTitle("Fake rate");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleMisidPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  //  ElectronEffPtBarrel.AddNumerator(GenPtBarrel);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(CutIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CutIdOnlyConvPtBarrel[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(LHIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(LHIdOnlyConvPtBarrel[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPtBarrel.AddNumerator(CiCIdPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIDPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIsoPtBarrel[icut]);
      ElectronEffPtBarrel.AddNumerator(CiCIdOnlyConvPtBarrel[icut]);
    }

  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("fake rate vs p_{T}");
  ElectronEffPtBarrel.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-EleMisidPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  //  ElectronEffPtEndcap.AddNumerator(GenPtEndcap);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(CutIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CutIdOnlyConvPtEndcap[icut]);
    }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(LHIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(LHIdOnlyConvPtEndcap[icut]);
    }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut)
    {
      ElectronEffPtEndcap.AddNumerator(CiCIdPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIDPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIsoPtEndcap[icut]);
      ElectronEffPtEndcap.AddNumerator(CiCIdOnlyConvPtEndcap[icut]);
    }

  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("fake rate vs p_{T}");
  ElectronEffPtEndcap.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();
}

void LikelihoodAnalysis::estimateFakeRateForHToWW_QCD(const char *outname) {

  // study vs eta
  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
  TH1F *FakeableJetsEta = new TH1F( "FakeableJetsEta", "fakeable jets #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEta         = new TH1F( "RecoEta",         "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaHighPt   = new TH1F( "RecoEtaHighPt",   "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *RecoEtaLowPt    = new TH1F( "RecoEtaLowPt",    "reconstructed #eta", nbinsEta, minEta, maxEta );

  std::vector<TH1F*> CutIdEta;
  std::vector<TH1F*> CutIdOnlyIDEta;
  std::vector<TH1F*> CutIdOnlyIsoEta;
  std::vector<TH1F*> CutIdOnlyConvEta;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta", "cut ID #eta",   nbinsEta, minEta, maxEta );
    CutIdEta.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta", "cut ID #eta",   nbinsEta, minEta, maxEta );
    CutIdOnlyIDEta.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEta", "cut ID #eta",  nbinsEta, minEta, maxEta );
    CutIdOnlyIsoEta.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
    CutIdOnlyConvEta.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdEta;
  std::vector<TH1F*> LHIdOnlyIDEta;
  std::vector<TH1F*> LHIdOnlyIsoEta;
  std::vector<TH1F*> LHIdOnlyConvEta;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta",   "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdEta.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta",   "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyIDEta.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEta",  "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyIsoEta.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyConvEta.push_back(aHisto);
  }

  std::vector<TH1F*> CiCIdEta;
  std::vector<TH1F*> CiCIdOnlyIDEta;
  std::vector<TH1F*> CiCIdOnlyIsoEta;
  std::vector<TH1F*> CiCIdOnlyConvEta;
  for (int i=0;i<EgammaCiCBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta",   "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdEta.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEta",   "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyIDEta.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEta",  "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyIsoEta.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyConvEta.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CutIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaHighPt", "cut ID #eta",   nbinsEta, minEta, maxEta );
    CutIdEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta",   nbinsEta, minEta, maxEta );
    CutIdOnlyIDEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta",  nbinsEta, minEta, maxEta );
    CutIdOnlyIsoEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
    CutIdOnlyConvEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIDEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> LHIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaHighPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaHighPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyIDEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaHighPt",  "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyIsoEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyConvEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> CiCIdEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaHighPt", "cut ID #eta",   nbinsEta, minEta, maxEta );
    CiCIdEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaHighPt", "cut ID #eta",   nbinsEta, minEta, maxEta );
    CiCIdOnlyIDEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaHighPt", "cut ID #eta",  nbinsEta, minEta, maxEta );
    CiCIdOnlyIsoEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyConvEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CutIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaLowPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    CutIdEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaLowPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    CutIdOnlyIDEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaLowPt",  "cut ID #eta", nbinsEta, minEta, maxEta );
    CutIdOnlyIsoEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
    CutIdOnlyConvEtaLowPt.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIDEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> LHIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaLowPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaLowPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyIDEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaLowPt",  "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyIsoEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
    LHIdOnlyConvEtaLowPt.push_back(aHisto);
  }

  std::vector<TH1F*> CiCIdEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CiCIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"EtaLowPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEtaLowPt",   "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyIDEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoEtaLowPt",  "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyIsoEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", nbinsEta, minEta, maxEta );
    CiCIdOnlyConvEtaLowPt.push_back(aHisto);
  }


  // study vs pT
  // int nbinsPt = 9;
  // float minPt = 10.;
  // float maxPt = 80.;
  int nbinsPt = 18;
  float minPt = 10.0;
  float maxPt = 100.;
  TH1F *FakeableJetsPt = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPt         = new TH1F( "RecoPt",          "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );

  std::vector<TH1F*> CutIdPt;
  std::vector<TH1F*> CutIdOnlyIDPt;
  std::vector<TH1F*> CutIdOnlyIsoPt;
  std::vector<TH1F*> CutIdOnlyConvPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt",   "cut ID #eta", nbinsPt, minPt, maxPt );
    CutIdPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPt",   "cut ID #eta", nbinsPt, minPt, maxPt );
    CutIdOnlyIDPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPt",  "cut ID #eta", nbinsPt, minPt, maxPt );
    CutIdOnlyIsoPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
    CutIdOnlyConvPt.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdPt;
  std::vector<TH1F*> LHIdOnlyIDPt;
  std::vector<TH1F*> LHIdOnlyIsoPt;
  std::vector<TH1F*> LHIdOnlyConvPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt", "cut ID #eta",   nbinsPt, minPt, maxPt );
    LHIdPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta",   nbinsPt, minPt, maxPt );
    LHIdOnlyIDPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta",  nbinsPt, minPt, maxPt );
    LHIdOnlyIsoPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
    LHIdOnlyConvPt.push_back(aHisto);
  }

  std::vector<TH1F*> CiCIdPt;
  std::vector<TH1F*> CiCIdOnlyIDPt;
  std::vector<TH1F*> CiCIdOnlyIsoPt;
  std::vector<TH1F*> CiCIdOnlyConvPt;
  for (int i=0;i<EgammaCiCBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CiCIdPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPt", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CiCIdOnlyIDPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPt", "cut ID #eta",  nbinsPt, minPt, maxPt );
    CiCIdOnlyIsoPt.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", nbinsPt, minPt, maxPt );
    CiCIdOnlyConvPt.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdPtBarrel;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CutIdPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CutIdOnlyIDPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta",  nbinsPt, minPt, maxPt );
    CutIdOnlyIsoPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
    CutIdOnlyConvPtBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdPtBarrel;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel", "cut ID #eta",   nbinsPt, minPt, maxPt );
    LHIdPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta",   nbinsPt, minPt, maxPt );
    LHIdOnlyIDPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta",  nbinsPt, minPt, maxPt );
    LHIdOnlyIsoPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
    LHIdOnlyConvPtBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> CiCIdPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIDPtBarrel;
  std::vector<TH1F*> CiCIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CiCIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCiCBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtBarrel", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CiCIdPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CiCIdOnlyIDPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta",  nbinsPt, minPt, maxPt );
    CiCIdOnlyIsoPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", nbinsPt, minPt, maxPt );
    CiCIdOnlyConvPtBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdPtEndcap;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CutIdPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CutIdOnlyIDPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta",  nbinsPt, minPt, maxPt );
    CutIdOnlyIsoPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
    CutIdOnlyConvPtEndcap.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdPtEndcap;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap", "cut ID #eta",   nbinsPt, minPt, maxPt );
    LHIdPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta",   nbinsPt, minPt, maxPt );
    LHIdOnlyIDPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta",  nbinsPt, minPt, maxPt );
    LHIdOnlyIsoPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
    LHIdOnlyConvPtEndcap.push_back(aHisto);
  }

  std::vector<TH1F*> CiCIdPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIDPtEndcap;
  std::vector<TH1F*> CiCIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CiCIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCiCBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PtEndcap", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CiCIdPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta",   nbinsPt, minPt, maxPt );
    CiCIdOnlyIDPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta",  nbinsPt, minPt, maxPt );
    CiCIdOnlyIsoPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", nbinsPt, minPt, maxPt );
    CiCIdOnlyConvPtEndcap.push_back(aHisto);
  }

  
  // trigger: QCD
  cout << "using di-jet triggers" << endl;
  // 2010 data
  // requiredTriggers.push_back("HLT_Jet50U");
  // requiredTriggers.push_back("HLT_Jet50U_v1");
  // requiredTriggers.push_back("HLT_Jet50U_v2");
  // requiredTriggers.push_back("HLT_Jet50U_v3");
  // 2011 data
  requiredTriggers.push_back("HLT_Jet80_v1");

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
    
    // all events
    myCounter.IncrVar("event",1);

    // event selection: trigger
    reloadTriggerMask(true);
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( !passedHLT ) continue;   
    myCounter.IncrVar("trigger",1);    

    // electrons passing the denominator selection to reduce W and Z contamination
    std::vector<int> denomElectrons;
    for(int iele=0; iele<nEle; iele++) {
      
      bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[iele], bits::isEcalDriven);
      if(!ecalDriven) continue;   
      
      // combined isol 
      TVector3 p3Ele(pxEle[iele], pyEle[iele], pzEle[iele]);
      float pt = p3Ele.Pt();
      float myCombIso = (fabs(etaEle[iele])<1.479) ? (dr03TkSumPtEle[iele] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[iele]-1.0) + dr03HcalTowerSumEtFullConeEle[iele])
      	: ( dr03TkSumPtEle[iele] + dr03EcalRecHitSumEtEle[iele] + dr03HcalTowerSumEtFullConeEle[iele]);
      float combIso = (myCombIso - rhoFastjet*TMath::Pi()*0.3*0.3) / pt;

      if(combIso>0.15) continue;    
      
      // H/E
      if(hOverEEle[iele]>0.15) continue;
      
      // sigmaIetaIeta
      bool isBarrelSc;
      int sc = superClusterIndexEle[iele];
      if ( sc < 0 ) continue;
      if ( fabs(etaSC[sc]) <  1.479 ) isBarrelSc = true;
      if ( fabs(etaSC[sc]) >= 1.479 ) isBarrelSc = false;
      if ( isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.014 ) continue;   
      if (!isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.034 ) continue;   
      
      // spikes 
      float theE1 = eMaxSC[sc];
      float theE4SwissCross = e4SwissCrossSC[sc];
      float theSpikeSC = 1.0 - (theE4SwissCross/theE1);
      if (theSpikeSC>0.95) continue;
      
      denomElectrons.push_back(iele);
    } 

    // event selection: at least one candidate for denominator
    if (denomElectrons.size()==0) continue;
    myCounter.IncrVar("denom",1);    

    // event selection: met cut to reduce W->enu
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    if( p3Met.Pt() > 20 ) continue;       
    myCounter.IncrVar("met",1);    

    // event selection: mT cut to reduce W->enu [highest pT denominator candidate used]
    std::pair<int,int> possibleWcand = getBestGoodElePair(denomElectrons);
    int theWEle1(possibleWcand.first);
    TLorentzVector tlvWEle1;
    TVector3 tv3Ele1;
    tlvWEle1.SetXYZT(pxEle[theWEle1],pyEle[theWEle1],pzEle[theWEle1],energyEle[theWEle1]);
    tv3Ele1.SetXYZ(pxEle[theWEle1],pyEle[theWEle1],pzEle[theWEle1]);
    float WmT = sqrt(2*tlvWEle1.Pt()*p3Met.Pt()*(1-cos(tv3Ele1.Angle(p3Met))) );      
    if (WmT > 25. ) continue;
    myCounter.IncrVar("trasvMass",1);    

    // event selection: invariant mass between two electrons passing the denominator selection to reduce Z->ee
    std::pair<int,int> possibleZcand = getBestGoodElePair(denomElectrons);
    int theZEle1(possibleZcand.first);
    int theZEle2(possibleZcand.second);
    TLorentzVector tlvZEle1, tlvZEle2;
    tlvZEle1.SetXYZT(pxEle[theZEle1],pyEle[theZEle1],pzEle[theZEle1],energyEle[theZEle1]);
    tlvZEle2.SetXYZT(pxEle[theZEle2],pyEle[theZEle2],pzEle[theZEle2],energyEle[theZEle2]);
    double theInvMass = (tlvZEle1+tlvZEle2).M();
    if (theInvMass>60. && theInvMass<120.) continue;
    myCounter.IncrVar("Zmass",1);    

    // look for the leading jet
    // (not considered as fakeable object to remove trigger bias)
    float maxEt = -1.;
    int leading = -1;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.Eta()) < maxEta && p3Jet.Pt() > maxEt) {  
	maxEt   = p3Jet.Pt();
	leading = jet;
      }
    }
    
    // need a triggering object 
    if ( leading < 0 ) continue;
    myCounter.IncrVar("leadingExist",1);    

    // minimal cut on leading jet ET
    TVector3 p3LeadingCut(pxAK5PFPUcorrJet[leading],pyAK5PFPUcorrJet[leading],pzAK5PFPUcorrJet[leading]);
    if (p3LeadingCut.Pt()<55) continue;

    myCounter.IncrVar("leadingPT",1);    

    
    // consider the reco electrons not matching the leading jet as fakes (denominator)
    // with the following requirements:
    // Gsf Electron (Ecal driven)
    // - relative Iso < 0.1
    // - H/E < 0.15
    // - sigma ieta ieta < 0.014 (0.034) 
    // cleaning for spikes

    // denominator and numerator
    for(int iele=0; iele<nEle; iele++) {

      TVector3 p3Ele(pxEle[iele], pyEle[iele], pzEle[iele]);

      // not matching the leading object
      float dr;
      TVector3 p3LeadingJet(pxAK5PFPUcorrJet[leading],pyAK5PFPUcorrJet[leading],pzAK5PFPUcorrJet[leading]);
      dr=p3LeadingJet.DeltaR(p3Ele);
      
      // denominator selection
      bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[iele], bits::isEcalDriven);
      if(!ecalDriven) continue;
      
      // combined isol 
      float pt = p3Ele.Pt();
      float myCombIso = (fabs(etaEle[iele])<1.479) ? (dr03TkSumPtEle[iele] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[iele]-1.0) + dr03HcalTowerSumEtFullConeEle[iele])
      	: ( dr03TkSumPtEle[iele] + dr03EcalRecHitSumEtEle[iele] + dr03HcalTowerSumEtFullConeEle[iele]);
      float combIso = (myCombIso - rhoFastjet*TMath::Pi()*0.3*0.3) / pt;

      if(combIso>0.1) continue;   
      
      // H/E
      if(hOverEEle[iele]>0.15) continue;  
      
      // sigmaIetaIeta
      bool isBarrelSc;
      int sc = superClusterIndexEle[iele];
      if ( sc < 0 ) continue;
      if ( fabs(etaSC[sc]) <  1.479 ) isBarrelSc = true;
      if ( fabs(etaSC[sc]) >= 1.479 ) isBarrelSc = false;
      if ( isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.014 ) continue;   
      if (!isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.034 ) continue;   
      
      // spikes 
      float theE1 = eMaxSC[sc];
      float theE4SwissCross = e4SwissCrossSC[sc];
      float theSpikeSC = 1.0 - (theE4SwissCross/theE1);
      if (theSpikeSC>0.95) continue;

      // end denominator selection

      // exclude the object corresponding to (probably) triggering jet: avoid trigger bias
      if( dr < 0.3 ) continue;
      
      // only electrons within the acceptance
      if( fabs(p3Ele.Eta()) > maxEta ) continue;
      if( p3Ele.Pt() < minPt ) continue;

      // fill the denominator
      float etaFake = p3Ele.Eta();
      float etFake  = p3Ele.Pt();
      bool isInEB   = anaUtils.fiducialFlagECAL(fiducialFlagsEle[iele], isEB);
      bool isInEE   = anaUtils.fiducialFlagECAL(fiducialFlagsEle[iele], isEE);
      bool highPt   = (etFake>20.);
      bool lowPt    = (etFake<=20.);
      RecoEta         -> Fill(etaFake);
      RecoPt          -> Fill(etFake);
      FakeableJetsEta -> Fill( etaFake );
      FakeableJetsPt  -> Fill( etFake );
      if (highPt) RecoEtaHighPt->Fill(etaFake);
      if (lowPt)  RecoEtaLowPt ->Fill(etaFake);
      if (isInEB) RecoPtBarrel ->Fill(etFake);
      if (isInEE) RecoPtEndcap ->Fill(etFake);
      
      
      // numerator: simple cut based
      for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut) {
	
        bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
        isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
        isEleID(&EgammaCutBasedID[icut],iele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
	
        // for the smurf selection, add further cut for pT<20 GeV. --> move to WP70
        if(TString(EgammaCutBasedIDWPs[icut].c_str()).Contains("Smurf") && etFake<20.) {
          // apply the VBTF70 smurfs
          isEleID(&EgammaCutBasedID[6],iele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
          isEleIDCutBased = isEleIDCutBased && (fbremEle[iele]>0.15 || (fabs(etaEle[iele])<1.0 && eSuperClusterOverPEle[iele]>0.95));
        }

        if ( isEleIDCutBased ) {
          CutIdOnlyIDEta[icut]->Fill(etaFake);
          CutIdOnlyIDPt[icut] ->Fill(etFake);
          if (highPt) CutIdOnlyIDEtaHighPt[icut]->Fill(etaFake);
          if (lowPt)  CutIdOnlyIDEtaLowPt[icut] ->Fill(etaFake);
          if (isInEB) CutIdOnlyIDPtBarrel[icut] ->Fill(etFake);
          if (isInEE) CutIdOnlyIDPtEndcap[icut] ->Fill(etFake);
	}
	
	if ( isIsolCutBased ) {
          CutIdOnlyIsoEta[icut]->Fill(etaFake);
          CutIdOnlyIsoPt[icut]->Fill(etFake);
          if (highPt) CutIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
          if (lowPt)  CutIdOnlyIsoEtaLowPt[icut] ->Fill(etaFake);
          if (isInEB) CutIdOnlyIsoPtBarrel[icut] ->Fill(etFake);
          if (isInEE) CutIdOnlyIsoPtEndcap[icut] ->Fill(etFake);
        }

        if ( isConvRejCutBased ) {
          CutIdOnlyConvEta[icut]->Fill(etaFake);
          CutIdOnlyConvPt[icut]->Fill(etFake);
          if (highPt) CutIdOnlyConvEtaHighPt[icut]->Fill(etaFake);
          if (lowPt)  CutIdOnlyConvEtaLowPt[icut] ->Fill(etaFake);
          if (isInEB) CutIdOnlyConvPtBarrel[icut] ->Fill(etFake);
          if (isInEE) CutIdOnlyConvPtEndcap[icut] ->Fill(etFake);
        }

        if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
          CutIdEta[icut]->Fill(etaFake);
          CutIdPt[icut] ->Fill(etFake);
          if (highPt) CutIdEtaHighPt[icut]->Fill(etaFake);
          if (lowPt)  CutIdEtaLowPt[icut] ->Fill(etaFake);
          if (isInEB) CutIdPtBarrel[icut] ->Fill(etFake);
          if (isInEE) CutIdPtEndcap[icut] ->Fill(etFake);
        }
      }

      // numerator: likelihood based
      for (int icut=0;icut<EgammaLHBasedIDWPs.size();++icut) {

        bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
        isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
        isEleID(&EgammaLHBasedID[icut],iele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

        if ( isEleIDCutBased ) {
          LHIdOnlyIDEta[icut]->Fill(etaFake);
          LHIdOnlyIDPt[icut] ->Fill(etFake);
          if (highPt) LHIdOnlyIDEtaHighPt[icut]-> Fill(etaFake);
          if (lowPt)  LHIdOnlyIDEtaLowPt[icut] -> Fill(etaFake);
          if (isInEB) LHIdOnlyIDPtBarrel[icut] -> Fill(etFake);
          if (isInEE) LHIdOnlyIDPtEndcap[icut] -> Fill(etFake);
        }

        if ( isIsolCutBased ) {
          LHIdOnlyIsoEta[icut]->Fill(etaFake);
          LHIdOnlyIsoPt[icut]->Fill(etFake);
          if (highPt) LHIdOnlyIsoEtaHighPt[icut]-> Fill(etaFake);
          if (lowPt)  LHIdOnlyIsoEtaLowPt[icut] -> Fill(etaFake);
          if (isInEB) LHIdOnlyIsoPtBarrel[icut] -> Fill(etFake);
          if (isInEE) LHIdOnlyIsoPtEndcap[icut] -> Fill(etFake);
        }

	if ( isConvRejCutBased ) {
          LHIdOnlyConvEta[icut]->Fill(etaFake);
          LHIdOnlyConvPt[icut] ->Fill(etFake);
          if (highPt) LHIdOnlyConvEtaHighPt[icut]-> Fill(etaFake);
          if (lowPt)  LHIdOnlyConvEtaLowPt[icut] -> Fill(etaFake);
          if (isInEB) LHIdOnlyConvPtBarrel[icut] -> Fill(etFake);
          if (isInEE) LHIdOnlyConvPtEndcap[icut] -> Fill(etFake);
        }

        if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
          LHIdEta[icut]->Fill(etaFake);
          LHIdPt[icut] ->Fill(etFake);
          if (highPt) LHIdEtaHighPt[icut]->Fill(etaFake);
          if (lowPt)  LHIdEtaLowPt[icut] ->Fill(etaFake);
          if (isInEB) LHIdPtBarrel[icut] ->Fill(etFake);
          if (isInEE) LHIdPtEndcap[icut] ->Fill(etFake);
        }
      }

     // numerator: CIC based
      for (int icut=0;icut<EgammaCiCBasedIDWPs.size();++icut) {

        bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
        isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
        isEleID(&EgammaCiCBasedID[icut],iele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

        if ( isEleIDCutBased ) {
          CiCIdOnlyIDEta[icut]->Fill(etaFake);
          CiCIdOnlyIDPt[icut] ->Fill(etFake);
          if (highPt) CiCIdOnlyIDEtaHighPt[icut]-> Fill(etaFake);
          if (lowPt)  CiCIdOnlyIDEtaLowPt[icut] -> Fill(etaFake);
          if (isInEB) CiCIdOnlyIDPtBarrel[icut] -> Fill(etFake);
          if (isInEE) CiCIdOnlyIDPtEndcap[icut] -> Fill(etFake);
        }

        if ( isIsolCutBased ) {
          CiCIdOnlyIsoEta[icut]->Fill(etaFake);
          CiCIdOnlyIsoPt[icut] ->Fill(etFake);
          if (highPt) CiCIdOnlyIsoEtaHighPt[icut]->Fill(etaFake);
          if (lowPt)  CiCIdOnlyIsoEtaLowPt[icut] ->Fill(etaFake);
          if (isInEB) CiCIdOnlyIsoPtBarrel[icut] ->Fill(etFake);
          if (isInEE) CiCIdOnlyIsoPtEndcap[icut] ->Fill(etFake);
        }
 
	if ( isConvRejCutBased ) {
          CiCIdOnlyConvEta[icut]->Fill(etaFake);
          CiCIdOnlyConvPt[icut] ->Fill(etFake);
          if (highPt) CiCIdOnlyConvEtaHighPt[icut] -> Fill(etaFake);
          if (lowPt)  CiCIdOnlyConvEtaLowPt[icut]  -> Fill(etaFake);
          if (isInEB) CiCIdOnlyConvPtBarrel[icut]  -> Fill(etFake);
          if (isInEE) CiCIdOnlyConvPtEndcap[icut]  -> Fill(etFake);
        }

        if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
          CiCIdEta[icut]->Fill(etaFake);
          CiCIdPt[icut] ->Fill(etFake);
          if (highPt) CiCIdEtaHighPt[icut] -> Fill(etaFake);
          if (lowPt)  CiCIdEtaLowPt[icut]  -> Fill(etaFake);
          if (isInEB) CiCIdPtBarrel[icut]  -> Fill(etFake);
          if (isInEE) CiCIdPtEndcap[icut]  -> Fill(etFake);
        }
      
      }

    } // loop ele

  } // loop events

  
  // saving the counters
  char filename[200];
  sprintf(filename,"%sCounters.root",outname);
  myCounter.Save(filename,"recreate");

  // saving efficiency histos
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  ElectronEffEta.AddNumerator(FakeableJetsEta);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffEta.AddNumerator(CutIdEta[icut]);
    ElectronEffEta.AddNumerator(CutIdOnlyIDEta[icut]);
    ElectronEffEta.AddNumerator(CutIdOnlyIsoEta[icut]);
    ElectronEffEta.AddNumerator(CutIdOnlyConvEta[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut){
    ElectronEffEta.AddNumerator(LHIdEta[icut]);
    ElectronEffEta.AddNumerator(LHIdOnlyIDEta[icut]);
    ElectronEffEta.AddNumerator(LHIdOnlyIsoEta[icut]);
    ElectronEffEta.AddNumerator(LHIdOnlyConvEta[icut]);
  }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut) {
    ElectronEffEta.AddNumerator(CiCIdEta[icut]);
    ElectronEffEta.AddNumerator(CiCIdOnlyIDEta[icut]);
    ElectronEffEta.AddNumerator(CiCIdOnlyIsoEta[icut]);
    ElectronEffEta.AddNumerator(CiCIdOnlyConvEta[icut]);
  }

  ElectronEffEta.SetDenominator(FakeableJetsEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("fake rate vs #eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("Fake rate");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleMisidEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffEtaHighPt.AddNumerator(CutIdEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIDEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIsoEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CutIdOnlyConvEtaHighPt[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffEtaHighPt.AddNumerator(LHIdEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIDEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIsoEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(LHIdOnlyConvEtaHighPt[icut]);
  }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut) {
    ElectronEffEtaHighPt.AddNumerator(CiCIdEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIDEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyIsoEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CiCIdOnlyConvEtaHighPt[icut]);
  }

  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaHighPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffEtaLowPt.AddNumerator(CutIdEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIDEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIsoEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CutIdOnlyConvEtaLowPt[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffEtaLowPt.AddNumerator(LHIdEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIDEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIsoEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(LHIdOnlyConvEtaLowPt[icut]);
  }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut) {
    ElectronEffEtaLowPt.AddNumerator(CiCIdEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIDEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyIsoEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CiCIdOnlyConvEtaLowPt[icut]);
  }
  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  ElectronEffPt.AddNumerator(FakeableJetsPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPt.AddNumerator(CutIdPt[icut]);
    ElectronEffPt.AddNumerator(CutIdOnlyIDPt[icut]);
    ElectronEffPt.AddNumerator(CutIdOnlyIsoPt[icut]);
    ElectronEffPt.AddNumerator(CutIdOnlyConvPt[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPt.AddNumerator(LHIdPt[icut]);
    ElectronEffPt.AddNumerator(LHIdOnlyIDPt[icut]);
    ElectronEffPt.AddNumerator(LHIdOnlyIsoPt[icut]);
    ElectronEffPt.AddNumerator(LHIdOnlyConvPt[icut]);
  }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut) {
    ElectronEffPt.AddNumerator(CiCIdPt[icut]);
    ElectronEffPt.AddNumerator(CiCIdOnlyIDPt[icut]);
    ElectronEffPt.AddNumerator(CiCIdOnlyIsoPt[icut]);
    ElectronEffPt.AddNumerator(CiCIdOnlyConvPt[icut]);
  }
  ElectronEffPt.SetDenominator(FakeableJetsPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("fake rate vs p_{T}");
  ElectronEffPt.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPt.SetYaxisTitle("Fake rate");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleMisidPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtBarrel.AddNumerator(CutIdPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CutIdOnlyIDPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CutIdOnlyIsoPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CutIdOnlyConvPtBarrel[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtBarrel.AddNumerator(LHIdPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(LHIdOnlyIDPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(LHIdOnlyIsoPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(LHIdOnlyConvPtBarrel[icut]);
  }
  for (int icut=0;icut<EgammaCiCBasedID.size();++icut) {
    ElectronEffPtBarrel.AddNumerator(CiCIdPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIDPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CiCIdOnlyIsoPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CiCIdOnlyConvPtBarrel[icut]);
  }
  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("fake rate vs p_{T}");
  ElectronEffPtBarrel.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-EleMisidPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtEndcap.AddNumerator(CutIdPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CutIdOnlyIDPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CutIdOnlyIsoPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CutIdOnlyConvPtEndcap[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtEndcap.AddNumerator(LHIdPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(LHIdOnlyIDPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(LHIdOnlyIsoPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(LHIdOnlyConvPtEndcap[icut]);
  }
 for (int icut=0;icut<EgammaCiCBasedID.size();++icut) {
    ElectronEffPtEndcap.AddNumerator(CiCIdPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIDPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CiCIdOnlyIsoPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CiCIdOnlyConvPtEndcap[icut]);
  }
  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("fake rate vs p_{T}");
  ElectronEffPtEndcap.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();
}

void LikelihoodAnalysis::estimateFakeRateForHToWW_EGAMMA(const char *outname) {

  // to study which ET cut we want to apply on the jet (for Wjets HWW studies)
  TH1F *FakeJetForThresholdPT = new TH1F("FakeJetForThresholdPT", "fakeable jets p_{T}", 20, 0., 100.);


  // -----------------------------------------------------------------------
  // study vs eta
  float minEta = -2.5;
  float maxEta =  2.5;

  Float_t LowerEta[5];
  LowerEta[0]=0.0;
  LowerEta[1]=1.0;
  LowerEta[2]=1.479;  
  LowerEta[3]=2.0;
  LowerEta[4]=2.5;
  TH1F *FakeableJetsEta = new TH1F( "FakeableJetsEta",  "fakeable jets #eta", 4, LowerEta);
  TH1F *RecoEta         = new TH1F( "RecoEta",          "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaHighPt   = new TH1F( "RecoEtaHighPt",    "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaLowPt    = new TH1F( "RecoEtaLowPt",     "reconstructed #eta", 4, LowerEta);
  
  std::vector<TH1F*> CutIdEta;
  std::vector<TH1F*> CutIdOnlyIDEta;
  std::vector<TH1F*> CutIdOnlyIsoEta;
  std::vector<TH1F*> CutIdOnlyConvEta;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta",   "cut ID #eta", 4, LowerEta);
    CutIdEta.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta",   "cut ID #eta", 4, LowerEta);   
    CutIdOnlyIDEta.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEta",  "cut ID #eta", 4, LowerEta);  
    CutIdOnlyIsoEta.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEta", "cut ID #eta", 4, LowerEta); 
    CutIdOnlyConvEta.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdEta;
  std::vector<TH1F*> LHIdOnlyIDEta;
  std::vector<TH1F*> LHIdOnlyIsoEta;
  std::vector<TH1F*> LHIdOnlyConvEta;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta",   "like ID #eta", 4, LowerEta); 
    LHIdEta.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta",   "like ID #eta", 4, LowerEta); 
    LHIdOnlyIDEta.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEta",  "like ID #eta", 4, LowerEta); 
    LHIdOnlyIsoEta.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEta", "like ID #eta", 4, LowerEta); 
    LHIdOnlyConvEta.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIDEtaHighPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> CutIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaHighPt",   "cut ID #eta", 4, LowerEta);   
    CutIdEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaHighPt",   "cut ID #eta", 4, LowerEta);   
    CutIdOnlyIDEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaHighPt",  "cut ID #eta", 4, LowerEta);  
    CutIdOnlyIsoEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaHighPt", "cut ID #eta", 4, LowerEta); 
    CutIdOnlyConvEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIDEtaHighPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaHighPt;
  std::vector<TH1F*> LHIdOnlyConvEtaHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaHighPt",   "like ID #eta", 4, LowerEta); 
    LHIdEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaHighPt",   "like ID #eta", 4, LowerEta); 
    LHIdOnlyIDEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaHighPt",  "like ID #eta", 4, LowerEta); 
    LHIdOnlyIsoEtaHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaHighPt", "like ID #eta", 4, LowerEta); 
    LHIdOnlyConvEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIDEtaLowPt;
  std::vector<TH1F*> CutIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> CutIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"EtaLowPt",   "cut ID #eta", 4, LowerEta);  
    CutIdEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEtaLowPt",   "cut ID #eta", 4, LowerEta); 
    CutIdOnlyIDEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoEtaLowPt",  "cut ID #eta", 4, LowerEta); 
    CutIdOnlyIsoEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvEtaLowPt", "cut ID #eta", 4, LowerEta); 
    CutIdOnlyConvEtaLowPt.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIDEtaLowPt;
  std::vector<TH1F*> LHIdOnlyIsoEtaLowPt;
  std::vector<TH1F*> LHIdOnlyConvEtaLowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"EtaLowPt",   "like ID #eta", 4, LowerEta);
    LHIdEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEtaLowPt",   "like ID #eta", 4, LowerEta);
    LHIdOnlyIDEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoEtaLowPt",  "like ID #eta", 4, LowerEta);
    LHIdOnlyIsoEtaLowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvEtaLowPt", "like ID #eta", 4, LowerEta);
    LHIdOnlyConvEtaLowPt.push_back(aHisto);
  }


  // -----------------------------------------------------------------------
  // study vs pT
  Float_t LowerPt[6];
  LowerPt[0]=10;
  LowerPt[1]=15;
  LowerPt[2]=20;  
  LowerPt[3]=25;
  LowerPt[4]=30;  
  LowerPt[5]=35;  
  TH1F *FakeableJetsPt = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)", 5, LowerPt);
  TH1F *RecoPt         = new TH1F( "RecoPt",          "reconstructed p_{T} (GeV)", 5, LowerPt);
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", 5, LowerPt);
  TH1F *RecoPtBarrel1  = new TH1F( "RecoPtBarrel1",   "reconstructed p_{T} (GeV)", 5, LowerPt);    // this is done to subsplit barrel and endcap in 2 
  TH1F *RecoPtBarrel2  = new TH1F( "RecoPtBarrel2",   "reconstructed p_{T} (GeV)", 5, LowerPt);
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", 5, LowerPt);
  TH1F *RecoPtEndcap1  = new TH1F( "RecoPtEndcap1",   "reconstructed p_{T} (GeV)", 5, LowerPt);
  TH1F *RecoPtEndcap2  = new TH1F( "RecoPtEndcap2",   "reconstructed p_{T} (GeV)", 5, LowerPt);

  std::vector<TH1F*> CutIdPt;
  std::vector<TH1F*> CutIdOnlyIDPt;
  std::vector<TH1F*> CutIdOnlyIsoPt;
  std::vector<TH1F*> CutIdOnlyConvPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt",   "cut ID #eta", 5, LowerPt );
    CutIdPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPt",   "cut ID #eta", 5, LowerPt );
    CutIdOnlyIDPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPt",  "cut ID #eta", 5, LowerPt );
    CutIdOnlyIsoPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPt", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPt.push_back(aHisto);
  }
  
  std::vector<TH1F*> LHIdPt;
  std::vector<TH1F*> LHIdOnlyIDPt;
  std::vector<TH1F*> LHIdOnlyIsoPt;
  std::vector<TH1F*> LHIdOnlyConvPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt", "like ID #eta",   5, LowerPt );
    LHIdPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPt", "like ID #eta",   5, LowerPt );
    LHIdOnlyIDPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPt", "like ID #eta",  5, LowerPt );
    LHIdOnlyIsoPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPt", "like ID #eta", 5, LowerPt );
    LHIdOnlyConvPt.push_back(aHisto);
  }

  // to have the full picture
  std::vector<TH1F*> CutIdPtBarrel;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel", "cut ID #eta",   5, LowerPt );
    CutIdPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel", "cut ID #eta",   5, LowerPt );
    CutIdOnlyIDPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel", "cut ID #eta",  5, LowerPt );
    CutIdOnlyIsoPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPtBarrel.push_back(aHisto);
  }

  // to have the full picture
  std::vector<TH1F*> LHIdPtBarrel;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel", "like ID #eta",   5, LowerPt );
    LHIdPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel", "like ID #eta",   5, LowerPt );
    LHIdOnlyIDPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel", "like ID #eta",  5, LowerPt );
    LHIdOnlyIsoPtBarrel.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel", "like ID #eta", 5, LowerPt );
    LHIdOnlyConvPtBarrel.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins
  std::vector<TH1F*> CutIdPtBarrel1;
  std::vector<TH1F*> CutIdPtBarrel2;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel1;
  std::vector<TH1F*> CutIdOnlyIDPtBarrel2;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel1;
  std::vector<TH1F*> CutIdOnlyIsoPtBarrel2;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel1;
  std::vector<TH1F*> CutIdOnlyConvPtBarrel2;

  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel1", "cut ID #eta",   5, LowerPt );
    CutIdPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel1", "cut ID #eta",   5, LowerPt );
    CutIdOnlyIDPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel1", "cut ID #eta",  5, LowerPt );
    CutIdOnlyIsoPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel1", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtBarrel2", "cut ID #eta",   5, LowerPt );
    CutIdPtBarrel2.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtBarrel2", "cut ID #eta",   5, LowerPt );
    CutIdOnlyIDPtBarrel2.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtBarrel2", "cut ID #eta",  5, LowerPt );
    CutIdOnlyIsoPtBarrel2.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtBarrel2", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPtBarrel2.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins
  std::vector<TH1F*> LHIdPtBarrel1;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel1;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel1;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel1;
  std::vector<TH1F*> LHIdPtBarrel2;
  std::vector<TH1F*> LHIdOnlyIDPtBarrel2;
  std::vector<TH1F*> LHIdOnlyIsoPtBarrel2;
  std::vector<TH1F*> LHIdOnlyConvPtBarrel2;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel1", "like ID #eta",   5, LowerPt );
    LHIdPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel1", "like ID #eta",   5, LowerPt );
    LHIdOnlyIDPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel1", "like ID #eta",  5, LowerPt );
    LHIdOnlyIsoPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel1", "like ID #eta", 5, LowerPt );
    LHIdOnlyConvPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtBarrel2", "like ID #eta",   5, LowerPt );
    LHIdPtBarrel2.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtBarrel2", "like ID #eta",   5, LowerPt );
    LHIdOnlyIDPtBarrel2.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtBarrel2", "like ID #eta",  5, LowerPt );
    LHIdOnlyIsoPtBarrel2.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtBarrel2", "like ID #eta", 5, LowerPt );
    LHIdOnlyConvPtBarrel2.push_back(aHisto);
  }

  // to have the full picture
  std::vector<TH1F*> CutIdPtEndcap;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap", "cut ID #eta",   5, LowerPt );
    CutIdPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap", "cut ID #eta",   5, LowerPt );
    CutIdOnlyIDPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap", "cut ID #eta",  5, LowerPt );
    CutIdOnlyIsoPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPtEndcap.push_back(aHisto);
  }

  // to have the full picture
  std::vector<TH1F*> LHIdPtEndcap;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap", "like ID #eta",   5, LowerPt);
    LHIdPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap", "like ID #eta",   5, LowerPt);
    LHIdOnlyIDPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap", "like ID #eta",  5, LowerPt);
    LHIdOnlyIsoPtEndcap.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap", "like ID #eta", 5, LowerPt);
    LHIdOnlyConvPtEndcap.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins
  std::vector<TH1F*> CutIdPtEndcap1;
  std::vector<TH1F*> CutIdPtEndcap2;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap1;
  std::vector<TH1F*> CutIdOnlyIDPtEndcap2;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap1;
  std::vector<TH1F*> CutIdOnlyIsoPtEndcap2;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap1;
  std::vector<TH1F*> CutIdOnlyConvPtEndcap2;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap1", "cut ID #eta",   5, LowerPt );
    CutIdPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap1", "cut ID #eta",   5, LowerPt );
    CutIdOnlyIDPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap1", "cut ID #eta",  5, LowerPt );
    CutIdOnlyIsoPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap1", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PtEndcap2", "cut ID #eta",   5, LowerPt );
    CutIdPtEndcap2.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPtEndcap2", "cut ID #eta",   5, LowerPt );
    CutIdOnlyIDPtEndcap2.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPtEndcap2", "cut ID #eta",  5, LowerPt );
    CutIdOnlyIsoPtEndcap2.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPtEndcap2", "cut ID #eta", 5, LowerPt );
    CutIdOnlyConvPtEndcap2.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins
  std::vector<TH1F*> LHIdPtEndcap1;
  std::vector<TH1F*> LHIdPtEndcap2;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap1;
  std::vector<TH1F*> LHIdOnlyIDPtEndcap2;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap1;
  std::vector<TH1F*> LHIdOnlyIsoPtEndcap2;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap1;
  std::vector<TH1F*> LHIdOnlyConvPtEndcap2;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap1", "like ID #eta",   5, LowerPt);
    LHIdPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap1", "like ID #eta",   5, LowerPt);
    LHIdOnlyIDPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap1", "like ID #eta",  5, LowerPt);
    LHIdOnlyIsoPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap1", "like ID #eta", 5, LowerPt);
    LHIdOnlyConvPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PtEndcap2", "like ID #eta",   5, LowerPt);
    LHIdPtEndcap2.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPtEndcap2", "like ID #eta",   5, LowerPt);
    LHIdOnlyIDPtEndcap2.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPtEndcap2", "like ID #eta",  5, LowerPt);
    LHIdOnlyIsoPtEndcap2.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPtEndcap2", "like ID #eta", 5, LowerPt);
    LHIdOnlyConvPtEndcap2.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs PU
  int nbinsPU = 26;
  float minPU = 0;
  float maxPU = 25;
  
  TH1F *GenPU          = new TH1F( "GenPU",          "generated nPU",     nbinsPU, minPU, maxPU );
  TH1F *RecoPU         = new TH1F( "RecoPU",         "reconstructed nPU", nbinsPU, minPU, maxPU );
  TH1F *FakeableJetsPU = new TH1F( "FakeableJetsPU", "fakeable jets nPU", nbinsPU, minPU, maxPU );
  TH1F *RecoPUHighPt   = new TH1F( "RecoPUHighPt",   "reconstructed nPU", nbinsPU, minPU, maxPU );
  TH1F *RecoPULowPt    = new TH1F( "RecoPULowPt",    "reconstructed nPU", nbinsPU, minPU, maxPU );

  std::vector<TH1F*> CutIdPU;
  std::vector<TH1F*> CutIdOnlyIDPU;
  std::vector<TH1F*> CutIdOnlyIsoPU;
  std::vector<TH1F*> CutIdOnlyConvPU;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PU", "cut ID nPU",   nbinsPU, minPU, maxPU );
    CutIdPU.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPU", "cut ID nPU",   nbinsPU, minPU, maxPU );
    CutIdOnlyIDPU.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPU", "cut ID nPU",  nbinsPU, minPU, maxPU );
    CutIdOnlyIsoPU.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPU", "cut ID nPU", nbinsPU, minPU, maxPU );
    CutIdOnlyConvPU.push_back(aHisto);
  }
  
  std::vector<TH1F*> LHIdPU;
  std::vector<TH1F*> LHIdOnlyIDPU;
  std::vector<TH1F*> LHIdOnlyIsoPU;
  std::vector<TH1F*> LHIdOnlyConvPU;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PU", "like ID nPU",   nbinsPU, minPU, maxPU );
    LHIdPU.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPU", "like ID nPU",   nbinsPU, minPU, maxPU );
    LHIdOnlyIDPU.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPU", "like ID nPU",  nbinsPU, minPU, maxPU );
    LHIdOnlyIsoPU.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPU", "like ID nPU", nbinsPU, minPU, maxPU );
    LHIdOnlyConvPU.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdPUHighPt;
  std::vector<TH1F*> CutIdOnlyIDPUHighPt;
  std::vector<TH1F*> CutIdOnlyIsoPUHighPt;
  std::vector<TH1F*> CutIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PUHighPt", "cut ID nPU",   nbinsPU, minPU, maxPU );
    CutIdPUHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPUHighPt", "cut ID nPU",   nbinsPU, minPU, maxPU );
    CutIdOnlyIDPUHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPUHighPt", "cut ID nPU",  nbinsPU, minPU, maxPU );
    CutIdOnlyIsoPUHighPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPUHighPt", "cut ID nPU", nbinsPU, minPU, maxPU );
    CutIdOnlyConvPUHighPt.push_back(aHisto);
  }
  
  std::vector<TH1F*> LHIdPUHighPt;
  std::vector<TH1F*> LHIdOnlyIDPUHighPt;
  std::vector<TH1F*> LHIdOnlyIsoPUHighPt;
  std::vector<TH1F*> LHIdOnlyConvPUHighPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PUHighPt", "like ID nPU",   nbinsPU, minPU, maxPU );
    LHIdPUHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPUHighPt", "like ID nPU",   nbinsPU, minPU, maxPU );
    LHIdOnlyIDPUHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPUHighPt", "like ID nPU",  nbinsPU, minPU, maxPU );
    LHIdOnlyIsoPUHighPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPUHighPt", "like ID nPU", nbinsPU, minPU, maxPU );
    LHIdOnlyConvPUHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> CutIdPULowPt;
  std::vector<TH1F*> CutIdOnlyIDPULowPt;
  std::vector<TH1F*> CutIdOnlyIsoPULowPt;
  std::vector<TH1F*> CutIdOnlyConvPULowPt;
  for (int i=0;i<EgammaCutBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"PULowPt", "cut ID nPU",   nbinsPU, minPU, maxPU );
    CutIdPULowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPULowPt", "cut ID nPU",   nbinsPU, minPU, maxPU );
    CutIdOnlyIDPULowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIsoPULowPt", "cut ID nPU",  nbinsPU, minPU, maxPU );
    CutIdOnlyIsoPULowPt.push_back(aHisto);
    aHisto = new TH1F( "CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyConvPULowPt", "cut ID nPU", nbinsPU, minPU, maxPU );
    CutIdOnlyConvPULowPt.push_back(aHisto);
  }

  std::vector<TH1F*> LHIdPULowPt;
  std::vector<TH1F*> LHIdOnlyIDPULowPt;
  std::vector<TH1F*> LHIdOnlyIsoPULowPt;
  std::vector<TH1F*> LHIdOnlyConvPULowPt;
  for (int i=0;i<EgammaLHBasedID.size();++i) { 
    TH1F* aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"PULowPt", "like ID nPU",   nbinsPU, minPU, maxPU );
    LHIdPULowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPULowPt", "like ID nPU",   nbinsPU, minPU, maxPU );
    LHIdOnlyIDPULowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIsoPULowPt", "like ID nPU",  nbinsPU, minPU, maxPU );
    LHIdOnlyIsoPULowPt.push_back(aHisto);
    aHisto = new TH1F( "LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyConvPULowPt", "like ID nPU", nbinsPU, minPU, maxPU );
    LHIdOnlyConvPULowPt.push_back(aHisto);
  }

  // to study the event selection
  char filename[200];
  sprintf(filename,"%sKineTree.root",outname);
  myOutTree = new FakeTree(filename);
  sprintf(filename,"%sKineTreePassed.root",outname);
  myOutTreePassed = new FakeTree(filename);
  
  // trigger: electrons - chiara
  cout << "using electrons triggers" << endl;
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v1");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v3");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v4");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v5");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v6");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v7");
  // 
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v1");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v3");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v4");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v5");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v6");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v7");
  // 
  // for mc Winter10 and for spring11
  // requiredTriggers.push_back("HLT_Ele10_SW_L1R_v2");
  // requiredTriggers.push_back("HLT_Ele17_SW_L1R_v2");


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
    
    // all events
    myCounter.IncrVar("event",1);

    // event selection: trigger
    reloadTriggerMask(true);
    bool passedHLT = hasPassedHLT();
    if ( !passedHLT ) continue;      
    myCounter.IncrVar("trigger",1);    

    // electrons passing the denominator selection to reduce W and Z contamination
    std::vector<int> denomElectrons;
    for(int iele=0; iele<nEle; iele++) {
      bool isGoodDenom = isDenomFake_HwwEgamma(iele);    
      // bool isGoodDenom = isDenomFake_smurfs(iele);
      if (!isGoodDenom) continue;
      denomElectrons.push_back(iele);
    } 
    
    // event selection: at least one candidate for denominator
    if (denomElectrons.size()==0) continue;
    myCounter.IncrVar("denom",1);    

    // denominators: kinematics
    std::pair<int,int> possibleDenom = getBestGoodElePair(denomElectrons);     
    int theDenom1(possibleDenom.first);
    int theDenom2(possibleDenom.second);
    TLorentzVector tlvDenom1, tlvDenom2;
    tlvDenom1.SetXYZT(pxEle[theDenom1],pyEle[theDenom1],pzEle[theDenom1],energyEle[theDenom1]);
    tlvDenom2.SetXYZT(pxEle[theDenom2],pyEle[theDenom2],pzEle[theDenom2],energyEle[theDenom2]);
    TVector3 tv3Denom1, tv3Denom2;
    tv3Denom1.SetXYZ(pxEle[theDenom1],pyEle[theDenom1],pzEle[theDenom1]);
    tv3Denom2.SetXYZ(pxEle[theDenom2],pyEle[theDenom2],pzEle[theDenom2]);
    if (theDenom1==-1) cout << "sanity check: impossibile!" << endl;

    // look for the leading jet not matching the HLT object - to further reduce the W contribution
    float maxEt = -1.;
    int leadingJet = -1;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      bool HLTmatch = triggerMatch(p3Jet.Eta(),p3Jet.Phi(),0.2);
      if (HLTmatch) continue;
      if ( fabs(p3Jet.Eta()) < maxEta && p3Jet.Pt() > maxEt) {   
	maxEt = p3Jet.Pt();
	leadingJet = jet;
      }
    }

    // variables used in the selection
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    TVector3 p3LeadingCut(pxAK5PFPUcorrJet[leadingJet],pyAK5PFPUcorrJet[leadingJet],pzAK5PFPUcorrJet[leadingJet]);
    float WmT = sqrt(2*tlvDenom1.Pt()*p3Met.Pt()*(1-cos(tv3Denom1.Angle(p3Met))) );      
    float theInvMass = -10.;
    if (theDenom1>-1 && theDenom2>-1) theInvMass = (tlvDenom1+tlvDenom2).M();
    float deltaPhi = fabs(p3LeadingCut.DeltaPhi(tv3Denom1));

    // fill the kine tree - after HLT and denominator
    myOutTree -> fill( p3Met.Pt(), WmT, theInvMass, maxEt, tv3Denom1.Pt(), deltaPhi);
    myOutTree -> store();
    
    // event selection: met cut to reduce W->enu
    if( p3Met.Pt() > 20 ) continue;       
    myCounter.IncrVar("met",1);    

    // event selection: mT cut to reduce W->enu [highest pT denominator candidate used]
    if (WmT > 25. ) continue;
    myCounter.IncrVar("trasvMass",1);    

    // event selection: invariant mass between two electrons passing the denominator selection to reduce Z->ee and low resonances
    if (theInvMass>60. && theInvMass<120.) continue;
    if (theInvMass>0.1 && theInvMass<12.)  continue;  
    myCounter.IncrVar("Zmass",1);    
    
    // need at least one reco jet
    if ( leadingJet < 0 ) continue;
    myCounter.IncrVar("leadingExist",1);    

    // the leading jet must be back-to-back wrt the fake candidate
    if ( deltaPhi<1 ) continue;
    myCounter.IncrVar("deltaPhi",1);  
    
    // minimal cut on leading jet or photon ET
    // if (p3LeadingCut.Pt()<15) continue;   
    // if (p3LeadingCut.Pt()<20) continue;   
    if (p3LeadingCut.Pt()<30) continue;        // chiara
    // if (p3LeadingCut.Pt()<35) continue;   
    // if (p3LeadingCut.Pt()<50) continue;   
    myCounter.IncrVar("leadingPT",1);    
    
    // fill the kine tree - after full selection  
    myOutTreePassed -> fill( p3Met.Pt(), WmT, theInvMass, maxEt, tv3Denom1.Pt(), deltaPhi);
    myOutTreePassed -> store();
    
    // consider as denominator all the reco electrons matching the HLT candidate 
    // passing some requirements on the following variables:
    // Gsf Electron (Ecal driven)
    // - ecal, hcal and racker Isolation as in the trigger
    // - H/E
    // - sigma ieta ieta
    // - deltaEta
    // - deltaPhi 
    // conversion rejection ?? 

    // weight to the events to reweight for jet pT
    float theWeight = 1.;

    // fill the denominator: take only the highest pT denominator candidate
    Utils anaUtils;
    float etaFake = fabs(tv3Denom1.Eta()); 
    float etFake  = tv3Denom1.Pt();
    bool isInEB   = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theDenom1], isEB);
    bool isInEE   = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theDenom1], isEE);
    bool highPt   = (etFake>20.);
    bool lowPt    = (etFake<=20.);
    // to split in 4 eta regions
    int etaRegion = -1;
    if (fabs(etaFake)>=0. && fabs(etaFake)<1.)    etaRegion = 1; 
    if (fabs(etaFake)>=1. && fabs(etaFake)<1.479) etaRegion = 2;
    if (fabs(etaFake)>=1.479 && fabs(etaFake)<2.) etaRegion = 3; 
    if (fabs(etaFake)>=2. && fabs(etaFake)<2.5)   etaRegion = 4; 
     
    // filling
    RecoEta         -> Fill(etaFake, theWeight);
    RecoPt          -> Fill(etFake, theWeight);
    RecoPU          -> Fill(nPV, theWeight);
    FakeableJetsEta -> Fill(etaFake, theWeight);
    FakeableJetsPt  -> Fill(etFake, theWeight);
    FakeableJetsPU  -> Fill(nPV, theWeight);
    if (highPt) { 
      RecoEtaHighPt -> Fill(etaFake, theWeight);
      RecoPUHighPt  -> Fill(nPV, theWeight);
    }
    if (lowPt)  { 
      RecoEtaLowPt -> Fill(etaFake, theWeight);
      RecoPULowPt  -> Fill(nPV, theWeight);
    }
    if (isInEB) RecoPtBarrel -> Fill(etFake, theWeight);
    if (isInEE) RecoPtEndcap -> Fill(etFake, theWeight);
    if (etaRegion==1) RecoPtBarrel1 -> Fill(etFake, theWeight);      
    if (etaRegion==2) RecoPtBarrel2 -> Fill(etFake, theWeight);      
    if (etaRegion==3) RecoPtEndcap1 -> Fill(etFake, theWeight);      
    if (etaRegion==4) RecoPtEndcap2 -> Fill(etFake, theWeight);      

    // for the systematics of HWW - init
    float dREleJet_min = 1000;
    int closestJet=-1;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 10.0 ) {          
	float dREleJet = p3Jet.DeltaR(tv3Denom1);
	if(dREleJet<dREleJet_min) {
	  closestJet=jet;
	  dREleJet_min=dREleJet;
	}}}
    if(closestJet > -1) {
      TVector3 p3ClosestJet(pxAK5PFPUcorrJet[closestJet],pyAK5PFPUcorrJet[closestJet],pzAK5PFPUcorrJet[closestJet]);
      if ( etFake>10.0 ) FakeJetForThresholdPT ->Fill(p3ClosestJet.Pt(), theWeight);
    }
    // for the systematics of HWW - end

    // does this denominator pass the IP cut as for H->WW ?  // chiara, hardcoded
    bool isDenomIP = true;
    int gsfTrack = gsfTrackIndexEle[theDenom1];
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[gsfTrack];   
    if ( fabs(dxyEle)>0.02 ) isDenomIP = false;
    if ( fabs(dzEle)>0.10 )  isDenomIP = false;
    // does this denominator pass the IP cut as for H->WW ?  
    
    // fill the numerator: simple cut based
    for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut) {
      
      bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
      isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
      isEleID(&EgammaCutBasedID[icut],theDenom1,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

      // for the smurf selection, add further cut for pT<20 GeV --> move to WP70
      // this way for 5 and 6 we superimpose in any case the 6 at low pT
      if(TString(EgammaCutBasedIDWPs[icut].c_str()).Contains("Smurf") && etFake<20.) {
        isEleID(&EgammaCutBasedID[6],theDenom1,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
        Utils anaUtils;
        int smurfsc;
        bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theDenom1], bits::isEcalDriven);
        float scEta = -1.;
        if ( ecalDriven) { smurfsc = superClusterIndexEle[theDenom1];   scEta = etaSC[smurfsc]; }
        if (!ecalDriven) { smurfsc = PFsuperClusterIndexEle[theDenom1]; scEta = etaPFSC[smurfsc]; }
	isEleIDCutBased = isEleIDCutBased && (fbremEle[theDenom1]>0.15 || (fabs(scEta)<1.0 && eSuperClusterOverPEle[theDenom1]>0.95));   
      }

      if ( isEleIDCutBased ) {
	CutIdOnlyIDEta[icut]-> Fill(etaFake, theWeight);
	CutIdOnlyIDPt[icut] -> Fill(etFake, theWeight);
	CutIdOnlyIDPU[icut] -> Fill(nPV, theWeight);
	if (highPt) { 
	  CutIdOnlyIDEtaHighPt[icut]-> Fill(etaFake, theWeight);
	  CutIdOnlyIDPUHighPt[icut] -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  CutIdOnlyIDEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  CutIdOnlyIDPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) CutIdOnlyIDPtBarrel[icut]  -> Fill(etFake, theWeight);
	if (isInEE) CutIdOnlyIDPtEndcap[icut]  -> Fill(etFake, theWeight);
	if (etaRegion==1) CutIdOnlyIDPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) CutIdOnlyIDPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) CutIdOnlyIDPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) CutIdOnlyIDPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
      
      if ( isIsolCutBased ) {
	CutIdOnlyIsoEta[icut]-> Fill(etaFake, theWeight);
	CutIdOnlyIsoPt[icut] -> Fill(etFake, theWeight);
	CutIdOnlyIsoPU[icut] -> Fill(nPV, theWeight);
	if (highPt) { 
	  CutIdOnlyIsoEtaHighPt[icut] -> Fill(etaFake, theWeight);
	  CutIdOnlyIsoPUHighPt[icut]  -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  CutIdOnlyIsoEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  CutIdOnlyIsoPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) CutIdOnlyIsoPtBarrel[icut]  -> Fill(etFake, theWeight);
	if (isInEE) CutIdOnlyIsoPtEndcap[icut]  -> Fill(etFake, theWeight);
	if (etaRegion==1) CutIdOnlyIsoPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) CutIdOnlyIsoPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) CutIdOnlyIsoPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) CutIdOnlyIsoPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
      
      if ( isConvRejCutBased ) {
	CutIdOnlyConvEta[icut]-> Fill(etaFake, theWeight);
	CutIdOnlyConvPt[icut] -> Fill(etFake, theWeight);
	CutIdOnlyConvPU[icut] -> Fill(nPV, theWeight);
	if (highPt) { 
	  CutIdOnlyConvEtaHighPt[icut] -> Fill(etaFake, theWeight);
	  CutIdOnlyConvPUHighPt[icut]  -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  CutIdOnlyConvEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  CutIdOnlyConvPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) CutIdOnlyConvPtBarrel[icut] ->Fill(etFake, theWeight);
	if (isInEE) CutIdOnlyConvPtEndcap[icut] ->Fill(etFake, theWeight);
	if (etaRegion==1) CutIdOnlyConvPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) CutIdOnlyConvPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) CutIdOnlyConvPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) CutIdOnlyConvPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
      
      if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased && isDenomIP) {  // chiara: IP cut added
	CutIdEta[icut]->Fill(etaFake, theWeight);
	CutIdPt[icut] ->Fill(etFake, theWeight);
	CutIdPU[icut] ->Fill(nPV, theWeight);
	if (highPt) { 
	  CutIdEtaHighPt[icut]-> Fill(etaFake, theWeight);
	  CutIdPUHighPt[icut] -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  CutIdEtaLowPt[icut]-> Fill(etaFake, theWeight);
	  CutIdPULowPt[icut] -> Fill(nPV, theWeight);
	}
	if (isInEB) CutIdPtBarrel[icut]-> Fill(etFake, theWeight);
	if (isInEE) CutIdPtEndcap[icut]-> Fill(etFake, theWeight);
	if (etaRegion==1) CutIdPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) CutIdPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) CutIdPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) CutIdPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
    }
    
    // numerator: likelihood based
    for (int icut=0;icut<EgammaLHBasedIDWPs.size();++icut) {

      bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
      isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
      isEleID(&EgammaLHBasedID[icut],theDenom1,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
      
      // to combine high pT and low pT likelihood
      // this way for 5 and 6 we superimpose in any case the 6 at low pT
      if(TString(EgammaLHBasedIDWPs[icut].c_str()).Contains("PFIso") && etFake<20.) 
        isEleID(&EgammaLHBasedID[6],theDenom1,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);    // chiara: ema e' quel che vuoi?

      if ( isEleIDCutBased ) {
	LHIdOnlyIDEta[icut]->Fill(etaFake, theWeight);
	LHIdOnlyIDPt[icut] ->Fill(etFake, theWeight);
	LHIdOnlyIDPU[icut] ->Fill(nPV, theWeight);
	if (highPt) { 
	  LHIdOnlyIDEtaHighPt[icut]-> Fill(etaFake, theWeight);
	  LHIdOnlyIDPUHighPt[icut] -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  LHIdOnlyIDEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  LHIdOnlyIDPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) LHIdOnlyIDPtBarrel[icut] -> Fill(etFake, theWeight);
	if (isInEE) LHIdOnlyIDPtEndcap[icut] -> Fill(etFake, theWeight);
	if (etaRegion==1) LHIdOnlyIDPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) LHIdOnlyIDPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) LHIdOnlyIDPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) LHIdOnlyIDPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
      
      if ( isIsolCutBased ) {
	LHIdOnlyIsoEta[icut]-> Fill(etaFake, theWeight);
	LHIdOnlyIsoPt[icut] -> Fill(etFake, theWeight);
	LHIdOnlyIsoPU[icut] -> Fill(nPV, theWeight);
	if (highPt) { 
	  LHIdOnlyIsoEtaHighPt[icut]-> Fill(etaFake, theWeight);
	  LHIdOnlyIsoPUHighPt[icut] -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  LHIdOnlyIsoEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  LHIdOnlyIsoPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) LHIdOnlyIsoPtBarrel[icut] -> Fill(etFake, theWeight);
	if (isInEE) LHIdOnlyIsoPtEndcap[icut] -> Fill(etFake, theWeight);
	if (etaRegion==1) LHIdOnlyIsoPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) LHIdOnlyIsoPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) LHIdOnlyIsoPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) LHIdOnlyIsoPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
      
      if ( isConvRejCutBased ) {
	LHIdOnlyConvEta[icut]->Fill(etaFake, theWeight);
	LHIdOnlyConvPt[icut] ->Fill(etFake, theWeight);
	LHIdOnlyConvPU[icut] ->Fill(nPV, theWeight);
	if (highPt) { 
	  LHIdOnlyConvEtaHighPt[icut]-> Fill(etaFake, theWeight);
	  LHIdOnlyConvPUHighPt[icut] -> Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  LHIdOnlyConvEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  LHIdOnlyConvPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) LHIdOnlyConvPtBarrel[icut] -> Fill(etFake, theWeight);
	if (isInEE) LHIdOnlyConvPtEndcap[icut] -> Fill(etFake, theWeight);
	if (etaRegion==1) LHIdOnlyConvPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) LHIdOnlyConvPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) LHIdOnlyConvPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) LHIdOnlyConvPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
      
      if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased && isDenomIP) {  // chiara: IP added
	LHIdEta[icut]->Fill(etaFake, theWeight);
	LHIdPt[icut] ->Fill(etFake, theWeight);
	LHIdPU[icut] ->Fill(nPV, theWeight);
	if (highPt) { 
	  LHIdEtaHighPt[icut]->Fill(etaFake, theWeight);
	  LHIdPUHighPt[icut] ->Fill(nPV, theWeight);
	}
	if (lowPt)  { 
	  LHIdEtaLowPt[icut] -> Fill(etaFake, theWeight);
	  LHIdPULowPt[icut]  -> Fill(nPV, theWeight);
	}
	if (isInEB) LHIdPtBarrel[icut] ->Fill(etFake, theWeight);
	if (isInEE) LHIdPtEndcap[icut] ->Fill(etFake, theWeight);
	if (etaRegion==1) LHIdPtBarrel1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==2) LHIdPtBarrel2[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==3) LHIdPtEndcap1[icut] -> Fill(etFake, theWeight);      
	if (etaRegion==4) LHIdPtEndcap2[icut] -> Fill(etFake, theWeight);      
      }
    }
  } // loop events

  // saving the counters
  sprintf(filename,"%sCounters.root",outname);
  myCounter.Save(filename,"recreate");
  
  // saving the output tree
  myOutTree       -> save();
  myOutTreePassed -> save();

  // saving histo to choose the ET cut 
  sprintf(filename,"%s-EleMisidPt-Threshold.root",outname);
  TFile fileThreshold(filename,"RECREATE");
  fileThreshold.cd();
  FakeJetForThresholdPT->Write();

  // saving efficiency histos
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  ElectronEffEta.AddNumerator(FakeableJetsEta);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffEta.AddNumerator(CutIdEta[icut]);
    ElectronEffEta.AddNumerator(CutIdOnlyIDEta[icut]);
    ElectronEffEta.AddNumerator(CutIdOnlyIsoEta[icut]);
    ElectronEffEta.AddNumerator(CutIdOnlyConvEta[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut){
    ElectronEffEta.AddNumerator(LHIdEta[icut]);
    ElectronEffEta.AddNumerator(LHIdOnlyIDEta[icut]);
    ElectronEffEta.AddNumerator(LHIdOnlyIsoEta[icut]);
    ElectronEffEta.AddNumerator(LHIdOnlyConvEta[icut]);
  }
  ElectronEffEta.SetDenominator(FakeableJetsEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("fake rate vs #eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("Fake rate");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleMisidEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffEtaHighPt.AddNumerator(CutIdEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIDEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CutIdOnlyIsoEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(CutIdOnlyConvEtaHighPt[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffEtaHighPt.AddNumerator(LHIdEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIDEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(LHIdOnlyIsoEtaHighPt[icut]);
    ElectronEffEtaHighPt.AddNumerator(LHIdOnlyConvEtaHighPt[icut]);
  }
  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaHighPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffEtaLowPt.AddNumerator(CutIdEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIDEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CutIdOnlyIsoEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(CutIdOnlyConvEtaLowPt[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffEtaLowPt.AddNumerator(LHIdEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIDEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(LHIdOnlyIsoEtaLowPt[icut]);
    ElectronEffEtaLowPt.AddNumerator(LHIdOnlyConvEtaLowPt[icut]);
  }
  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("fake rate vs #eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  ElectronEffPt.AddNumerator(FakeableJetsPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPt.AddNumerator(CutIdPt[icut]);
    ElectronEffPt.AddNumerator(CutIdOnlyIDPt[icut]);
    ElectronEffPt.AddNumerator(CutIdOnlyIsoPt[icut]);
    ElectronEffPt.AddNumerator(CutIdOnlyConvPt[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPt.AddNumerator(LHIdPt[icut]);
    ElectronEffPt.AddNumerator(LHIdOnlyIDPt[icut]);
    ElectronEffPt.AddNumerator(LHIdOnlyIsoPt[icut]);
    ElectronEffPt.AddNumerator(LHIdOnlyConvPt[icut]);
  }
  ElectronEffPt.SetDenominator(FakeableJetsPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("fake rate vs p_{T}");
  ElectronEffPt.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPt.SetYaxisTitle("Fake rate");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleMisidPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtBarrel.AddNumerator(CutIdPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CutIdOnlyIDPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CutIdOnlyIsoPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(CutIdOnlyConvPtBarrel[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtBarrel.AddNumerator(LHIdPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(LHIdOnlyIDPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(LHIdOnlyIsoPtBarrel[icut]);
    ElectronEffPtBarrel.AddNumerator(LHIdOnlyConvPtBarrel[icut]);
  }
  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("fake rate vs p_{T}");
  ElectronEffPtBarrel.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-EleMisidPtBarrel1.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel1(filename);
  ElectronEffPtBarrel1.AddNumerator(RecoPtBarrel1);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtBarrel1.AddNumerator(CutIdPtBarrel1[icut]);
    ElectronEffPtBarrel1.AddNumerator(CutIdOnlyIDPtBarrel1[icut]);
    ElectronEffPtBarrel1.AddNumerator(CutIdOnlyIsoPtBarrel1[icut]);
    ElectronEffPtBarrel1.AddNumerator(CutIdOnlyConvPtBarrel1[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtBarrel1.AddNumerator(LHIdPtBarrel1[icut]);
    ElectronEffPtBarrel1.AddNumerator(LHIdOnlyIDPtBarrel1[icut]);
    ElectronEffPtBarrel1.AddNumerator(LHIdOnlyIsoPtBarrel1[icut]);
    ElectronEffPtBarrel1.AddNumerator(LHIdOnlyConvPtBarrel1[icut]);
  }
  ElectronEffPtBarrel1.SetDenominator(RecoPtBarrel1);
  ElectronEffPtBarrel1.ComputeEfficiencies();
  ElectronEffPtBarrel1.SetTitle("fake rate vs p_{T}");
  ElectronEffPtBarrel1.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel1.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel1.SetYaxisMin(0.0);
  ElectronEffPtBarrel1.Write();

  sprintf(filename,"%s-EleMisidPtBarrel2.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel2(filename);
  ElectronEffPtBarrel2.AddNumerator(RecoPtBarrel2);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtBarrel2.AddNumerator(CutIdPtBarrel2[icut]);
    ElectronEffPtBarrel2.AddNumerator(CutIdOnlyIDPtBarrel2[icut]);
    ElectronEffPtBarrel2.AddNumerator(CutIdOnlyIsoPtBarrel2[icut]);
    ElectronEffPtBarrel2.AddNumerator(CutIdOnlyConvPtBarrel2[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtBarrel2.AddNumerator(LHIdPtBarrel2[icut]);
    ElectronEffPtBarrel2.AddNumerator(LHIdOnlyIDPtBarrel2[icut]);
    ElectronEffPtBarrel2.AddNumerator(LHIdOnlyIsoPtBarrel2[icut]);
    ElectronEffPtBarrel2.AddNumerator(LHIdOnlyConvPtBarrel2[icut]);
  }
  ElectronEffPtBarrel2.SetDenominator(RecoPtBarrel2);
  ElectronEffPtBarrel2.ComputeEfficiencies();
  ElectronEffPtBarrel2.SetTitle("fake rate vs p_{T}");
  ElectronEffPtBarrel2.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtBarrel2.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel2.SetYaxisMin(0.0);
  ElectronEffPtBarrel2.Write();

  sprintf(filename,"%s-EleMisidPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtEndcap.AddNumerator(CutIdPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CutIdOnlyIDPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CutIdOnlyIsoPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(CutIdOnlyConvPtEndcap[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtEndcap.AddNumerator(LHIdPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(LHIdOnlyIDPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(LHIdOnlyIsoPtEndcap[icut]);
    ElectronEffPtEndcap.AddNumerator(LHIdOnlyConvPtEndcap[icut]);
  }
  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("fake rate vs p_{T}");
  ElectronEffPtEndcap.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();

  sprintf(filename,"%s-EleMisidPtEndcap1.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap1(filename);
  ElectronEffPtEndcap1.AddNumerator(RecoPtEndcap1);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtEndcap1.AddNumerator(CutIdPtEndcap1[icut]);
    ElectronEffPtEndcap1.AddNumerator(CutIdOnlyIDPtEndcap1[icut]);
    ElectronEffPtEndcap1.AddNumerator(CutIdOnlyIsoPtEndcap1[icut]);
    ElectronEffPtEndcap1.AddNumerator(CutIdOnlyConvPtEndcap1[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtEndcap1.AddNumerator(LHIdPtEndcap1[icut]);
    ElectronEffPtEndcap1.AddNumerator(LHIdOnlyIDPtEndcap1[icut]);
    ElectronEffPtEndcap1.AddNumerator(LHIdOnlyIsoPtEndcap1[icut]);
    ElectronEffPtEndcap1.AddNumerator(LHIdOnlyConvPtEndcap1[icut]);
  }
  ElectronEffPtEndcap1.SetDenominator(RecoPtEndcap1);
  ElectronEffPtEndcap1.ComputeEfficiencies();
  ElectronEffPtEndcap1.SetTitle("fake rate vs p_{T}");
  ElectronEffPtEndcap1.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap1.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap1.SetYaxisMin(0.0);
  ElectronEffPtEndcap1.Write();

  sprintf(filename,"%s-EleMisidPtEndcap2.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap2(filename);
  ElectronEffPtEndcap2.AddNumerator(RecoPtEndcap2);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut) {
    ElectronEffPtEndcap2.AddNumerator(CutIdPtEndcap2[icut]);
    ElectronEffPtEndcap2.AddNumerator(CutIdOnlyIDPtEndcap2[icut]);
    ElectronEffPtEndcap2.AddNumerator(CutIdOnlyIsoPtEndcap2[icut]);
    ElectronEffPtEndcap2.AddNumerator(CutIdOnlyConvPtEndcap2[icut]);
  }
  for (int icut=0;icut<EgammaLHBasedID.size();++icut) {
    ElectronEffPtEndcap2.AddNumerator(LHIdPtEndcap2[icut]);
    ElectronEffPtEndcap2.AddNumerator(LHIdOnlyIDPtEndcap2[icut]);
    ElectronEffPtEndcap2.AddNumerator(LHIdOnlyIsoPtEndcap2[icut]);
    ElectronEffPtEndcap2.AddNumerator(LHIdOnlyConvPtEndcap2[icut]);
  }
  ElectronEffPtEndcap2.SetDenominator(RecoPtEndcap2);
  ElectronEffPtEndcap2.ComputeEfficiencies();
  ElectronEffPtEndcap2.SetTitle("fake rate vs p_{T}");
  ElectronEffPtEndcap2.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtEndcap2.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap2.SetYaxisMin(0.0);
  ElectronEffPtEndcap2.Write();
}



void LikelihoodAnalysis::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxPFSC[sc];
      e4SwissCross = e4SwissCrossPFSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagPFSC[sc];
      seedTime = timePFSC[sc];
      seedChi2 = chi2PFSC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  bool isEleEB= anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex], isEB);
  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->applyElectronIDOnPFlowElectrons(true);
  selector->SetHOverE( HoE );
  selector->SetS9S25( s9s25 );
  selector->SetNBrem( nbremsEle[eleIndex] );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetDPhiOut( dphiout );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetSigmaPhiPhi( spp );
  selector->SetEOverPout( eopout );
  selector->SetEOverPin( eop );
  selector->SetElectronClass ( classificationEle[eleIndex] );
  selector->SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  selector->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  selector->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetTrkIsolation ( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );

  float combinedIso = 0.0;
  if (isEleEB) combinedIso = dr03TkSumPtEle[eleIndex] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else combinedIso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  selector->SetCombinedIsolation( (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / pEle.Pt() ); 

  selector->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex]) / pEle.Pt() );

  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );
  selector->SetHasMatchedConversion ( hasMatchedConversionEle[eleIndex] );

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();
  
}

void LikelihoodAnalysis::isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eseedopin, eop, eseed,  sce, scet,sceta, rawe;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
    eseed = seedEnergySC[sc];
    sce = energySC[sc];
    rawe = rawEnergySC[sc];
    scet = energySC[sc]/TMath::CosH(etaSC[sc]);
    sceta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxPFSC[sc];
      e4SwissCross = e4SwissCrossPFSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagPFSC[sc];
      seedTime = timePFSC[sc];
      seedChi2 = chi2PFSC[sc];
      eseed = seedEnergyPFSC[sc];
      sce = energyPFSC[sc];
      rawe = rawEnergyPFSC[sc];
      scet = energyPFSC[sc]/TMath::CosH(etaPFSC[sc]);
      sceta = etaPFSC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  eseedopin = eSeedOverPoutEle[eleIndex]*(1-fbremEle[eleIndex]);
  
  selector->reset();
  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->SetSCEt(scet);
  selector->SetSCEta(sceta);
  selector->SetHOverE( HoE );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetEOverPin( eop );
  selector->SetESeedOverPin( eseedopin );
  // selector->SetLikelihood( eleIdLikelihoodEle[eleIndex] );

  selector->SetEcalIsolation( dr04EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4 );
  selector->SetTrkIsolation( dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3 );
  selector->SetHcalIsolation( dr04HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4);
  // selector->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex]) / pt );

  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );
  // selector->SetHasMatchedConversion ( hasMatchedConversionEle[eleIndex] );
  
  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputEleId();
  *isolOutput = selector->outputIso();
  *convRejOutput = selector->outputConv();

}

/*
/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float LikelihoodAnalysis::SigmaiEiE(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

/// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
float LikelihoodAnalysis::SigmaiPiP(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}
*/


// two highest pT electrons - passing the denominator selection  
std::pair<int,int> LikelihoodAnalysis::getBestGoodElePair(std::vector<int> goodElectrons) {
  
  int theEle1=-1;
  int theEle2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;
  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
    float thisPt=pEle.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theEle2 = theEle1; theEle1 = eleIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theEle2 = eleIndex; }
  }
  return make_pair(theEle1,theEle2);
}

// denominator for fake rate: for HtoWW, egamma triggers
bool LikelihoodAnalysis::isDenomFake_HwwEgamma(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);

  // match with the HLT firing candidates
  bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  if (!HLTmatch) isGoodDenom = false;
  
  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // barrel or endcap
  bool isEleEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theEle], isEB);    

  // taking shower shape                                                                                                             
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  if (ecalDriven) {
    sc = superClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]);
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]);
  }
  if ( sc<0 ) { isGoodDenom = false; }

  // sigmaIetaIeta                                                                                                                   
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; }

  // H/E
  if ( isEleEB && hOverEEle[theEle]>0.12) isGoodDenom = false;   
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) isGoodDenom = false;
  
  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) isGoodDenom = false;

  // isolation 
  float ecalIsol    = (dr03EcalRecHitSumEtEle[theEle])/p3Ele.Pt();
  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if(ecalIsol>0.2)    isGoodDenom = false;
  if(hcalIsol>0.2)    isGoodDenom = false;
  if(trackerIsol>0.2) isGoodDenom = false;                                
      
  return isGoodDenom;
}

// denominator for fake rate: for HtoWW, egamma triggers, same as smurfs
bool LikelihoodAnalysis::isDenomFake_smurfs(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // match with the HLT firing candidates
  // bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  // if (!HLTmatch) isGoodDenom = false;
  
  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // taking shower shape                                                                                                             
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  float scEta = -1.;                                                               
  if (ecalDriven) {
    sc = superClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]);
    scEta = etaSC[sc];
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]);
    scEta = etaPFSC[sc];
  }
  if ( sc<0 ) { isGoodDenom = false; }

  // barrel or endcap
  bool isEleEB = false;
  if (fabs(scEta)<1.479) isEleEB = true;   
  
  // sigmaIetaIeta                                                                                                                   
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; }

  // isolation
  float ecalIsolAbs = 0.0;
  if ( isEleEB ) ecalIsolAbs = max(0.0,dr03EcalRecHitSumEtEle[theEle]-1.0);
  else ecalIsolAbs = dr03EcalRecHitSumEtEle[theEle];
  float ecalIsol = ecalIsolAbs/p3Ele.Pt(); 
  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if(ecalIsol>0.2)    isGoodDenom = false;
  if(hcalIsol>0.2)    isGoodDenom = false;
  if(trackerIsol>0.2) isGoodDenom = false;                                

  // H/E
  if ( isEleEB && hOverEEle[theEle]>0.12) isGoodDenom = false;   
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) isGoodDenom = false;
  
  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) isGoodDenom = false;
  
  // full conversion rejection 
  int gsf = gsfTrackIndexEle[theEle];
  int missHits = expInnerLayersGsfTrack[gsf];
  bool matchConv = hasMatchedConversionEle[theEle];
  if (missHits>0 || matchConv) isGoodDenom = false;

  // impact parameter cuts 
  float dxyEle = transvImpactParGsfTrack[gsf];
  float dzEle  = PVzPV[0] - trackVzGsfTrack[gsf];
  if (fabs(dxyEle)>0.02) isGoodDenom = false;
  if (fabs(dzEle)>0.10)  isGoodDenom = false;

  return isGoodDenom;
}

