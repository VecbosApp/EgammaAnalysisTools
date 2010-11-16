#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"

using namespace bits;
using namespace std;

LikelihoodAnalysis::LikelihoodAnalysis(TTree *tree)
  : Egamma(tree) {

  _isData = false;

  EgammaCutBasedIDWPs.push_back("WP95");
  EgammaCutBasedIDWPs.push_back("WP90");
  EgammaCutBasedIDWPs.push_back("WP80");
  EgammaCutBasedIDWPs.push_back("WP70");

  EgammaCiCBasedIDWPs.push_back("CiCVeryLoose");
  EgammaCiCBasedIDWPs.push_back("CiCLoose");
  EgammaCiCBasedIDWPs.push_back("CiCMedium");
  EgammaCiCBasedIDWPs.push_back("CiCTight");
  EgammaCiCBasedIDWPs.push_back("CiCSuperTight");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight2");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight3");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight4");


  EgammaLHBasedIDWPs.push_back("LHVeryLoose");
  EgammaLHBasedIDWPs.push_back("LHLoose");
  EgammaLHBasedIDWPs.push_back("LHMedium");
  EgammaLHBasedIDWPs.push_back("LHTight");
  EgammaLHBasedIDWPs.push_back("LHHyperTight");


  // single electron efficiency
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

  // single electron efficiency
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
  TDirectory *EBlt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EBgt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;

  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;

  LH = new ElectronLikelihood(&(*EBlt15dir), &(*EElt15dir), &(*EBgt15dir), &(*EEgt15dir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);
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
//   TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
//   TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
//   TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );

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
//   TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
//   TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );


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

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  //  Long64_t nentries = 10000;
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

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
      if (highPt) RecoEtaHighPt->Fill(mcEta);
      if (lowPt) RecoEtaLowPt->Fill(mcEta);
      if (isInEB) RecoPtBarrel->Fill(mcPt);
      if (isInEE) RecoPtEndcap->Fill(mcPt);
      // Golden means 0 brem clusters, showering otherwise
//       if ( nbremsEle[matchedRecoEle] == 0 ) {
//         GoldenEta->Fill(mcEta);
//         GoldenPt->Fill(mcPt);
//       } else {
//         ShoweringEta->Fill(mcEta);
//         ShoweringPt->Fill(mcPt);
//       }

      for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCutBasedID[icut],matchedRecoEle,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

	  if ( isEleIDCutBased ) {
	    CutIdOnlyIDEta[icut]->Fill(mcEta);
	    CutIdOnlyIDPt[icut]->Fill(mcPt);
	    if (highPt) CutIdOnlyIDEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CutIdOnlyIDEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) CutIdOnlyIDPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdOnlyIDPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isIsolCutBased ) {
	    CutIdOnlyIsoEta[icut]->Fill(mcEta);
	    CutIdOnlyIsoPt[icut]->Fill(mcPt);
	    if (highPt) CutIdOnlyIsoEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CutIdOnlyIsoEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) CutIdOnlyIsoPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdOnlyIsoPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isConvRejCutBased ) {
	    CutIdOnlyConvEta[icut]->Fill(mcEta);
	    CutIdOnlyConvPt[icut]->Fill(mcPt);
	    if (highPt) CutIdOnlyConvEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CutIdOnlyConvEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) CutIdOnlyConvPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CutIdOnlyConvPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CutIdEta[icut]->Fill(mcEta);
	    CutIdPt[icut]->Fill(mcPt);
	    if (highPt) CutIdEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CutIdEtaLowPt[icut]->Fill(mcEta);
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
	    if (highPt) LHIdOnlyIDEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) LHIdOnlyIDEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) LHIdOnlyIDPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdOnlyIDPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isIsolCutBased ) {
	    LHIdOnlyIsoEta[icut]->Fill(mcEta);
	    LHIdOnlyIsoPt[icut]->Fill(mcPt);
	    if (highPt) LHIdOnlyIsoEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) LHIdOnlyIsoEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) LHIdOnlyIsoPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdOnlyIsoPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isConvRejCutBased ) {
	    LHIdOnlyConvEta[icut]->Fill(mcEta);
	    LHIdOnlyConvPt[icut]->Fill(mcPt);
	    if (highPt) LHIdOnlyConvEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) LHIdOnlyConvEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) LHIdOnlyConvPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) LHIdOnlyConvPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    LHIdEta[icut]->Fill(mcEta);
	    LHIdPt[icut]->Fill(mcPt);
	    if (highPt) LHIdEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) LHIdEtaLowPt[icut]->Fill(mcEta);
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
	    if (highPt) CiCIdOnlyIDEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CiCIdOnlyIDEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) CiCIdOnlyIDPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdOnlyIDPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isIsolCutBased ) {
	    CiCIdOnlyIsoEta[icut]->Fill(mcEta);
	    CiCIdOnlyIsoPt[icut]->Fill(mcPt);
	    if (highPt) CiCIdOnlyIsoEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CiCIdOnlyIsoEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) CiCIdOnlyIsoPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdOnlyIsoPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isConvRejCutBased ) {
	    CiCIdOnlyConvEta[icut]->Fill(mcEta);
	    CiCIdOnlyConvPt[icut]->Fill(mcPt);
	    if (highPt) CiCIdOnlyConvEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CiCIdOnlyConvEtaLowPt[icut]->Fill(mcEta);
	    if (isInEB) CiCIdOnlyConvPtBarrel[icut]->Fill(mcPt);
	    if (isInEE) CiCIdOnlyConvPtEndcap[icut]->Fill(mcPt);
	  }
	  if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
	    CiCIdEta[icut]->Fill(mcEta);
	    CiCIdPt[icut]->Fill(mcPt);
	    if (highPt) CiCIdEtaHighPt[icut]->Fill(mcEta);
	    if (lowPt) CiCIdEtaLowPt[icut]->Fill(mcEta);
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
//   ElectronEffEta.AddNumerator(GoldenEta);
//   ElectronEffEta.AddNumerator(ShoweringEta);
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
//   ElectronEffEtaHighPt.AddNumerator(GoldenEtaHighPt);
//   ElectronEffEtaHighPt.AddNumerator(ShoweringEtaHighPt);
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
//   ElectronEffEtaLowPt.AddNumerator(GoldenEtaLowPt);
//   ElectronEffEtaLowPt.AddNumerator(ShoweringEtaLowPt);
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
//   ElectronEffPt.AddNumerator(GoldenPt);
//   ElectronEffPt.AddNumerator(ShoweringPt);
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
//   ElectronEffPtBarrel.AddNumerator(GoldenPtBarrel);
//   ElectronEffPtBarrel.AddNumerator(ShoweringPtBarrel);
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
//   ElectronEffPtEndcap.AddNumerator(GoldenPtEndcap);
//   ElectronEffPtEndcap.AddNumerator(ShoweringPtEndcap);
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

}



void LikelihoodAnalysis::estimateFakeRate(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *FakeableJetsEta   = new TH1F( "FakeableJetsEta",  "fakeable jets #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
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
  TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
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

//   TH1F *CutIdWP80Pt = new TH1F( "CutIdWP80Pt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
//   TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH ID tight p_{T} (GeV)",        nbinsPt, minPt, maxPt );
//   TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH ID loose p_{T} (GeV)",        nbinsPt, minPt, maxPt );

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

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

    for ( int jet=0; jet<nAK5PFJet; jet++ ) {

      TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);

      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 20.0 ) {
        
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

    for ( int ele=0; ele<nEle; ele++ ) {

      TVector3 p3Ele(pxEle[ele], pyEle[ele], pzEle[ele]);
      
      if ( fabs(etaEle[ele]) < 2.5 && p3Ele.Pt() > 10.0 ) {

      
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Ele.DeltaR(mcEle);
        
        if(deltaR<0.3) continue;

        float dREleJet_min = 1000;
        int closestJet=-1;

        for ( int jet=0; jet<nAK5PFJet; jet++ ) {

          TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);

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

          TVector3 p3ClosestJet(pxAK5PFJet[closestJet],pyAK5PFJet[closestJet],pzAK5PFJet[closestJet]);

	  //!!//
          float etFake=p3ClosestJet.Pt();
          float etaFake=p3ClosestJet.Eta();

          RecoEta->Fill(etaFake);
          RecoPt->Fill(etFake);

          if ( nbremsEle[ele] == 0 ) {
            GoldenEta->Fill(etaFake);
            GoldenPt->Fill(etFake);
          } else {
            ShoweringEta->Fill(etaFake);
            ShoweringPt->Fill(etFake);
          }

	  for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
	    {
	      bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	      isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	      isEleID(&EgammaCutBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
	      
	      if ( isEleIDCutBased ) {
		CutIdOnlyIDEta[icut]->Fill(etaFake);
		CutIdOnlyIDPt[icut]->Fill(etFake);
	      }
	      if ( isIsolCutBased ) {
		CutIdOnlyIsoEta[icut]->Fill(etaFake);
		CutIdOnlyIsoPt[icut]->Fill(etFake);
	      }
	      if ( isConvRejCutBased ) {
		CutIdOnlyConvEta[icut]->Fill(etaFake);
		CutIdOnlyConvPt[icut]->Fill(etFake);
	      }
	      if ( isEleIDCutBased && isIsolCutBased && isConvRejCutBased ) {
		CutIdEta[icut]->Fill(etaFake);
		CutIdPt[icut]->Fill(etFake);
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
  EfficiencyEvaluator ElectronFakeRateEta(filename);
  //  ElectronFakeRateEta.AddNumerator(FakeableJetsEta);
  ElectronFakeRateEta.AddNumerator(GoldenEta);
  ElectronFakeRateEta.AddNumerator(ShoweringEta);
  ElectronFakeRateEta.AddNumerator(RecoEta);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronFakeRateEta.AddNumerator(CutIdEta[icut]);
      ElectronFakeRateEta.AddNumerator(CutIdOnlyIDEta[icut]);
      ElectronFakeRateEta.AddNumerator(CutIdOnlyIsoEta[icut]);
      ElectronFakeRateEta.AddNumerator(CutIdOnlyConvEta[icut]);
    }

//   ElectronFakeRateEta.AddNumerator(CutIdWP80Eta);
//   ElectronFakeRateEta.AddNumerator(CutIdTightCICEta);
//   ElectronFakeRateEta.AddNumerator(CutIdSuperTightCICEta);
//   ElectronFakeRateEta.AddNumerator(LHIdLooseEta);
//   ElectronFakeRateEta.AddNumerator(LHIdTightEta);
  ElectronFakeRateEta.SetDenominator(RecoEta);
  ElectronFakeRateEta.ComputeEfficiencies();
  ElectronFakeRateEta.SetTitle("jet fake probability vs #eta");
  ElectronFakeRateEta.SetXaxisTitle("#eta of closest jet");
  ElectronFakeRateEta.SetYaxisTitle("jet #rightarrow fake e probability");
  ElectronFakeRateEta.SetYaxisMin(0.0);
  ElectronFakeRateEta.Write();

  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronFakeRatePt(filename);
  //  ElectronFakeRatePt.AddNumerator(FakeableJetsPt);
  ElectronFakeRatePt.AddNumerator(GoldenPt);
  ElectronFakeRatePt.AddNumerator(ShoweringPt);
  ElectronFakeRatePt.AddNumerator(RecoPt);
  for (int icut=0;icut<EgammaCutBasedID.size();++icut)
    {
      ElectronFakeRatePt.AddNumerator(CutIdPt[icut]);
      ElectronFakeRatePt.AddNumerator(CutIdOnlyIDPt[icut]);
      ElectronFakeRatePt.AddNumerator(CutIdOnlyIsoPt[icut]);
      ElectronFakeRatePt.AddNumerator(CutIdOnlyConvPt[icut]);
    }
//   ElectronFakeRatePt.AddNumerator(CutIdWP80Pt);
//   ElectronFakeRatePt.AddNumerator(CutIdTightCICPt);
//   ElectronFakeRatePt.AddNumerator(CutIdSuperTightCICPt);
//   ElectronFakeRatePt.AddNumerator(LHIdLoosePt);
//   ElectronFakeRatePt.AddNumerator(LHIdTightPt);
  ElectronFakeRatePt.SetDenominator(RecoPt);
  ElectronFakeRatePt.ComputeEfficiencies();
  ElectronFakeRatePt.SetTitle("jet fake probability vs p_{T}");
  ElectronFakeRatePt.SetXaxisTitle("p_{T} of closest jet [GeV]");
  ElectronFakeRatePt.SetYaxisTitle("jet #rightarrow fake e probability");
  ElectronFakeRatePt.SetYaxisMin(0.0);
  ElectronFakeRatePt.Write();

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
//   TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
//   TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
//   TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );

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
//   TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
//   TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );


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
    for ( int jet=0; jet<nAK5PFJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 25.0 && p3Jet.Pt() > maxEt) {
        maxEt = p3Jet.Pt();
        leadingJet = jet;
      }
    }

    // consider the other as fakes (denominator)
    for ( int jet=0; jet<nAK5PFJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);
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

        for ( int jet=0; jet<nAK5PFJet; jet++ ) {
          TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);
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

          TVector3 p3ClosestJet(pxAK5PFJet[closestJet],pyAK5PFJet[closestJet],pzAK5PFJet[closestJet]);

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
      // Golden means 0 brem clusters, showering otherwise
//       if ( nbremsEle[ele] == 0 ) {
//         GoldenEta->Fill(etaFake);
//         GoldenPt->Fill(etFake);
//       } else {
//         ShoweringEta->Fill(etaFake);
//         ShoweringPt->Fill(etFake);
//       }

      for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
	{
	  bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
	  isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
	  isEleID(&EgammaCutBasedID[icut],ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

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
//   ElectronEffEta.AddNumerator(GoldenEta);
//   ElectronEffEta.AddNumerator(ShoweringEta);
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
//   ElectronEffEtaHighPt.AddNumerator(GoldenEtaHighPt);
//   ElectronEffEtaHighPt.AddNumerator(ShoweringEtaHighPt);
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
//   ElectronEffEtaLowPt.AddNumerator(GoldenEtaLowPt);
//   ElectronEffEtaLowPt.AddNumerator(ShoweringEtaLowPt);
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
//   ElectronEffPt.AddNumerator(GoldenPt);
//   ElectronEffPt.AddNumerator(ShoweringPt);
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
//   ElectronEffPtBarrel.AddNumerator(GoldenPtBarrel);
//   ElectronEffPtBarrel.AddNumerator(ShoweringPtBarrel);
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
//   ElectronEffPtEndcap.AddNumerator(GoldenPtEndcap);
//   ElectronEffPtEndcap.AddNumerator(ShoweringPtEndcap);
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


void LikelihoodAnalysis::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);

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

  float lh=likelihoodRatio(eleIndex,*LH);
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
  selector->SetLikelihood( lh );
  selector->SetEcalIsolation( dr03EcalRecHitSumEtEle[eleIndex]/pEle.Pt() );
  selector->SetTrkIsolation( dr03TkSumPtEle[eleIndex]/pEle.Pt() );
  selector->SetHcalIsolation( dr03HcalTowerSumEtEle[eleIndex]/pEle.Pt() );
  if (isEleEB)
    selector->SetCombinedIsolation( (dr03TkSumPtEle[eleIndex] + 
				     TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + 
				     dr03HcalTowerSumEtEle[eleIndex]) / pEle.Pt() );
  else
    selector->SetCombinedIsolation( (dr03TkSumPtEle[eleIndex] + 
				     dr03EcalRecHitSumEtEle[eleIndex] + 
				     dr03HcalTowerSumEtEle[eleIndex]) / pEle.Pt() );
  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
//   selector->m_cleaner->SetE1(e1);
//   selector->m_cleaner->SetE4SwissCross(e4SwissCross);
//   selector->m_cleaner->SetFiducialFlag(fidFlagSC);
//   selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
//   selector->m_cleaner->SetSeedTime(seedTime);
//   selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();

}

void LikelihoodAnalysis::isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);

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
  //  std::cout << eop << "," <<  eSeedOverPoutEle[eleIndex] << ","<< eseedopin*sce/rawe << "," << fbremEle[eleIndex] << "," << nbremsEle[eleIndex] << std::endl;;
  
  //  float lh=likelihoodRatio(eleIndex,*LH);
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
  selector->SetEcalIsolation( dr04EcalRecHitSumEtEle[eleIndex] );
  selector->SetTrkIsolation( dr03TkSumPtEle[eleIndex] );
  selector->SetHcalIsolation( dr04HcalTowerSumEtEle[eleIndex] );
  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
//   selector->m_cleaner->SetE1(e1);
//   selector->m_cleaner->SetE4SwissCross(e4SwissCross);
//   selector->m_cleaner->SetFiducialFlag(fidFlagSC);
//   selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
//   selector->m_cleaner->SetSeedTime(seedTime);
//   selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputEleId();
  *isolOutput = selector->outputIso();
  *convRejOutput = selector->outputConv();

}

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

float LikelihoodAnalysis::likelihoodRatio(int eleIndex, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  
  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  measurements.pt = pEle.Pt();
  measurements.subdet = (fabs(etaEle[eleIndex])<1.479) ? 0 : 1;
  measurements.deltaPhi = deltaPhiAtVtxEle[eleIndex];
  measurements.deltaEta = deltaEtaAtVtxEle[eleIndex];
  measurements.eSeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  measurements.eSuperClusterOverP = eSuperClusterOverPEle[eleIndex];
  measurements.hadronicOverEm = hOverEEle[eleIndex];
  measurements.sigmaIEtaIEta = SigmaiEiE(eleIndex);
  measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
  measurements.fBrem = fbremEle[eleIndex];
  measurements.nBremClusters = nbremsEle[eleIndex];
  return lh.result(measurements);
}

