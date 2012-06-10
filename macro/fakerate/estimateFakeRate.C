#define estimateFakeRate_cxx
#include "estimateFakeRate.h"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "EgammaAnalysisTools/include/HZZEleIDSelector.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <iostream>

using namespace std;

void estimateFakeRate::Loop(const char *outname)
{
//   In a ROOT session, you can do:
//      Root > .L estimateFakeRate.C
//      Root > estimateFakeRate t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  // friend tree with isolation
  fChain->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", _isofriend);
  Float_t combPFIsoHZZ, combDetIsoHZZ;
  fChain->SetBranchAddress("combPFIsoHZZ",&combPFIsoHZZ);
  fChain->SetBranchAddress("combDetIsoHZZ",&combDetIsoHZZ);

  // friend tree with id
  fChain->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", _idbitsfriend);
  Int_t hzzMvaLoose, hzzMvaTight, newhwwWP;
  fChain->SetBranchAddress("hzzMvaLoose", &hzzMvaLoose);
  fChain->SetBranchAddress("hzzMvaTight", &hzzMvaTight);
  fChain->SetBranchAddress("newhwwWP", &newhwwWP);

  std::vector<TString> EgammaTrgEleID, EgammaNoTrgEleID;
  // this is for triggering electrons
  EgammaTrgEleID.push_back("WP80");   // 2011 def
  EgammaTrgEleID.push_back("WP80EA"); // 2011 def with EAs
  EgammaTrgEleID.push_back("IsoWP80");
  EgammaTrgEleID.push_back("IsoWP80EA");
  EgammaTrgEleID.push_back("IsoWP80ZZNoEA");
  EgammaTrgEleID.push_back("IdWP80");
  EgammaTrgEleID.push_back("newWP80"); // 2012 MVA + EAs iso
  EgammaTrgEleID.push_back("newWP70"); // 2012 MVA + EAs iso
  EgammaTrgEleID.push_back("newWPHWW"); // same eff as 2011 HWW WP

  // no need to repeat the iso only: non triggering eles
  EgammaNoTrgEleID.push_back("WP80");
  EgammaNoTrgEleID.push_back("WP95"); 
  EgammaNoTrgEleID.push_back("CiCMedium"); // using id and iso cuts from CiC[3]
  EgammaNoTrgEleID.push_back("CiCMediumIso1"); 
  EgammaNoTrgEleID.push_back("CiCMediumIso2"); 
  EgammaNoTrgEleID.push_back("newWPHZZ");  // same fake rate as 2011 HZZ WP (isolation optimized too)
  EgammaNoTrgEleID.push_back("newWPHZZIso1");  // same fake rate as 2011 HZZ WP for ID, iso bins
  EgammaNoTrgEleID.push_back("newWPHZZIso2");  // same fake rate as 2011 HZZ WP for ID, iso bins
  EgammaNoTrgEleID.push_back("hzzMvaLoose");  // duncan's optimization for HZZ
  EgammaNoTrgEleID.push_back("hzzMvaTight");  // duncan's optimization for HZZ

  // -----------------------------------------------------------------------
  // study vs eta

  Float_t LowerEta[5];
  LowerEta[0]=0.0;
  LowerEta[1]=1.0;
  LowerEta[2]=1.479;  
  LowerEta[3]=2.0;
  LowerEta[4]=2.5;
  TH1F *RecoEta         = new TH1F( "RecoEta",          "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaHighPt   = new TH1F( "RecoEtaHighPt",    "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaLowPt    = new TH1F( "RecoEtaLowPt",     "reconstructed #eta", 4, LowerEta);
  
  std::vector<TH1F*> TrgEleEta;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"Eta",   "HWW BDT ID #eta", 4, LowerEta);
    TrgEleEta.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgEleEta;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"Eta",   "HZZ BDT ID #eta", 4, LowerEta); 
    NoTrgEleEta.push_back(aHisto);
  }

  // eta, high pT
  std::vector<TH1F*> TrgEleEtaHighPt;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"EtaHighPt",   "HWW BDT ID #eta", 4, LowerEta);   
    TrgEleEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgEleEtaHighPt;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"EtaHighPt",   "HZZ BDT ID #eta", 4, LowerEta); 
    NoTrgEleEtaHighPt.push_back(aHisto);
  }

  // eta, low pT
  std::vector<TH1F*> TrgEleEtaLowPt;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"EtaLowPt",   "HWW BDT ID #eta", 4, LowerEta);  
    TrgEleEtaLowPt.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgEleEtaLowPt;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"EtaLowPt",   "HZZ BDT ID #eta", 4, LowerEta);
    NoTrgEleEtaLowPt.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs pT
  Float_t LowerPt[11] = {0.0,5.0,7.0,10.0,14.0,20.0,25.0,30.0,40.0,50.0,80.0};
  TH1F *RecoPt         = new TH1F( "RecoPt",          "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtBarrel1  = new TH1F( "RecoPtBarrel1",   "reconstructed p_{T} (GeV)", 10, LowerPt);    // this is done to subsplit 
  TH1F *RecoPtBarrel2  = new TH1F( "RecoPtBarrel2",   "reconstructed p_{T} (GeV)", 10, LowerPt);    // barrel and endcap in 2 parts 
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtEndcap1  = new TH1F( "RecoPtEndcap1",   "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtEndcap2  = new TH1F( "RecoPtEndcap2",   "reconstructed p_{T} (GeV)", 10, LowerPt);

  std::vector<TH1F*> TrgElePt;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"Pt",   "HWW BDT ID #eta", 10, LowerPt );
    TrgElePt.push_back(aHisto);
  }
  
  std::vector<TH1F*> NoTrgElePt;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"Pt", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElePt.push_back(aHisto);
  }


  // to have the full picture in the barrel
  std::vector<TH1F*> TrgElePtBarrel;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PtBarrel", "HWW BDT ID #eta",   10, LowerPt );
    TrgElePtBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgElePtBarrel;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PtBarrel", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElePtBarrel.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins ( barrel )
  std::vector<TH1F*> TrgElePtBarrel1;
  std::vector<TH1F*> TrgElePtBarrel2;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PtBarrel1", "HWW BDT ID #eta",   10, LowerPt );
    TrgElePtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PtBarrel2", "HWW BDT ID #eta",   10, LowerPt );
    TrgElePtBarrel2.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgElePtBarrel1;
  std::vector<TH1F*> NoTrgElePtBarrel2;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PtBarrel1", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElePtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PtBarrel2", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElePtBarrel2.push_back(aHisto);
  }

  // to have the full picture in the endcap
  std::vector<TH1F*> TrgElePtEndcap;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PtEndcap", "HWW BDT ID #eta",   10, LowerPt );
    TrgElePtEndcap.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgElePtEndcap;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PtEndcap", "HZZ BDT ID #eta",   10, LowerPt);
    NoTrgElePtEndcap.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins (endcap)
  std::vector<TH1F*> TrgElePtEndcap1;
  std::vector<TH1F*> TrgElePtEndcap2;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PtEndcap1", "HWW BDT ID #eta",   10, LowerPt );
    TrgElePtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PtEndcap2", "HWW BDT ID #eta",   10, LowerPt );
    TrgElePtEndcap2.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgElePtEndcap1;
  std::vector<TH1F*> NoTrgElePtEndcap2;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PtEndcap1", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElePtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PtEndcap2", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElePtEndcap2.push_back(aHisto);
  }


  // -----------------------------------------------------------------------
  // study vs PU
  Float_t LowerPU[11];
  // 2011 binning
  // LowerPU[0] = 1;
  // LowerPU[1] = 4;
  // for(int i=2;i<6;++i) LowerPU[i]=i+3;
  // LowerPU[6] = 10;
  // LowerPU[7] = 15;
  // LowerPU[8] = 25;
  // LowerPU[9] = 30;
  // LowerPU[10] = 35;

  LowerPU[0] = 1;
  LowerPU[1] = 6;
  LowerPU[2] = 10;
  LowerPU[3] = 13;
  LowerPU[4] = 15;
  LowerPU[5] = 17;
  LowerPU[6] = 19;
  LowerPU[7] = 21;
  LowerPU[8] = 26;
  LowerPU[9] = 36;
  LowerPU[10] = 50;

  TH1F *RecoPU         = new TH1F( "RecoPU",         "reconstructed nPU", 10, LowerPU );
  TH1F *RecoPUBarrel   = new TH1F( "RecoPUBarrel",   "reconstructed nPU", 10, LowerPU);    // this is done to subsplit 
  TH1F *RecoPUBarrel1  = new TH1F( "RecoPUBarrel1",  "reconstructed nPU", 10, LowerPU);    // this is done to subsplit 
  TH1F *RecoPUBarrel2  = new TH1F( "RecoPUBarrel2",  "reconstructed nPU", 10, LowerPU);    // barrel and endcap in 2 parts 
  TH1F *RecoPUEndcap   = new TH1F( "RecoPUEndcap",   "reconstructed nPU", 10, LowerPU);
  TH1F *RecoPUEndcap1  = new TH1F( "RecoPUEndcap1",  "reconstructed nPU", 10, LowerPU);
  TH1F *RecoPUEndcap2  = new TH1F( "RecoPUEndcap2",  "reconstructed nPU", 10, LowerPU);

  std::vector<TH1F*> TrgElePU;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PU", "HWW BDT ID nPU",  10, LowerPU );
    TrgElePU.push_back(aHisto);
  }
  
  std::vector<TH1F*> NoTrgElePU;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PU", "HZZ BDT ID nPU",  10, LowerPU );
    NoTrgElePU.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins ( barrel )
  std::vector<TH1F*> TrgElePUBarrel, TrgElePUBarrel1, TrgElePUBarrel2;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PUBarrel1", "HWW BDT ID #eta",   10, LowerPU );
    TrgElePUBarrel1.push_back(aHisto);
    aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PUBarrel2", "HWW BDT ID #eta",   10, LowerPU );
    TrgElePUBarrel2.push_back(aHisto);
    aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PUBarrel", "HWW BDT ID #eta",   10, LowerPU );
    TrgElePUBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgElePUBarrel, NoTrgElePUBarrel1, NoTrgElePUBarrel2;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PUBarrel1", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElePUBarrel1.push_back(aHisto);
    aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PUBarrel2", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElePUBarrel2.push_back(aHisto);
    aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PUBarrel", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElePUBarrel.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins (endcap)
  std::vector<TH1F*> TrgElePUEndcap, TrgElePUEndcap1, TrgElePUEndcap2;
  for (int i=0;i<(int)EgammaTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PUEndcap1", "HWW BDT ID #eta",   10, LowerPU );
    TrgElePUEndcap1.push_back(aHisto);
    aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PUEndcap2", "HWW BDT ID #eta",   10, LowerPU );
    TrgElePUEndcap2.push_back(aHisto);
    aHisto = new TH1F( "TrgEle"+TString(EgammaTrgEleID[i])+"PUEndcap", "HWW BDT ID #eta",   10, LowerPU );
    TrgElePUEndcap.push_back(aHisto);
  }

  std::vector<TH1F*> NoTrgElePUEndcap, NoTrgElePUEndcap1, NoTrgElePUEndcap2;
  for (int i=0;i<(int)EgammaNoTrgEleID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PUEndcap1", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElePUEndcap1.push_back(aHisto);
    aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PUEndcap2", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElePUEndcap2.push_back(aHisto);
    aHisto = new TH1F( "NoTrgEle"+TString(EgammaNoTrgEleID[i])+"PUEndcap", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElePUEndcap.push_back(aHisto);
  }

  // this is the utils class to apply the 2012 WPs
  HZZEleIDSelector aSel;

  // loop on events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // NOT EWK CORRECTED: STOP AT 35 GEV!!!
    if(!TString(outname).Contains("zee1fake") && (pt>35 || !DenomFakeSmurf)) continue;
    // JET PT threshold
    if(leadJetPt<35) continue;

    // fill the denominator: take only the highest pT denominator candidate
    float etaFake = fabs(eta);
    float etFake  = pt;
    bool isInEB   = fabs(eta)<1.479;
    bool isInEE   = !isInEB;
    bool highPt   = (pt>20.);
    bool lowPt    = (pt<=20.);
    // to split in 4 eta regions
    int etaRegion = -1;
    if (fabs(etaFake)>=0. && fabs(etaFake)<1.)    etaRegion = 1; 
    if (fabs(etaFake)>=1. && fabs(etaFake)<1.479) etaRegion = 2;
    if (fabs(etaFake)>=1.479 && fabs(etaFake)<2.) etaRegion = 3; 
    if (fabs(etaFake)>=2. && fabs(etaFake)<2.5)   etaRegion = 4; 

    // filling
    RecoEta         -> Fill(etaFake);   // , theWeight);
    RecoPt          -> Fill(etFake);    // , theWeight);
    RecoPU          -> Fill(vertices);  // , theWeight);

    if (highPt) RecoEtaHighPt -> Fill(etaFake);  //, theWeight); 
    if (lowPt)  RecoEtaLowPt  -> Fill(etaFake);  // , theWeight);

    if (isInEB) {
      RecoPtBarrel -> Fill(etFake);  //, theWeight);
      RecoPUBarrel -> Fill(vertices);  //, theWeight);
    }
    if (isInEE) {
      RecoPtEndcap -> Fill(etFake);  //, theWeight);
      RecoPUEndcap -> Fill(vertices);  //, theWeight);
    }

    if (etaRegion==1) {
      RecoPtBarrel1 -> Fill(etFake);  //, theWeight);      
      RecoPUBarrel1 -> Fill(vertices);  //, theWeight); 
    }
    if (etaRegion==2) {
      RecoPtBarrel2 -> Fill(etFake);  //, theWeight);      
      RecoPUBarrel2 -> Fill(vertices);  //, theWeight); 
    }
    if (etaRegion==3) {
      RecoPtEndcap1 -> Fill(etFake);  //, theWeight);      
      RecoPUEndcap1 -> Fill(vertices);  //, theWeight);  
    }
    if (etaRegion==4) {
      RecoPtEndcap2 -> Fill(etFake);  //, theWeight);      
      RecoPUEndcap2 -> Fill(vertices);  //, theWeight); 
    }

    // fill the numerator(s)
    // === WP80 HWW 2011: BDT + ISO ===
    if(passID(kIsoHWW2011) && passID(kBDTHWW2011_withIP)) {
      TrgEleEta[kWP80]->Fill(etaFake);
      TrgElePt[kWP80] ->Fill(etFake);
      TrgElePU[kWP80] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[kWP80]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[kWP80] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[kWP80] ->Fill(etFake);
	TrgElePUBarrel[kWP80] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[kWP80] ->Fill(etFake);
	TrgElePUEndcap[kWP80] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[kWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[kWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[kWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[kWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[kWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[kWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[kWP80] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[kWP80] -> Fill(vertices); //, theWeight);     
      }
    }
    // === WP80 HWW BDT + EA corr Iso ===
    if(passID(kIsoEACorr) && passID(kBDTHWW2011_withIP)) {
      TrgEleEta[kWP80EA]->Fill(etaFake);
      TrgElePt[kWP80EA] ->Fill(etFake);
      TrgElePU[kWP80EA] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[kWP80EA]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[kWP80EA] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[kWP80EA] ->Fill(etFake);
	TrgElePUBarrel[kWP80EA] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[kWP80EA] ->Fill(etFake);
	TrgElePUEndcap[kWP80EA] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[kWP80EA] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[kWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[kWP80EA] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[kWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[kWP80EA] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[kWP80EA] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[kWP80EA] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[kWP80EA] -> Fill(vertices); //, theWeight);     
      }
    }
    // === HWW 2011 Iso ===
    if(passID(kIsoHWW2011)) {
      TrgEleEta[kIsoWP80]->Fill(etaFake);
      TrgElePt[kIsoWP80] ->Fill(etFake);
      TrgElePU[kIsoWP80] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[kIsoWP80]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[kIsoWP80] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[kIsoWP80] ->Fill(etFake);
	TrgElePUBarrel[kIsoWP80] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[kIsoWP80] ->Fill(etFake);
	TrgElePUEndcap[kIsoWP80] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[kIsoWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[kIsoWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[kIsoWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[kIsoWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[kIsoWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[kIsoWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[kIsoWP80] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[kIsoWP80] -> Fill(vertices); //, theWeight);     
      }
    }
    // === EA corr Iso ===
    if(passID(kIsoEACorr)) {
      TrgEleEta[kIsoWP80EA]->Fill(etaFake);
      TrgElePt[kIsoWP80EA] ->Fill(etFake);
      TrgElePU[kIsoWP80EA] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[kIsoWP80EA]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[kIsoWP80EA] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[kIsoWP80EA] ->Fill(etFake);
	TrgElePUBarrel[kIsoWP80EA] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[kIsoWP80EA] ->Fill(etFake);
	TrgElePUEndcap[kIsoWP80EA] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[kIsoWP80EA] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[kIsoWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[kIsoWP80EA] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[kIsoWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[kIsoWP80EA] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[kIsoWP80EA] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[kIsoWP80EA] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[kIsoWP80EA] -> Fill(vertices); //, theWeight);     
      }
    }
    // === HZZ, but _not EA corr_ Iso ===
    if(passID(kIso)) {
      TrgEleEta[kIsoWP80HZZNoEA]->Fill(etaFake);
      TrgElePt[kIsoWP80HZZNoEA] ->Fill(etFake);
      TrgElePU[kIsoWP80HZZNoEA] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[kIsoWP80HZZNoEA]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[kIsoWP80HZZNoEA] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[kIsoWP80HZZNoEA] ->Fill(etFake);
	TrgElePUBarrel[kIsoWP80HZZNoEA] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[kIsoWP80HZZNoEA] ->Fill(etFake);
	TrgElePUEndcap[kIsoWP80HZZNoEA] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[kIsoWP80HZZNoEA] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[kIsoWP80HZZNoEA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[kIsoWP80HZZNoEA] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[kIsoWP80HZZNoEA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[kIsoWP80HZZNoEA] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[kIsoWP80HZZNoEA] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[kIsoWP80HZZNoEA] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[kIsoWP80HZZNoEA] -> Fill(vertices); //, theWeight);     
      }
    }

    // === HWW 2011 BDT ===
    if(passID(kBDTHWW2011_withIP)) {
      TrgEleEta[kIdWP80]->Fill(etaFake);
      TrgElePt[kIdWP80] ->Fill(etFake);
      TrgElePU[kIdWP80] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[kIdWP80]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[kIdWP80] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[kIdWP80] ->Fill(etFake);
	TrgElePUBarrel[kIdWP80] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[kIdWP80] ->Fill(etFake);
	TrgElePUEndcap[kIdWP80] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[kIdWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[kIdWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[kIdWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[kIdWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[kIdWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[kIdWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[kIdWP80] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[kIdWP80] -> Fill(vertices); //, theWeight);     
      }
    }

    // === HWW 2012 WP80 BDT ===
    if(aSel.output(etFake,etaFake,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVABiased)) {
      TrgEleEta[knewWP80]->Fill(etaFake);
      TrgElePt[knewWP80] ->Fill(etFake);
      TrgElePU[knewWP80] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[knewWP80]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[knewWP80] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[knewWP80] ->Fill(etFake);
	TrgElePUBarrel[knewWP80] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[knewWP80] ->Fill(etFake);
	TrgElePUEndcap[knewWP80] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[knewWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[knewWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[knewWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[knewWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[knewWP80] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[knewWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[knewWP80] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[knewWP80] -> Fill(vertices); //, theWeight);     
      }
    }

    // === HWW 2012 WP70 BDT ===
    if(aSel.output(etFake,etaFake,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVABiased)) {
      TrgEleEta[knewWP70]->Fill(etaFake);
      TrgElePt[knewWP70] ->Fill(etFake);
      TrgElePU[knewWP70] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[knewWP70]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[knewWP70] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[knewWP70] ->Fill(etFake);
	TrgElePUBarrel[knewWP70] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[knewWP70] ->Fill(etFake);
	TrgElePUEndcap[knewWP70] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[knewWP70] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[knewWP70] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[knewWP70] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[knewWP70] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[knewWP70] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[knewWP70] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[knewWP70] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[knewWP70] -> Fill(vertices); //, theWeight);     
      }
    }

    // === HWW 2012 ~WP70x80 BDT ===
    if(newhwwWP) {
      TrgEleEta[knewWPHWW]->Fill(etaFake);
      TrgElePt[knewWPHWW] ->Fill(etFake);
      TrgElePU[knewWPHWW] ->Fill(vertices);
      if (highPt) TrgEleEtaHighPt[knewWPHWW]->Fill(etaFake);
      if (lowPt)  TrgEleEtaLowPt[knewWPHWW] ->Fill(etaFake);
      if (isInEB) { 
	TrgElePtBarrel[knewWPHWW] ->Fill(etFake);
	TrgElePUBarrel[knewWPHWW] ->Fill(vertices);
      }
      if (isInEE) {
	TrgElePtEndcap[knewWPHWW] ->Fill(etFake);
	TrgElePUEndcap[knewWPHWW] ->Fill(vertices);
      }
      if (etaRegion==1) {
	TrgElePtBarrel1[knewWPHWW] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel1[knewWPHWW] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	TrgElePtBarrel2[knewWPHWW] -> Fill(etFake);   //, theWeight);      
	TrgElePUBarrel2[knewWPHWW] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	TrgElePtEndcap1[knewWPHWW] -> Fill(etFake);   //, theWeight);      
	TrgElePUEndcap1[knewWPHWW] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	TrgElePtEndcap2[knewWPHWW] -> Fill(etFake);   //, theWeight);     
	TrgElePUEndcap2[knewWPHWW] -> Fill(vertices); //, theWeight);     
      }
    }

    // === WP80 HZZ reoptimized cut ===
    if(aSel.output(etFake,etaFake,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVAUnbiased)) {
      NoTrgEleEta[kHzzWP80]->Fill(etaFake);
      NoTrgElePt[kHzzWP80] ->Fill(etFake);
      NoTrgElePU[kHzzWP80] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[kHzzWP80]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[kHzzWP80] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[kHzzWP80] ->Fill(etFake);
	NoTrgElePUBarrel[kHzzWP80] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[kHzzWP80] ->Fill(etFake);
	NoTrgElePUEndcap[kHzzWP80] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[kHzzWP80] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[kHzzWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[kHzzWP80] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[kHzzWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[kHzzWP80] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[kHzzWP80] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[kHzzWP80] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[kHzzWP80] -> Fill(vertices); //, theWeight);
      }

    }

    // === WP95 HZZ reoptimized cut ===
    if(aSel.output(etFake,etaFake,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVAUnbiased)) {
      NoTrgEleEta[kHzzWP95]->Fill(etaFake);
      NoTrgElePt[kHzzWP95] ->Fill(etFake);
      NoTrgElePU[kHzzWP95] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[kHzzWP95]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[kHzzWP95] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[kHzzWP95] ->Fill(etFake);
	NoTrgElePUBarrel[kHzzWP95] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[kHzzWP95] ->Fill(etFake);
	NoTrgElePUEndcap[kHzzWP95] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[kHzzWP95] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[kHzzWP95] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[kHzzWP95] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[kHzzWP95] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[kHzzWP95] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[kHzzWP95] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[kHzzWP95] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[kHzzWP95] -> Fill(vertices); //, theWeight);
      }
    }

    // === CiC used in HZZ in 2011: iso from CiC[3] ===
    if(missHits<=1 && cicid(3) && ciciso(3)) {
      NoTrgEleEta[kHzzCiCMedium]->Fill(etaFake);
      NoTrgElePt[kHzzCiCMedium] ->Fill(etFake);
      NoTrgElePU[kHzzCiCMedium] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[kHzzCiCMedium]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[kHzzCiCMedium] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[kHzzCiCMedium] ->Fill(etFake);
	NoTrgElePUBarrel[kHzzCiCMedium] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[kHzzCiCMedium] ->Fill(etFake);
	NoTrgElePUEndcap[kHzzCiCMedium] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[kHzzCiCMedium] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[kHzzCiCMedium] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[kHzzCiCMedium] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[kHzzCiCMedium] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[kHzzCiCMedium] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[kHzzCiCMedium] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[kHzzCiCMedium] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[kHzzCiCMedium] -> Fill(vertices); //, theWeight);
      }
    }

    // === CiC used in HZZ in 2011: iso bin 1 ===
    if(missHits<=1 && cicid(3) && combDetIsoHZZ/pt<0.175) {
      NoTrgEleEta[kHzzCiCMediumIso1]->Fill(etaFake);
      NoTrgElePt[kHzzCiCMediumIso1] ->Fill(etFake);
      NoTrgElePU[kHzzCiCMediumIso1] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[kHzzCiCMediumIso1]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[kHzzCiCMediumIso1] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[kHzzCiCMediumIso1] ->Fill(etFake);
	NoTrgElePUBarrel[kHzzCiCMediumIso1] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[kHzzCiCMediumIso1] ->Fill(etFake);
	NoTrgElePUEndcap[kHzzCiCMediumIso1] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[kHzzCiCMediumIso1] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[kHzzCiCMediumIso1] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[kHzzCiCMediumIso1] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[kHzzCiCMediumIso1] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[kHzzCiCMediumIso1] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[kHzzCiCMediumIso1] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[kHzzCiCMediumIso1] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[kHzzCiCMediumIso1] -> Fill(vertices); //, theWeight);
      }
    }

    // === CiC used in HZZ in 2011: iso bin 2 ===
    if(missHits<=1 && cicid(3) && combDetIsoHZZ/pt>0.175 && combDetIsoHZZ/pt<0.35) {
      NoTrgEleEta[kHzzCiCMediumIso2]->Fill(etaFake);
      NoTrgElePt[kHzzCiCMediumIso2] ->Fill(etFake);
      NoTrgElePU[kHzzCiCMediumIso2] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[kHzzCiCMediumIso2]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[kHzzCiCMediumIso2] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[kHzzCiCMediumIso2] ->Fill(etFake);
	NoTrgElePUBarrel[kHzzCiCMediumIso2] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[kHzzCiCMediumIso2] ->Fill(etFake);
	NoTrgElePUEndcap[kHzzCiCMediumIso2] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[kHzzCiCMediumIso2] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[kHzzCiCMediumIso2] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[kHzzCiCMediumIso2] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[kHzzCiCMediumIso2] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[kHzzCiCMediumIso2] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[kHzzCiCMediumIso2] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[kHzzCiCMediumIso2] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[kHzzCiCMediumIso2] -> Fill(vertices); //, theWeight);
      }
    }
    
    // === WP with same fake rate as CiC used in HZZ in 2011 (cross-check of tuning) ===
    if(missHits<=1 && aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWPHZZ,HZZEleIDSelector::kMVAUnbiased)) {
      NoTrgEleEta[knewWPHZZ]->Fill(etaFake);
      NoTrgElePt[knewWPHZZ] ->Fill(etFake);
      NoTrgElePU[knewWPHZZ] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[knewWPHZZ]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[knewWPHZZ] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[knewWPHZZ] ->Fill(etFake);
	NoTrgElePUBarrel[knewWPHZZ] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[knewWPHZZ] ->Fill(etFake);
	NoTrgElePUEndcap[knewWPHZZ] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[knewWPHZZ] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[knewWPHZZ] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[knewWPHZZ] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[knewWPHZZ] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[knewWPHZZ] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[knewWPHZZ] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[knewWPHZZ] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[knewWPHZZ] -> Fill(vertices); //, theWeight);
      }
    }

    // === MVA for ID, iso cut: iso bin 1 ===
    if(missHits<=1 && aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWPHZZ,HZZEleIDSelector::kMVAUnbiased) && combDetIsoHZZ/pt<0.175) {
      NoTrgEleEta[knewWPHZZIso1]->Fill(etaFake);
      NoTrgElePt[knewWPHZZIso1] ->Fill(etFake);
      NoTrgElePU[knewWPHZZIso1] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[knewWPHZZIso1]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[knewWPHZZIso1] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[knewWPHZZIso1] ->Fill(etFake);
	NoTrgElePUBarrel[knewWPHZZIso1] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[knewWPHZZIso1] ->Fill(etFake);
	NoTrgElePUEndcap[knewWPHZZIso1] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[knewWPHZZIso1] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[knewWPHZZIso1] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[knewWPHZZIso1] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[knewWPHZZIso1] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[knewWPHZZIso1] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[knewWPHZZIso1] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[knewWPHZZIso1] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[knewWPHZZIso1] -> Fill(vertices); //, theWeight);
      }
    }

    // === MVA for ID, iso cut: iso bin 2 ===
    if(missHits<=1 && aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWPHZZ,HZZEleIDSelector::kMVAUnbiased) && combDetIsoHZZ/pt>0.175 && combDetIsoHZZ/pt<0.35) {
      NoTrgEleEta[knewWPHZZIso2]->Fill(etaFake);
      NoTrgElePt[knewWPHZZIso2] ->Fill(etFake);
      NoTrgElePU[knewWPHZZIso2] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[knewWPHZZIso2]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[knewWPHZZIso2] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[knewWPHZZIso2] ->Fill(etFake);
	NoTrgElePUBarrel[knewWPHZZIso2] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[knewWPHZZIso2] ->Fill(etFake);
	NoTrgElePUEndcap[knewWPHZZIso2] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[knewWPHZZIso2] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[knewWPHZZIso2] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[knewWPHZZIso2] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[knewWPHZZIso2] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[knewWPHZZIso2] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[knewWPHZZIso2] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[knewWPHZZIso2] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[knewWPHZZIso2] -> Fill(vertices); //, theWeight);
      }
    }

    // === MVA for ID, MVA for Iso ===
    if(missHits<=1 && hzzMvaLoose) {
      NoTrgEleEta[khzzMvaLoose]->Fill(etaFake);
      NoTrgElePt[khzzMvaLoose] ->Fill(etFake);
      NoTrgElePU[khzzMvaLoose] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[khzzMvaLoose]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[khzzMvaLoose] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[khzzMvaLoose] ->Fill(etFake);
	NoTrgElePUBarrel[khzzMvaLoose] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[khzzMvaLoose] ->Fill(etFake);
	NoTrgElePUEndcap[khzzMvaLoose] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[khzzMvaLoose] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[khzzMvaLoose] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[khzzMvaLoose] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[khzzMvaLoose] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[khzzMvaLoose] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[khzzMvaLoose] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[khzzMvaLoose] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[khzzMvaLoose] -> Fill(vertices); //, theWeight);
      }
    }

    // === MVA for ID, MVA for Iso ===
    if(missHits<=1 && hzzMvaTight) {
      NoTrgEleEta[khzzMvaTight]->Fill(etaFake);
      NoTrgElePt[khzzMvaTight] ->Fill(etFake);
      NoTrgElePU[khzzMvaTight] ->Fill(vertices);
      if (highPt) NoTrgEleEtaHighPt[khzzMvaTight]->Fill(etaFake);
      if (lowPt)  NoTrgEleEtaLowPt[khzzMvaTight] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElePtBarrel[khzzMvaTight] ->Fill(etFake);
	NoTrgElePUBarrel[khzzMvaTight] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgElePtEndcap[khzzMvaTight] ->Fill(etFake);
	NoTrgElePUEndcap[khzzMvaTight] ->Fill(vertices);
      }
      if (etaRegion==1) {
	NoTrgElePtBarrel1[khzzMvaTight] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel1[khzzMvaTight] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	NoTrgElePtBarrel2[khzzMvaTight] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUBarrel2[khzzMvaTight] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	NoTrgElePtEndcap1[khzzMvaTight] -> Fill(etFake);   //, theWeight);      
	NoTrgElePUEndcap1[khzzMvaTight] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	NoTrgElePtEndcap2[khzzMvaTight] -> Fill(etFake);   //, theWeight);     
	NoTrgElePUEndcap2[khzzMvaTight] -> Fill(vertices); //, theWeight);
      }
    }

  } // loop over events


  // saving efficiency histos
  // === as a function of eta ===
  char filename[500];
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  ElectronEffEta.AddNumerator(RecoEta);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffEta.AddNumerator(TrgEleEta[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffEta.AddNumerator(NoTrgEleEta[icut]);
  }

  ElectronEffEta.SetDenominator(RecoEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("fake rate vs eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("Fake rate");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleMisidEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(TrgEleEtaHighPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(NoTrgEleEtaHighPt[icut]);
  }

  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs pT");
  ElectronEffEtaHighPt.SetXaxisTitle("electron p_{T}");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(TrgEleEtaLowPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(NoTrgEleEtaLowPt[icut]);
  }

  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("fake rate vs eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  // === as a function of pt ===
  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  ElectronEffPt.AddNumerator(RecoPt);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPt.AddNumerator(TrgElePt[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPt.AddNumerator(NoTrgElePt[icut]);
  }

  ElectronEffPt.SetDenominator(RecoPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("fake rate vs pT");
  ElectronEffPt.SetXaxisTitle("electron pT");
  ElectronEffPt.SetYaxisTitle("Fake rate");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleMisidPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPtBarrel.AddNumerator(TrgElePtBarrel[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPtBarrel.AddNumerator(NoTrgElePtBarrel[icut]);
  }

  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-EleMisidPtBarrel1.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel1(filename);
  ElectronEffPtBarrel1.AddNumerator(RecoPtBarrel1);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPtBarrel1.AddNumerator(TrgElePtBarrel1[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPtBarrel1.AddNumerator(NoTrgElePtBarrel1[icut]);
  }

  ElectronEffPtBarrel1.SetDenominator(RecoPtBarrel1);
  ElectronEffPtBarrel1.ComputeEfficiencies();
  ElectronEffPtBarrel1.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel1.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel1.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel1.SetYaxisMin(0.0);
  ElectronEffPtBarrel1.Write();

  sprintf(filename,"%s-EleMisidPtBarrel2.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel2(filename);
  ElectronEffPtBarrel2.AddNumerator(RecoPtBarrel2);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPtBarrel2.AddNumerator(TrgElePtBarrel2[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPtBarrel2.AddNumerator(NoTrgElePtBarrel2[icut]);
  }

  ElectronEffPtBarrel2.SetDenominator(RecoPtBarrel2);
  ElectronEffPtBarrel2.ComputeEfficiencies();
  ElectronEffPtBarrel2.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel2.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel2.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel2.SetYaxisMin(0.0);
  ElectronEffPtBarrel2.Write();

  sprintf(filename,"%s-EleMisidPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPtEndcap.AddNumerator(TrgElePtEndcap[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPtEndcap.AddNumerator(NoTrgElePtEndcap[icut]);
  }

  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("fake rate vs pT");
  ElectronEffPtEndcap.SetXaxisTitle("electron pT");
  ElectronEffPtEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();

  sprintf(filename,"%s-EleMisidPtEndcap1.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap1(filename);
  ElectronEffPtEndcap1.AddNumerator(RecoPtEndcap1);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPtEndcap1.AddNumerator(TrgElePtEndcap1[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPtEndcap1.AddNumerator(NoTrgElePtEndcap1[icut]);
  }

  ElectronEffPtEndcap1.SetDenominator(RecoPtEndcap1);
  ElectronEffPtEndcap1.ComputeEfficiencies();
  ElectronEffPtEndcap1.SetTitle("fake rate vs pT");
  ElectronEffPtEndcap1.SetXaxisTitle("electron pT");
  ElectronEffPtEndcap1.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap1.SetYaxisMin(0.0);
  ElectronEffPtEndcap1.Write();

  sprintf(filename,"%s-EleMisidPtEndcap2.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap2(filename);
  ElectronEffPtEndcap2.AddNumerator(RecoPtEndcap2);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPtEndcap2.AddNumerator(TrgElePtEndcap2[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPtEndcap2.AddNumerator(NoTrgElePtEndcap2[icut]);
  }

  ElectronEffPtEndcap2.SetDenominator(RecoPtEndcap2);
  ElectronEffPtEndcap2.ComputeEfficiencies();
  ElectronEffPtEndcap2.SetTitle("fake rate vs pT");
  ElectronEffPtEndcap2.SetXaxisTitle("electron pT");
  ElectronEffPtEndcap2.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap2.SetYaxisMin(0.0);
  ElectronEffPtEndcap2.Write();

  // === as a function of PU ===
  sprintf(filename,"%s-EleMisidPU.root",outname);
  EfficiencyEvaluator ElectronEffPU(filename);
  ElectronEffPU.AddNumerator(RecoPU);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPU.AddNumerator(TrgElePU[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPU.AddNumerator(NoTrgElePU[icut]);
  }

  ElectronEffPU.SetDenominator(RecoPU);
  ElectronEffPU.ComputeEfficiencies();
  ElectronEffPU.SetTitle("fake rate vs PU");
  ElectronEffPU.SetXaxisTitle("# vertices");
  ElectronEffPU.SetYaxisTitle("Fake rate");
  ElectronEffPU.SetYaxisMin(0.0);
  ElectronEffPU.Write();

  sprintf(filename,"%s-EleMisidPUBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPUBarrel(filename);
  ElectronEffPUBarrel.AddNumerator(RecoPUBarrel);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPUBarrel.AddNumerator(TrgElePUBarrel[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPUBarrel.AddNumerator(NoTrgElePUBarrel[icut]);
  }

  ElectronEffPUBarrel.SetDenominator(RecoPUBarrel);
  ElectronEffPUBarrel.ComputeEfficiencies();
  ElectronEffPUBarrel.SetTitle("fake rate vs vertices");
  ElectronEffPUBarrel.SetXaxisTitle("# vertices");
  ElectronEffPUBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPUBarrel.SetYaxisMin(0.0);
  ElectronEffPUBarrel.Write();

  sprintf(filename,"%s-EleMisidPUBarrel1.root",outname);
  EfficiencyEvaluator ElectronEffPUBarrel1(filename);
  ElectronEffPUBarrel1.AddNumerator(RecoPUBarrel1);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPUBarrel1.AddNumerator(TrgElePUBarrel1[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPUBarrel1.AddNumerator(NoTrgElePUBarrel1[icut]);
  }

  ElectronEffPUBarrel1.SetDenominator(RecoPUBarrel1);
  ElectronEffPUBarrel1.ComputeEfficiencies();
  ElectronEffPUBarrel1.SetTitle("fake rate vs vertices");
  ElectronEffPUBarrel1.SetXaxisTitle("# vertices");
  ElectronEffPUBarrel1.SetYaxisTitle("Fake rate");
  ElectronEffPUBarrel1.SetYaxisMin(0.0);
  ElectronEffPUBarrel1.Write();

  sprintf(filename,"%s-EleMisidPUBarrel2.root",outname);
  EfficiencyEvaluator ElectronEffPUBarrel2(filename);
  ElectronEffPUBarrel2.AddNumerator(RecoPUBarrel2);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPUBarrel2.AddNumerator(TrgElePUBarrel2[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPUBarrel2.AddNumerator(NoTrgElePUBarrel2[icut]);
  }

  ElectronEffPUBarrel2.SetDenominator(RecoPUBarrel2);
  ElectronEffPUBarrel2.ComputeEfficiencies();
  ElectronEffPUBarrel2.SetTitle("fake rate vs vertices");
  ElectronEffPUBarrel2.SetXaxisTitle("# vertices");
  ElectronEffPUBarrel2.SetYaxisTitle("Fake rate");
  ElectronEffPUBarrel2.SetYaxisMin(0.0);
  ElectronEffPUBarrel2.Write();

  sprintf(filename,"%s-EleMisidPUEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPUEndcap(filename);
  ElectronEffPUEndcap.AddNumerator(RecoPUEndcap);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPUEndcap.AddNumerator(TrgElePUEndcap[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPUEndcap.AddNumerator(NoTrgElePUEndcap[icut]);
  }

  ElectronEffPUEndcap.SetDenominator(RecoPUEndcap);
  ElectronEffPUEndcap.ComputeEfficiencies();
  ElectronEffPUEndcap.SetTitle("fake rate vs vertices");
  ElectronEffPUEndcap.SetXaxisTitle("# vertices");
  ElectronEffPUEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPUEndcap.SetYaxisMin(0.0);
  ElectronEffPUEndcap.Write();

  sprintf(filename,"%s-EleMisidPUEndcap1.root",outname);
  EfficiencyEvaluator ElectronEffPUEndcap1(filename);
  ElectronEffPUEndcap1.AddNumerator(RecoPUEndcap1);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPUEndcap1.AddNumerator(TrgElePUEndcap1[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPUEndcap1.AddNumerator(NoTrgElePUEndcap1[icut]);
  }

  ElectronEffPUEndcap1.SetDenominator(RecoPUEndcap1);
  ElectronEffPUEndcap1.ComputeEfficiencies();
  ElectronEffPUEndcap1.SetTitle("fake rate vs vertices");
  ElectronEffPUEndcap1.SetXaxisTitle("# vertices");
  ElectronEffPUEndcap1.SetYaxisTitle("Fake rate");
  ElectronEffPUEndcap1.SetYaxisMin(0.0);
  ElectronEffPUEndcap1.Write();

  sprintf(filename,"%s-EleMisidPUEndcap2.root",outname);
  EfficiencyEvaluator ElectronEffPUEndcap2(filename);
  ElectronEffPUEndcap2.AddNumerator(RecoPUEndcap2);
  for (int icut=0;icut<(int)EgammaTrgEleID.size();++icut){
    ElectronEffPUEndcap2.AddNumerator(TrgElePUEndcap2[icut]);
  }
  for (int icut=0;icut<(int)EgammaNoTrgEleID.size();++icut){
    ElectronEffPUEndcap2.AddNumerator(NoTrgElePUEndcap2[icut]);
  }

  ElectronEffPUEndcap2.SetDenominator(RecoPUEndcap2);
  ElectronEffPUEndcap2.ComputeEfficiencies();
  ElectronEffPUEndcap2.SetTitle("fake rate vs vertices");
  ElectronEffPUEndcap2.SetXaxisTitle("# vertices");
  ElectronEffPUEndcap2.SetYaxisTitle("Fake rate");
  ElectronEffPUEndcap2.SetYaxisMin(0.0);
  ElectronEffPUEndcap2.Write();

  
}

bool estimateFakeRate::passID(estimateFakeRate::idType type) {
  if(type == kIsoHWW2011) {
    if(fabs(eta)<1.479) return (combPFIsoHWW/pt < 0.13);
    else return (combPFIsoHWW/pt < 0.09);
  }

  // THIS GIVES THE SAME FR AS 2011
  // if(type == kIsoEACorr) {
  //   float combIso=chaPFIso+neuPFIso+phoPFIso;

  //   if(fabs(eta) <  1.0) combIso -= 0.18 * rho;
  //   if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) combIso -= 0.19 * rho;
  //   if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) combIso -= 0.21 * rho;
  //   if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) combIso -= 0.38 * rho;
  //   if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) combIso -= 0.61 * rho;
  //   if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) combIso -= 0.73 * rho;
  //   if(fabs(eta) >=  2.4) combIso -= 0.78 * rho;

  //   if(pt>20) {
  //     if(fabs(eta) <  1.0) return (combIso/pt < 0.22); 
  //     if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.21);
  //     if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.12);
  //     if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.11);
  //     if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.074);
  //     if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.053);
  //     if(fabs(eta) >=  2.4) return (combIso/pt < 0.010);
  //   } else {
  //     if(fabs(eta) <  1.0) return (combIso/pt < 0.20); 
  //     if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.24);
  //     if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.13);
  //     if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.083);
  //     if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.01);
  //     if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.027);
  //     if(fabs(eta) >=  2.4) return (combIso/pt < -0.035);
  //   }
  // }

  // THIS GIVES THE SAME EFFICIENCY AS 2011
  if(type == kIsoEACorr) {
    float combIso=chaPFIso+neuPFIso+phoPFIso;

    if(fabs(eta) <  1.0) combIso -= 0.18 * rho;
    if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) combIso -= 0.19 * rho;
    if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) combIso -= 0.21 * rho;
    if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) combIso -= 0.38 * rho;
    if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) combIso -= 0.61 * rho;
    if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) combIso -= 0.73 * rho;
    if(fabs(eta) >=  2.4) combIso -= 0.78 * rho;

    if(pt>20) {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.14); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.15);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.09);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.058);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.016);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.0071);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.01);
    } else {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.13); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.14);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.079);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.028);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.031);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.059);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.069);
    }
  }

  if(type == kIso) {
    float combIso=chaPFIso+neuPFIso+phoPFIso;
    if(fabs(eta) < 1.479) return (combIso/pt < 0.29);
    else return (combIso/pt < 0.21);
  }

  if(type == kBDTHWW2011_withIP) {
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww[0] > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww[0] > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww[0] > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww[0] > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww[0] > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww[0] > 0.884);
  }

  if(type == kBDTHWW2011_noIP) {
    // warning: no reoptimized WP!
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww[1] > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww[1] > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww[1] > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww[1] > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww[1] > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww[1] > 0.884);
  }

  if(type == kBDTHWW2012_noIP) {
    // same EFFICIENCY as 2011
    if(pt < 20 && fabs(eta) < 1.0) return (newbdthww[3] > -0.578);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (newbdthww[3] > -0.146);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (newbdthww[3] > -0.182);
    if(pt >= 20 && fabs(eta) < 1.0) return (newbdthww[3] > 0.806);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (newbdthww[3] > 0.890);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (newbdthww[3] > 0.814);
  }

  return false;
}

