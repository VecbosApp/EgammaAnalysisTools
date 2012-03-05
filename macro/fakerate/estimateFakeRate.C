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
  fChain->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", "../results_data_fakes/merged_hzzisoFriend.root");
  Float_t combPFIsoHZZ;
  fChain->SetBranchAddress("combPFIsoHZZ",&combPFIsoHZZ);

  std::vector<TString> EgammaBdtHWWBasedID, EgammaBdtHZZBasedID;
  EgammaBdtHWWBasedID.push_back("WP80");
  EgammaBdtHWWBasedID.push_back("WP80EA");
  EgammaBdtHWWBasedID.push_back("IsoWP80");
  EgammaBdtHWWBasedID.push_back("IsoWP80EA");
  EgammaBdtHWWBasedID.push_back("IdWP80");

  // no need to repeat the iso only
  EgammaBdtHZZBasedID.push_back("WP80HWWFR");
  EgammaBdtHZZBasedID.push_back("WP80");
  EgammaBdtHZZBasedID.push_back("WP80ch"); // ch is with tracker isolation only
  EgammaBdtHZZBasedID.push_back("WP95"); 
  EgammaBdtHZZBasedID.push_back("WP95ch"); // ch is with tracker isolation only
  EgammaBdtHZZBasedID.push_back("WP70x80"); // low pT - high pT different WPs
  EgammaBdtHZZBasedID.push_back("WP70x80ch"); // low pT - high pT different WPs (ch isolation only)

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
  
  std::vector<TH1F*> BdtHWWEta;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"Eta",   "HWW BDT ID #eta", 4, LowerEta);
    BdtHWWEta.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZEta;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"Eta",   "HZZ BDT ID #eta", 4, LowerEta); 
    BdtHZZEta.push_back(aHisto);
  }

  // eta, high pT
  std::vector<TH1F*> BdtHWWEtaHighPt;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"EtaHighPt",   "HWW BDT ID #eta", 4, LowerEta);   
    BdtHWWEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZEtaHighPt;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"EtaHighPt",   "HZZ BDT ID #eta", 4, LowerEta); 
    BdtHZZEtaHighPt.push_back(aHisto);
  }

  // eta, low pT
  std::vector<TH1F*> BdtHWWEtaLowPt;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"EtaLowPt",   "HWW BDT ID #eta", 4, LowerEta);  
    BdtHWWEtaLowPt.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZEtaLowPt;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"EtaLowPt",   "HZZ BDT ID #eta", 4, LowerEta);
    BdtHZZEtaLowPt.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs pT
  Float_t LowerPt[9];
  LowerPt[0]=10;
  LowerPt[1]=15;
  LowerPt[2]=20;  
  LowerPt[3]=25;
  LowerPt[4]=30;  
  LowerPt[5]=35;  
  LowerPt[6]=40;  
  LowerPt[7]=45;  
  LowerPt[8]=50;  
  TH1F *RecoPt         = new TH1F( "RecoPt",          "reconstructed p_{T} (GeV)", 8, LowerPt);
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", 8, LowerPt);
  TH1F *RecoPtBarrel1  = new TH1F( "RecoPtBarrel1",   "reconstructed p_{T} (GeV)", 8, LowerPt);    // this is done to subsplit 
  TH1F *RecoPtBarrel2  = new TH1F( "RecoPtBarrel2",   "reconstructed p_{T} (GeV)", 8, LowerPt);    // barrel and endcap in 2 parts 
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", 8, LowerPt);
  TH1F *RecoPtEndcap1  = new TH1F( "RecoPtEndcap1",   "reconstructed p_{T} (GeV)", 8, LowerPt);
  TH1F *RecoPtEndcap2  = new TH1F( "RecoPtEndcap2",   "reconstructed p_{T} (GeV)", 8, LowerPt);

  std::vector<TH1F*> BdtHWWPt;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"Pt",   "HWW BDT ID #eta", 8, LowerPt );
    BdtHWWPt.push_back(aHisto);
  }
  
  std::vector<TH1F*> BdtHZZPt;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"Pt", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZPt.push_back(aHisto);
  }


  // to have the full picture in the barrel
  std::vector<TH1F*> BdtHWWPtBarrel;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PtBarrel", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWPtBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZPtBarrel;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PtBarrel", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZPtBarrel.push_back(aHisto);
  }


  // to compute the FR in 4 eta bins ( barrel )
  std::vector<TH1F*> BdtHWWPtBarrel1;
  std::vector<TH1F*> BdtHWWPtBarrel2;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PtBarrel1", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PtBarrel2", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWPtBarrel2.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZPtBarrel1;
  std::vector<TH1F*> BdtHZZPtBarrel2;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PtBarrel1", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PtBarrel2", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZPtBarrel2.push_back(aHisto);
  }

  // to have the full picture in the endcap
  std::vector<TH1F*> BdtHWWPtEndcap;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PtEndcap", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWPtEndcap.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZPtEndcap;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PtEndcap", "HZZ BDT ID #eta",   8, LowerPt);
    BdtHZZPtEndcap.push_back(aHisto);
  }


  // to compute the FR in 4 eta bins (endcap)
  std::vector<TH1F*> BdtHWWPtEndcap1;
  std::vector<TH1F*> BdtHWWPtEndcap2;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PtEndcap1", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PtEndcap2", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWPtEndcap2.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZPtEndcap1;
  std::vector<TH1F*> BdtHZZPtEndcap2;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PtEndcap1", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PtEndcap2", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZPtEndcap2.push_back(aHisto);
  }


  // -----------------------------------------------------------------------
  // study vs PU
  Float_t LowerPU[9];
  LowerPU[0] = 1;
  LowerPU[1] = 4;
  for(int i=2;i<6;++i) LowerPU[i]=i+3;
  LowerPU[6] = 10;
  LowerPU[7] = 15;
  LowerPU[8] = 25;

  TH1F *RecoPU         = new TH1F( "RecoPU",         "reconstructed nPU", 8, LowerPU );
  TH1F *RecoPUBarrel1  = new TH1F( "RecoPUBarrel1",  "reconstructed nPU", 8, LowerPU);    // this is done to subsplit 
  TH1F *RecoPUBarrel2  = new TH1F( "RecoPUBarrel2",  "reconstructed nPU", 8, LowerPU);    // barrel and endcap in 2 parts 
  TH1F *RecoPUEndcap1  = new TH1F( "RecoPUEndcap1",  "reconstructed nPU", 8, LowerPU);
  TH1F *RecoPUEndcap2  = new TH1F( "RecoPUEndcap2",  "reconstructed nPU", 8, LowerPU);

  std::vector<TH1F*> BdtHWWPU;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PU", "HWW BDT ID nPU",  8, LowerPU );
    BdtHWWPU.push_back(aHisto);
  }
  
  std::vector<TH1F*> BdtHZZPU;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PU", "HZZ BDT ID nPU",  8, LowerPU );
    BdtHZZPU.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins ( barrel )
  std::vector<TH1F*> BdtHWWPUBarrel1;
  std::vector<TH1F*> BdtHWWPUBarrel2;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PUBarrel1", "HWW BDT ID #eta",   8, LowerPU );
    BdtHWWPUBarrel1.push_back(aHisto);
    aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PUBarrel2", "HWW BDT ID #eta",   8, LowerPU );
    BdtHWWPUBarrel2.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZPUBarrel1;
  std::vector<TH1F*> BdtHZZPUBarrel2;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PUBarrel1", "HZZ BDT ID #eta",   8, LowerPU );
    BdtHZZPUBarrel1.push_back(aHisto);
    aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PUBarrel2", "HZZ BDT ID #eta",   8, LowerPU );
    BdtHZZPUBarrel2.push_back(aHisto);
  }

  // to compute the FR in 4 eta bins (endcap)
  std::vector<TH1F*> BdtHWWPUEndcap1;
  std::vector<TH1F*> BdtHWWPUEndcap2;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PUEndcap1", "HWW BDT ID #eta",   8, LowerPU );
    BdtHWWPUEndcap1.push_back(aHisto);
    aHisto = new TH1F( "BdtHWW"+TString(EgammaBdtHWWBasedID[i])+"PUEndcap2", "HWW BDT ID #eta",   8, LowerPU );
    BdtHWWPUEndcap2.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZPUEndcap1;
  std::vector<TH1F*> BdtHZZPUEndcap2;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PUEndcap1", "HZZ BDT ID #eta",   8, LowerPU );
    BdtHZZPUEndcap1.push_back(aHisto);
    aHisto = new TH1F( "BdtHZZ"+TString(EgammaBdtHZZBasedID[i])+"PUEndcap2", "HZZ BDT ID #eta",   8, LowerPU );
    BdtHZZPUEndcap2.push_back(aHisto);
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
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // NOT EWK CORRECTED: STOP AT 35 GEV!!!
    if(pt>35) continue;

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

    if (isInEB) RecoPtBarrel -> Fill(etFake);  //, theWeight);
    if (isInEE) RecoPtEndcap -> Fill(etFake);  //, theWeight);

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
      BdtHWWEta[kWP80]->Fill(etaFake);
      BdtHWWPt[kWP80] ->Fill(etFake);
      BdtHWWPU[kWP80] ->Fill(vertices);
      if (highPt) BdtHWWEtaHighPt[kWP80]->Fill(etaFake);
      if (lowPt)  BdtHWWEtaLowPt[kWP80] ->Fill(etaFake);
      if (isInEB) BdtHWWPtBarrel[kWP80] ->Fill(etFake);
      if (isInEE) BdtHWWPtEndcap[kWP80] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHWWPtBarrel1[kWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel1[kWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHWWPtBarrel2[kWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel2[kWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	BdtHWWPtEndcap1[kWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUEndcap1[kWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	BdtHWWPtEndcap2[kWP80] -> Fill(etFake);   //, theWeight);     
	BdtHWWPUEndcap2[kWP80] -> Fill(vertices); //, theWeight);     
      }
    }
    // === WP80 HWW BDT + EA corr Iso ===
    if(passID(kIsoEACorr) && passID(kBDTHWW2011_withIP)) {
      BdtHWWEta[kWP80EA]->Fill(etaFake);
      BdtHWWPt[kWP80EA] ->Fill(etFake);
      BdtHWWPU[kWP80EA] ->Fill(vertices);
      if (highPt) BdtHWWEtaHighPt[kWP80EA]->Fill(etaFake);
      if (lowPt)  BdtHWWEtaLowPt[kWP80EA] ->Fill(etaFake);
      if (isInEB) BdtHWWPtBarrel[kWP80EA] ->Fill(etFake);
      if (isInEE) BdtHWWPtEndcap[kWP80EA] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHWWPtBarrel1[kWP80EA] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel1[kWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHWWPtBarrel2[kWP80EA] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel2[kWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	BdtHWWPtEndcap1[kWP80EA] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUEndcap1[kWP80EA] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	BdtHWWPtEndcap2[kWP80EA] -> Fill(etFake);   //, theWeight);     
	BdtHWWPUEndcap2[kWP80EA] -> Fill(vertices); //, theWeight);     
      }
    }
    // === HWW 2011 Iso ===
    if(passID(kIsoHWW2011)) {
      BdtHWWEta[kIsoWP80]->Fill(etaFake);
      BdtHWWPt[kIsoWP80] ->Fill(etFake);
      BdtHWWPU[kIsoWP80] ->Fill(vertices);
      if (highPt) BdtHWWEtaHighPt[kIsoWP80]->Fill(etaFake);
      if (lowPt)  BdtHWWEtaLowPt[kIsoWP80] ->Fill(etaFake);
      if (isInEB) BdtHWWPtBarrel[kIsoWP80] ->Fill(etFake);
      if (isInEE) BdtHWWPtEndcap[kIsoWP80] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHWWPtBarrel1[kIsoWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel1[kIsoWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHWWPtBarrel2[kIsoWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel2[kIsoWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	BdtHWWPtEndcap1[kIsoWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUEndcap1[kIsoWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	BdtHWWPtEndcap2[kIsoWP80] -> Fill(etFake);   //, theWeight);     
	BdtHWWPUEndcap2[kIsoWP80] -> Fill(vertices); //, theWeight);     
      }
    }
    // === EA corr Iso ===
    if(passID(kIsoEACorr)) {
      BdtHWWEta[kIsoWP80EA]->Fill(etaFake);
      BdtHWWPt[kIsoWP80EA] ->Fill(etFake);
      BdtHWWPU[kIsoWP80EA] ->Fill(vertices);
      if (highPt) BdtHWWEtaHighPt[kIsoWP80EA]->Fill(etaFake);
      if (lowPt)  BdtHWWEtaLowPt[kIsoWP80EA] ->Fill(etaFake);
      if (isInEB) BdtHWWPtBarrel[kIsoWP80EA] ->Fill(etFake);
      if (isInEE) BdtHWWPtEndcap[kIsoWP80EA] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHWWPtBarrel1[kIsoWP80EA] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel1[kIsoWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHWWPtBarrel2[kIsoWP80EA] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel2[kIsoWP80EA] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	BdtHWWPtEndcap1[kIsoWP80EA] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUEndcap1[kIsoWP80EA] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	BdtHWWPtEndcap2[kIsoWP80EA] -> Fill(etFake);   //, theWeight);     
	BdtHWWPUEndcap2[kIsoWP80EA] -> Fill(vertices); //, theWeight);     
      }
    }
    // === HWW 2011 BDT ===
    if(passID(kBDTHWW2011_withIP)) {
      BdtHWWEta[kIdWP80]->Fill(etaFake);
      BdtHWWPt[kIdWP80] ->Fill(etFake);
      BdtHWWPU[kIdWP80] ->Fill(vertices);
      if (highPt) BdtHWWEtaHighPt[kIdWP80]->Fill(etaFake);
      if (lowPt)  BdtHWWEtaLowPt[kIdWP80] ->Fill(etaFake);
      if (isInEB) BdtHWWPtBarrel[kIdWP80] ->Fill(etFake);
      if (isInEE) BdtHWWPtEndcap[kIdWP80] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHWWPtBarrel1[kIdWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel1[kIdWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHWWPtBarrel2[kIdWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUBarrel2[kIdWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) { 
	BdtHWWPtEndcap1[kIdWP80] -> Fill(etFake);   //, theWeight);      
	BdtHWWPUEndcap1[kIdWP80] -> Fill(vertices); //, theWeight);      
      }
      if (etaRegion==4) {
	BdtHWWPtEndcap2[kIdWP80] -> Fill(etFake);   //, theWeight);     
	BdtHWWPUEndcap2[kIdWP80] -> Fill(vertices); //, theWeight);     
      }
    }

    // === WP80 HZZ BDT + EA corr Iso === : this is tuned to give the same fake-rate as HWW 2011
    if(passID(kIsoEACorr) && passID(kBDTHZZ_withIP)) {
      BdtHZZEta[kHzzWP80HWWFR]->Fill(etaFake);
      BdtHZZPt[kHzzWP80HWWFR] ->Fill(etFake);
      BdtHZZPU[kHzzWP80HWWFR] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP80HWWFR]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP80HWWFR] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP80HWWFR] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP80HWWFR] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP80HWWFR] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP80HWWFR] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP80HWWFR] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP80HWWFR] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP80HWWFR] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP80HWWFR] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP80HWWFR] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP80HWWFR] -> Fill(vertices); //, theWeight);
      }
    }

    // === WP80 HZZ reoptimized cut ===
    if(aSel.output(etFake,etaFake,bdthzz,combPFIsoHZZ/etFake,HZZEleIDSelector::kWP80)) {
      BdtHZZEta[kHzzWP80]->Fill(etaFake);
      BdtHZZPt[kHzzWP80] ->Fill(etFake);
      BdtHZZPU[kHzzWP80] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP80]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP80] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP80] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP80] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP80] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP80] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP80] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP80] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP80] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP80] -> Fill(vertices); //, theWeight);
      }
    }

    // === WP80 HZZ reoptimized cut (charged isolation) ===
    if(aSel.output(etFake,etaFake,bdthzz,chaPFIso/etFake,HZZEleIDSelector::kWP80ChIso)) {
      BdtHZZEta[kHzzWP80ch]->Fill(etaFake);
      BdtHZZPt[kHzzWP80ch] ->Fill(etFake);
      BdtHZZPU[kHzzWP80ch] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP80ch]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP80ch] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP80ch] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP80ch] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP80ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP80ch] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP80ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP80ch] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP80ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP80ch] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP80ch] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP80ch] -> Fill(vertices); //, theWeight);
      }
    }

    // === WP95 HZZ reoptimized cut ===
    if(aSel.output(etFake,etaFake,bdthzz,combPFIsoHZZ/etFake,HZZEleIDSelector::kWP95)) {
      BdtHZZEta[kHzzWP95]->Fill(etaFake);
      BdtHZZPt[kHzzWP95] ->Fill(etFake);
      BdtHZZPU[kHzzWP95] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP95]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP95] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP95] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP95] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP95] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP95] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP95] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP95] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP95] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP95] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP95] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP95] -> Fill(vertices); //, theWeight);
      }
    }

    // === WP95 HZZ reoptimized cut (charged isolation) ===
    if(aSel.output(etFake,etaFake,bdthzz,chaPFIso/etFake,HZZEleIDSelector::kWP95ChIso)) {
      BdtHZZEta[kHzzWP95ch]->Fill(etaFake);
      BdtHZZPt[kHzzWP95ch] ->Fill(etFake);
      BdtHZZPU[kHzzWP95ch] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP95ch]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP95ch] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP95ch] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP95ch] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP95ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP95ch] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP95ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP95ch] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP95ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP95ch] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP95ch] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP95ch] -> Fill(vertices); //, theWeight);
      }
    }

    // === WP70 (pt<20 GeV) && WP80 (pt>20 GeV) HZZ reoptimized cut ===
    if((etFake<20 && aSel.output(etFake,etaFake,bdthzz,combPFIsoHZZ/etFake,HZZEleIDSelector::kWP70)) ||
       (etFake>=20 && aSel.output(etFake,etaFake,bdthzz,combPFIsoHZZ/etFake,HZZEleIDSelector::kWP80))) {
      BdtHZZEta[kHzzWP70x80]->Fill(etaFake);
      BdtHZZPt[kHzzWP70x80] ->Fill(etFake);
      BdtHZZPU[kHzzWP70x80] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP70x80]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP70x80] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP70x80] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP70x80] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP70x80] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP70x80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP70x80] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP70x80] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP70x80] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP70x80] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP70x80] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP70x80] -> Fill(vertices); //, theWeight);
      }
    }

    // === WP70 (pt<20 GeV) && WP80 (pt>20 GeV) HZZ reoptimized cut (charged isolation) ===
    if((etFake<20 && aSel.output(etFake,etaFake,bdthzz,combPFIsoHZZ/etFake,HZZEleIDSelector::kWP70ChIso)) ||
       (etFake>=20 && aSel.output(etFake,etaFake,bdthzz,combPFIsoHZZ/etFake,HZZEleIDSelector::kWP80ChIso))) {
      BdtHZZEta[kHzzWP70x80ch]->Fill(etaFake);
      BdtHZZPt[kHzzWP70x80ch] ->Fill(etFake);
      BdtHZZPU[kHzzWP70x80ch] ->Fill(vertices);
      if (highPt) BdtHZZEtaHighPt[kHzzWP70x80ch]->Fill(etaFake);
      if (lowPt)  BdtHZZEtaLowPt[kHzzWP70x80ch] ->Fill(etaFake);
      if (isInEB) BdtHZZPtBarrel[kHzzWP70x80ch] ->Fill(etFake);
      if (isInEE) BdtHZZPtEndcap[kHzzWP70x80ch] ->Fill(etFake);
      if (etaRegion==1) {
	BdtHZZPtBarrel1[kHzzWP70x80ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel1[kHzzWP70x80ch] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==2) {
	BdtHZZPtBarrel2[kHzzWP70x80ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUBarrel2[kHzzWP70x80ch] -> Fill(vertices); //, theWeight);
      }
      if (etaRegion==3) {
	BdtHZZPtEndcap1[kHzzWP70x80ch] -> Fill(etFake);   //, theWeight);      
	BdtHZZPUEndcap1[kHzzWP70x80ch] -> Fill(vertices); //, theWeight); 
      }
      if (etaRegion==4) {
	BdtHZZPtEndcap2[kHzzWP70x80ch] -> Fill(etFake);   //, theWeight);     
	BdtHZZPUEndcap2[kHzzWP70x80ch] -> Fill(vertices); //, theWeight);
      }
    }


  } // loop over events


  // saving efficiency histos
  // === as a function of eta ===
  char filename[500];
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  ElectronEffEta.AddNumerator(RecoEta);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffEta.AddNumerator(BdtHWWEta[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffEta.AddNumerator(BdtHZZEta[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(BdtHWWEtaHighPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(BdtHZZEtaHighPt[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(BdtHWWEtaLowPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(BdtHZZEtaLowPt[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPt.AddNumerator(BdtHWWPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPt.AddNumerator(BdtHZZPt[icut]);
  }

  ElectronEffPt.SetDenominator(RecoPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("fake rate vs pT");
  ElectronEffPt.SetXaxisTitle("electron pT");
  ElectronEffPt.SetYaxisTitle("Fake rate");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

  sprintf(filename,"%s-EleMisidPtBarrel1.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel1(filename);
  ElectronEffPtBarrel1.AddNumerator(RecoPtBarrel1);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtBarrel1.AddNumerator(BdtHWWPtBarrel1[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtBarrel1.AddNumerator(BdtHZZPtBarrel1[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtBarrel2.AddNumerator(BdtHWWPtBarrel2[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtBarrel2.AddNumerator(BdtHZZPtBarrel2[icut]);
  }

  ElectronEffPtBarrel2.SetDenominator(RecoPtBarrel2);
  ElectronEffPtBarrel2.ComputeEfficiencies();
  ElectronEffPtBarrel2.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel2.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel2.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel2.SetYaxisMin(0.0);
  ElectronEffPtBarrel2.Write();

  sprintf(filename,"%s-EleMisidPtEndcap1.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap1(filename);
  ElectronEffPtEndcap1.AddNumerator(RecoPtEndcap1);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtEndcap1.AddNumerator(BdtHWWPtEndcap1[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtEndcap1.AddNumerator(BdtHZZPtEndcap1[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtEndcap2.AddNumerator(BdtHWWPtEndcap2[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtEndcap2.AddNumerator(BdtHZZPtEndcap2[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPU.AddNumerator(BdtHWWPU[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPU.AddNumerator(BdtHZZPU[icut]);
  }

  ElectronEffPU.SetDenominator(RecoPU);
  ElectronEffPU.ComputeEfficiencies();
  ElectronEffPU.SetTitle("fake rate vs PU");
  ElectronEffPU.SetXaxisTitle("# vertices");
  ElectronEffPU.SetYaxisTitle("Fake rate");
  ElectronEffPU.SetYaxisMin(0.0);
  ElectronEffPU.Write();

  sprintf(filename,"%s-EleMisidPUBarrel1.root",outname);
  EfficiencyEvaluator ElectronEffPUBarrel1(filename);
  ElectronEffPUBarrel1.AddNumerator(RecoPUBarrel1);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPUBarrel1.AddNumerator(BdtHWWPUBarrel1[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPUBarrel1.AddNumerator(BdtHZZPUBarrel1[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPUBarrel2.AddNumerator(BdtHWWPUBarrel2[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPUBarrel2.AddNumerator(BdtHZZPUBarrel2[icut]);
  }

  ElectronEffPUBarrel2.SetDenominator(RecoPUBarrel2);
  ElectronEffPUBarrel2.ComputeEfficiencies();
  ElectronEffPUBarrel2.SetTitle("fake rate vs vertices");
  ElectronEffPUBarrel2.SetXaxisTitle("# vertices");
  ElectronEffPUBarrel2.SetYaxisTitle("Fake rate");
  ElectronEffPUBarrel2.SetYaxisMin(0.0);
  ElectronEffPUBarrel2.Write();

  sprintf(filename,"%s-EleMisidPUEndcap1.root",outname);
  EfficiencyEvaluator ElectronEffPUEndcap1(filename);
  ElectronEffPUEndcap1.AddNumerator(RecoPUEndcap1);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPUEndcap1.AddNumerator(BdtHWWPUEndcap1[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPUEndcap1.AddNumerator(BdtHZZPUEndcap1[icut]);
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
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPUEndcap2.AddNumerator(BdtHWWPUEndcap2[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPUEndcap2.AddNumerator(BdtHZZPUEndcap2[icut]);
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
      if(fabs(eta) <  1.0) return (combIso/pt < 0.22); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.21);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.12);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.11);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.074);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.053);
      if(fabs(eta) >=  2.4) return (combIso/pt < 0.010);
    } else {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.20); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.24);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.13);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.083);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.01);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.027);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.035);
    }
  }

  if(type == kBDTHWW2011_withIP) {
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.884);
  }

  if(type == kBDTHWW2011_noIP) {
    // warning: no reoptimized WP!
    if(pt < 20 && fabs(eta) < 1.0) return (bdthwwnoip > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthwwnoip > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthwwnoip > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthwwnoip > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthwwnoip > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthwwnoip > 0.884);
  }

  if(type == kBDTHZZ_withIP) {
    // WP with same fake rate as HWW with IP
    if(pt < 20 && fabs(eta) < 1.0) return (bdthzz > 0.099);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.105);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.12);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthzz > 0.08);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.091);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.086);
  }

  if(type == kBDTHZZ_noIP) {
    // not optimized WP!
    if(pt < 20 && fabs(eta) < 1.0) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.091);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthzz > 0.064);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.071);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.067);
  }
  return false;
}
