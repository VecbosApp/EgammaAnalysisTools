#define estimateFakeRate_cxx
#include "estimateFakeRate.h"
#include "CommonTools/include/EfficiencyEvaluator.hh"
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


  std::vector<TString> EgammaBdtHWWBasedID, EgammaBdtHZZBasedID;
  EgammaBdtHWWBasedID.push_back("WP80");
  EgammaBdtHWWBasedID.push_back("WP80EA");
  EgammaBdtHWWBasedID.push_back("IsoWP80");
  EgammaBdtHWWBasedID.push_back("IsoWP80EA");
  EgammaBdtHWWBasedID.push_back("IdWP80");

  // no need to repeat the iso only
  EgammaBdtHZZBasedID.push_back("WP80");
  EgammaBdtHZZBasedID.push_back("WP80EA");
  EgammaBdtHZZBasedID.push_back("IsoWP80"); // same as WW, not filled. Dummy
  EgammaBdtHZZBasedID.push_back("IsoWP80EA"); // same as WW, not filled. Dummy
  EgammaBdtHZZBasedID.push_back("IdWP80");

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
  
  std::vector<TH1F*> BdtHWWIdEta;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"Eta",   "HWW BDT ID #eta", 4, LowerEta);
    BdtHWWIdEta.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdEta;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"Eta",   "HZZ BDT ID #eta", 4, LowerEta); 
    BdtHZZIdEta.push_back(aHisto);
  }

  // eta, high pT
  std::vector<TH1F*> BdtHWWIdEtaHighPt;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"EtaHighPt",   "HWW BDT ID #eta", 4, LowerEta);   
    BdtHWWIdEtaHighPt.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdEtaHighPt;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"EtaHighPt",   "HZZ BDT ID #eta", 4, LowerEta); 
    BdtHZZIdEtaHighPt.push_back(aHisto);
  }

  // eta, low pT
  std::vector<TH1F*> BdtHWWIdEtaLowPt;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"EtaLowPt",   "HWW BDT ID #eta", 4, LowerEta);  
    BdtHWWIdEtaLowPt.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdEtaLowPt;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"EtaLowPt",   "HZZ BDT ID #eta", 4, LowerEta);
    BdtHZZIdEtaLowPt.push_back(aHisto);
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

  std::vector<TH1F*> BdtHWWIdPt;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"Pt",   "HWW BDT ID #eta", 8, LowerPt );
    BdtHWWIdPt.push_back(aHisto);
  }
  
  std::vector<TH1F*> BdtHZZIdPt;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"Pt", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZIdPt.push_back(aHisto);
  }


  // to have the full picture in the barrel
  std::vector<TH1F*> BdtHWWIdPtBarrel;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PtBarrel", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWIdPtBarrel.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdPtBarrel;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PtBarrel", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZIdPtBarrel.push_back(aHisto);
  }


  // to compute the FR in 4 eta bins ( barrel )
  std::vector<TH1F*> BdtHWWIdPtBarrel1;
  std::vector<TH1F*> BdtHWWIdPtBarrel2;

  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PtBarrel1", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWIdPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PtBarrel2", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWIdPtBarrel2.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdPtBarrel1;
  std::vector<TH1F*> BdtHZZIdPtBarrel2;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PtBarrel1", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZIdPtBarrel1.push_back(aHisto);
    aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PtBarrel2", "HZZ BDT ID #eta",   8, LowerPt );
    BdtHZZIdPtBarrel2.push_back(aHisto);
  }

  // to have the full picture in the endcap
  std::vector<TH1F*> BdtHWWIdPtEndcap;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PtEndcap", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWIdPtEndcap.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdPtEndcap;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PtEndcap", "HZZ BDT ID #eta",   8, LowerPt);
    BdtHZZIdPtEndcap.push_back(aHisto);
  }


  // to compute the FR in 4 eta bins (endcap)
  std::vector<TH1F*> BdtHWWIdPtEndcap1;
  std::vector<TH1F*> BdtHWWIdPtEndcap2;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PtEndcap1", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWIdPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PtEndcap2", "HWW BDT ID #eta",   8, LowerPt );
    BdtHWWIdPtEndcap2.push_back(aHisto);
  }

  std::vector<TH1F*> BdtHZZIdPtEndcap1;
  std::vector<TH1F*> BdtHZZIdPtEndcap2;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PtEndcap1", "HZZ BDT ID #eta",   8, LowerPt);
    BdtHZZIdPtEndcap1.push_back(aHisto);
    aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PtEndcap2", "HZZ BDT ID #eta",   8, LowerPt);
    BdtHZZIdPtEndcap2.push_back(aHisto);
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

  std::vector<TH1F*> BdtHWWIdPU;
  for (int i=0;i<(int)EgammaBdtHWWBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHWWId"+TString(EgammaBdtHWWBasedID[i])+"PU", "HWW BDT ID nPU",  8, LowerPU );
    BdtHWWIdPU.push_back(aHisto);
  }
  
  std::vector<TH1F*> BdtHZZIdPU;
  for (int i=0;i<(int)EgammaBdtHZZBasedID.size();++i) {
    TH1F* aHisto = new TH1F( "BdtHZZId"+TString(EgammaBdtHZZBasedID[i])+"PU", "HZZ BDT ID nPU",  8, LowerPU );
    BdtHZZIdPU.push_back(aHisto);
  }

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
    float etaFake = eta;
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

    if (etaRegion==1) RecoPtBarrel1 -> Fill(etFake);  //, theWeight);      
    if (etaRegion==2) RecoPtBarrel2 -> Fill(etFake);  //, theWeight);      
    if (etaRegion==3) RecoPtEndcap1 -> Fill(etFake);  //, theWeight);      
    if (etaRegion==4) RecoPtEndcap2 -> Fill(etFake);  //, theWeight);      

    // fill the numerator(s)
    // === WP80 HWW 2011: BDT + ISO ===
    if(passID(kIsoHWW2011) && passID(kBDTHWW2011_withIP)) {
      BdtHWWIdEta[kWP80]->Fill(etaFake);
      BdtHWWIdPt[kWP80] ->Fill(etFake);
      BdtHWWIdPU[kWP80] ->Fill(vertices);
      if (highPt) BdtHWWIdEtaHighPt[kWP80]->Fill(etaFake);
      if (lowPt)  BdtHWWIdEtaLowPt[kWP80] ->Fill(etaFake);
      if (isInEB) BdtHWWIdPtBarrel[kWP80] ->Fill(etFake);
      if (isInEE) BdtHWWIdPtEndcap[kWP80] ->Fill(etFake);
      if (etaRegion==1) BdtHWWIdPtBarrel1[kWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHWWIdPtBarrel2[kWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHWWIdPtEndcap1[kWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHWWIdPtEndcap2[kWP80] -> Fill(etFake);  //, theWeight);     
    }
    // === WP80 HWW BDT + EA corr Iso ===
    if(passID(kIsoEACorr) && passID(kBDTHWW2011_withIP)) {
      BdtHWWIdEta[kWP80EA]->Fill(etaFake);
      BdtHWWIdPt[kWP80EA] ->Fill(etFake);
      BdtHWWIdPU[kWP80EA] ->Fill(vertices);
      if (highPt) BdtHWWIdEtaHighPt[kWP80EA]->Fill(etaFake);
      if (lowPt)  BdtHWWIdEtaLowPt[kWP80EA] ->Fill(etaFake);
      if (isInEB) BdtHWWIdPtBarrel[kWP80EA] ->Fill(etFake);
      if (isInEE) BdtHWWIdPtEndcap[kWP80EA] ->Fill(etFake);
      if (etaRegion==1) BdtHWWIdPtBarrel1[kWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHWWIdPtBarrel2[kWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHWWIdPtEndcap1[kWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHWWIdPtEndcap2[kWP80EA] -> Fill(etFake);  //, theWeight);     
    }
    // === HWW 2011 Iso ===
    if(passID(kIsoHWW2011)) {
      BdtHWWIdEta[kIsoWP80]->Fill(etaFake);
      BdtHWWIdPt[kIsoWP80] ->Fill(etFake);
      BdtHWWIdPU[kIsoWP80] ->Fill(vertices);
      if (highPt) BdtHWWIdEtaHighPt[kIsoWP80]->Fill(etaFake);
      if (lowPt)  BdtHWWIdEtaLowPt[kIsoWP80] ->Fill(etaFake);
      if (isInEB) BdtHWWIdPtBarrel[kIsoWP80] ->Fill(etFake);
      if (isInEE) BdtHWWIdPtEndcap[kIsoWP80] ->Fill(etFake);
      if (etaRegion==1) BdtHWWIdPtBarrel1[kIsoWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHWWIdPtBarrel2[kIsoWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHWWIdPtEndcap1[kIsoWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHWWIdPtEndcap2[kIsoWP80] -> Fill(etFake);  //, theWeight);     
    }
    // === EA corr Iso ===
    if(passID(kIsoEACorr)) {
      BdtHWWIdEta[kIsoWP80EA]->Fill(etaFake);
      BdtHWWIdPt[kIsoWP80EA] ->Fill(etFake);
      BdtHWWIdPU[kIsoWP80EA] ->Fill(vertices);
      if (highPt) BdtHWWIdEtaHighPt[kIsoWP80EA]->Fill(etaFake);
      if (lowPt)  BdtHWWIdEtaLowPt[kIsoWP80EA] ->Fill(etaFake);
      if (isInEB) BdtHWWIdPtBarrel[kIsoWP80EA] ->Fill(etFake);
      if (isInEE) BdtHWWIdPtEndcap[kIsoWP80EA] ->Fill(etFake);
      if (etaRegion==1) BdtHWWIdPtBarrel1[kIsoWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHWWIdPtBarrel2[kIsoWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHWWIdPtEndcap1[kIsoWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHWWIdPtEndcap2[kIsoWP80EA] -> Fill(etFake);  //, theWeight);     
    }
    // === HWW 2011 BDT ===
    if(passID(kBDTHWW2011_withIP)) {
      BdtHWWIdEta[kIdWP80]->Fill(etaFake);
      BdtHWWIdPt[kIdWP80] ->Fill(etFake);
      BdtHWWIdPU[kIdWP80] ->Fill(vertices);
      if (highPt) BdtHWWIdEtaHighPt[kIdWP80]->Fill(etaFake);
      if (lowPt)  BdtHWWIdEtaLowPt[kIdWP80] ->Fill(etaFake);
      if (isInEB) BdtHWWIdPtBarrel[kIdWP80] ->Fill(etFake);
      if (isInEE) BdtHWWIdPtEndcap[kIdWP80] ->Fill(etFake);
      if (etaRegion==1) BdtHWWIdPtBarrel1[kIdWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHWWIdPtBarrel2[kIdWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHWWIdPtEndcap1[kIdWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHWWIdPtEndcap2[kIdWP80] -> Fill(etFake);  //, theWeight);     
    }
    // === WP80 HZZ 2012 + HWW 2011 Iso ===
    if(passID(kIsoHWW2011) && passID(kBDTHZZ_withIP)) {
      BdtHZZIdEta[kWP80]->Fill(etaFake);
      BdtHZZIdPt[kWP80] ->Fill(etFake);
      BdtHZZIdPU[kWP80] ->Fill(vertices);
      if (highPt) BdtHZZIdEtaHighPt[kWP80]->Fill(etaFake);
      if (lowPt)  BdtHZZIdEtaLowPt[kWP80] ->Fill(etaFake);
      if (isInEB) BdtHZZIdPtBarrel[kWP80] ->Fill(etFake);
      if (isInEE) BdtHZZIdPtEndcap[kWP80] ->Fill(etFake);
      if (etaRegion==1) BdtHZZIdPtBarrel1[kWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHZZIdPtBarrel2[kWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHZZIdPtEndcap1[kWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHZZIdPtEndcap2[kWP80] -> Fill(etFake);  //, theWeight);     
    }
    // === WP80 HZZ BDT + EA corr Iso ===
    if(passID(kIsoEACorr) && passID(kBDTHZZ_withIP)) {
      BdtHZZIdEta[kWP80EA]->Fill(etaFake);
      BdtHZZIdPt[kWP80EA] ->Fill(etFake);
      BdtHZZIdPU[kWP80EA] ->Fill(vertices);
      if (highPt) BdtHZZIdEtaHighPt[kWP80EA]->Fill(etaFake);
      if (lowPt)  BdtHZZIdEtaLowPt[kWP80EA] ->Fill(etaFake);
      if (isInEB) BdtHZZIdPtBarrel[kWP80EA] ->Fill(etFake);
      if (isInEE) BdtHZZIdPtEndcap[kWP80EA] ->Fill(etFake);
      if (etaRegion==1) BdtHZZIdPtBarrel1[kWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHZZIdPtBarrel2[kWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHZZIdPtEndcap1[kWP80EA] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHZZIdPtEndcap2[kWP80EA] -> Fill(etFake);  //, theWeight);     
    }
    // === HZZ 2011 BDT ===
    if(passID(kBDTHZZ_withIP)) {
      BdtHZZIdEta[kIdWP80]->Fill(etaFake);
      BdtHZZIdPt[kIdWP80] ->Fill(etFake);
      BdtHZZIdPU[kIdWP80] ->Fill(vertices);
      if (highPt) BdtHZZIdEtaHighPt[kIdWP80]->Fill(etaFake);
      if (lowPt)  BdtHZZIdEtaLowPt[kIdWP80] ->Fill(etaFake);
      if (isInEB) BdtHZZIdPtBarrel[kIdWP80] ->Fill(etFake);
      if (isInEE) BdtHZZIdPtEndcap[kIdWP80] ->Fill(etFake);
      if (etaRegion==1) BdtHZZIdPtBarrel1[kIdWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==2) BdtHZZIdPtBarrel2[kIdWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==3) BdtHZZIdPtEndcap1[kIdWP80] -> Fill(etFake);  //, theWeight);      
      if (etaRegion==4) BdtHZZIdPtEndcap2[kIdWP80] -> Fill(etFake);  //, theWeight);     
    }
  }


  // saving efficiency histos
  // === as a function of eta ===
  char filename[500];
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  ElectronEffEta.AddNumerator(RecoEta);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffEta.AddNumerator(BdtHWWIdEta[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffEta.AddNumerator(BdtHZZIdEta[icut]);
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
  ElectronEffEtaHighPt.AddNumerator(RecoEta);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(BdtHWWIdEtaHighPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(BdtHZZIdEtaHighPt[icut]);
  }

  ElectronEffEtaHighPt.SetDenominator(RecoEta);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs pT");
  ElectronEffEtaHighPt.SetXaxisTitle("electron p_{T}");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-EleMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  ElectronEffEtaLowPt.AddNumerator(RecoEta);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(BdtHWWIdEtaLowPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(BdtHZZIdEtaLowPt[icut]);
  }

  ElectronEffEtaLowPt.SetDenominator(RecoEta);
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
    ElectronEffPt.AddNumerator(BdtHWWIdPt[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPt.AddNumerator(BdtHZZIdPt[icut]);
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
  ElectronEffPtBarrel1.AddNumerator(RecoPt);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtBarrel1.AddNumerator(BdtHWWIdPtBarrel1[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtBarrel1.AddNumerator(BdtHZZIdPtBarrel1[icut]);
  }

  ElectronEffPtBarrel1.SetDenominator(RecoPt);
  ElectronEffPtBarrel1.ComputeEfficiencies();
  ElectronEffPtBarrel1.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel1.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel1.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel1.SetYaxisMin(0.0);
  ElectronEffPtBarrel1.Write();

  sprintf(filename,"%s-EleMisidPtBarrel2.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel2(filename);
  ElectronEffPtBarrel2.AddNumerator(RecoPt);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtBarrel2.AddNumerator(BdtHWWIdPtBarrel2[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtBarrel2.AddNumerator(BdtHZZIdPtBarrel2[icut]);
  }

  ElectronEffPtBarrel2.SetDenominator(RecoPt);
  ElectronEffPtBarrel2.ComputeEfficiencies();
  ElectronEffPtBarrel2.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel2.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel2.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel2.SetYaxisMin(0.0);
  ElectronEffPtBarrel2.Write();

  sprintf(filename,"%s-EleMisidPtEndcap1.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap1(filename);
  ElectronEffPtEndcap1.AddNumerator(RecoPt);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtEndcap1.AddNumerator(BdtHWWIdPtEndcap1[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtEndcap1.AddNumerator(BdtHZZIdPtEndcap1[icut]);
  }

  ElectronEffPtEndcap1.SetDenominator(RecoPt);
  ElectronEffPtEndcap1.ComputeEfficiencies();
  ElectronEffPtEndcap1.SetTitle("fake rate vs pT");
  ElectronEffPtEndcap1.SetXaxisTitle("electron pT");
  ElectronEffPtEndcap1.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap1.SetYaxisMin(0.0);
  ElectronEffPtEndcap1.Write();

  sprintf(filename,"%s-EleMisidPtEndcap2.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap2(filename);
  ElectronEffPtEndcap2.AddNumerator(RecoPt);
  for (int icut=0;icut<(int)EgammaBdtHWWBasedID.size();++icut){
    ElectronEffPtEndcap2.AddNumerator(BdtHWWIdPtEndcap2[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPtEndcap2.AddNumerator(BdtHZZIdPtEndcap2[icut]);
  }

  ElectronEffPtEndcap2.SetDenominator(RecoPt);
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
    ElectronEffPU.AddNumerator(BdtHWWIdPU[icut]);
  }
  for (int icut=0;icut<(int)EgammaBdtHZZBasedID.size();++icut){
    ElectronEffPU.AddNumerator(BdtHZZIdPU[icut]);
  }

  ElectronEffPU.SetDenominator(RecoPU);
  ElectronEffPU.ComputeEfficiencies();
  ElectronEffPU.SetTitle("fake rate vs PU");
  ElectronEffPU.SetXaxisTitle("# vertices");
  ElectronEffPU.SetYaxisTitle("Fake rate");
  ElectronEffPU.SetYaxisMin(0.0);
  ElectronEffPU.Write();

  
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
    if(fabs(eta) >=  2.4) combIso -= 0.90 * rho;

    if(pt>20) {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.23); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.20);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.12);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.11);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.049);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.070);
      if(fabs(eta) >=  2.4) return (combIso/pt < 0.010);
    } else {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.20); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.21);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.13);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.10);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.04);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.03);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.03);
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
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.077);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.077);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.092);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.065);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.074);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.068);
  }

  if(type == kBDTHZZ_noIP) {
    // not optimized WP!
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.077);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.077);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.092);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.068);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.071);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.071);
  }
  return false;
}
