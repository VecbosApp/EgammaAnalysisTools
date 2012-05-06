#define estimateMuonFakeRateHzz4lTree_cxx
#include "estimateMuonFakeRateHzz4lTree.h"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void estimateMuonFakeRateHzz4lTree::Loop(const char *outname)
{
  if (fChain == 0) return;

  std::vector<TString> NoTrgMuonID;
  NoTrgMuonID.push_back("hzzPfIso");  // reference pf iso
  NoTrgMuonID.push_back("hzzMvaPfIso");  // duncan's optimization for HZZ

  // -----------------------------------------------------------------------
  // study vs eta

  Float_t LowerEta[5];
  LowerEta[0]=0.0;
  LowerEta[1]=1.0;
  LowerEta[2]=1.479;  
  LowerEta[3]=2.0;
  LowerEta[4]=2.5;
  TH1F *RecoEtaHighPt   = new TH1F( "RecoEtaHighPt",    "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaLowPt    = new TH1F( "RecoEtaLowPt",     "reconstructed #eta", 4, LowerEta);
  
  // eta, high pT
  std::vector<TH1F*> NoTrgMuonEtaHighPt;
  for (int i=0;i<(int)NoTrgMuonID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgMuon"+TString(NoTrgMuonID[i])+"EtaHighPt",   "HZZ BDT ID #eta", 4, LowerEta);     NoTrgMuonEtaHighPt.push_back(aHisto);
  }

  // eta, low pT
  std::vector<TH1F*> NoTrgMuonEtaLowPt;
  for (int i=0;i<(int)NoTrgMuonID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgMuon"+TString(NoTrgMuonID[i])+"EtaLowPt",   "HZZ BDT ID #eta", 4, LowerEta);
    NoTrgMuonEtaLowPt.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs pT
  Float_t LowerPt[11] = {0.0,5.0,7.0,10.0,14.0,20.0,25.0,30.0,40.0,50.0,80.0};
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", 10, LowerPt);

  // to have the full picture in the barrel
  std::vector<TH1F*> NoTrgMuonPtBarrel;
  for (int i=0;i<(int)NoTrgMuonID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgMuon"+TString(NoTrgMuonID[i])+"PtBarrel", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgMuonPtBarrel.push_back(aHisto);
  }

  // to have the full picture in the endcap
  std::vector<TH1F*> NoTrgMuonPtEndcap;
  for (int i=0;i<(int)NoTrgMuonID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgMuon"+TString(NoTrgMuonID[i])+"PtEndcap", "HZZ BDT ID #eta",   10, LowerPt);
    NoTrgMuonPtEndcap.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs PU
  Float_t LowerPU[11];
  LowerPU[0] = 1;
  LowerPU[1] = 4;
  for(int i=2;i<6;++i) LowerPU[i]=i+3;
  LowerPU[6] = 10;
  LowerPU[7] = 15;
  LowerPU[8] = 25;
  LowerPU[9] = 30;
  LowerPU[10] = 35;

  TH1F *RecoPUBarrel   = new TH1F( "RecoPUBarrel",   "reconstructed nPU", 10, LowerPU);
  TH1F *RecoPUEndcap   = new TH1F( "RecoPUEndcap",   "reconstructed nPU", 10, LowerPU);

  //  barrel
  std::vector<TH1F*> NoTrgMuonPUBarrel;
  for (int i=0;i<(int)NoTrgMuonID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgMuon"+TString(NoTrgMuonID[i])+"PUBarrel", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgMuonPUBarrel.push_back(aHisto);
  }

  // endcap
  std::vector<TH1F*> NoTrgMuonPUEndcap;
  for (int i=0;i<(int)NoTrgMuonID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgMuon"+TString(NoTrgMuonID[i])+"PUEndcap", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgMuonPUEndcap.push_back(aHisto);
  }

  // loop on events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // fill the denominator: take only the highest pT denominator candidate
    float etaFake = fabs(eta);
    float etFake  = pt;
    bool isInEB   = fabs(eta)<1.479;
    bool isInEE   = !isInEB;
    bool highPt   = (pt>10.);
    bool lowPt    = (pt<=10.);

    // pass the denominator object
    if(etFake<5 || etaFake>2.4 || fabs(sip3d)>100) continue;

    // filling
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


    // fill the numerator(s)
    // === PF-muon ID, combined Iso ===
    if(passRefMuSel()) {
      if (highPt) NoTrgMuonEtaHighPt[khzzPfIso]->Fill(etaFake);
      if (lowPt)  NoTrgMuonEtaLowPt[khzzPfIso] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgMuonPtBarrel[khzzPfIso] ->Fill(etFake);
	NoTrgMuonPUBarrel[khzzPfIso] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgMuonPtEndcap[khzzPfIso] ->Fill(etFake);
	NoTrgMuonPUEndcap[khzzPfIso] ->Fill(vertices);
      }
    }

    // === PF-muon ID, MVA Iso ===
    if(passMvaMuSel(false)) {
      if (highPt) NoTrgMuonEtaHighPt[khzzMvaPfIso]->Fill(etaFake);
      if (lowPt)  NoTrgMuonEtaLowPt[khzzMvaPfIso] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgMuonPtBarrel[khzzMvaPfIso] ->Fill(etFake);
	NoTrgMuonPUBarrel[khzzMvaPfIso] ->Fill(vertices);
      }
      if (isInEE) {
	NoTrgMuonPtEndcap[khzzMvaPfIso] ->Fill(etFake);
	NoTrgMuonPUEndcap[khzzMvaPfIso] ->Fill(vertices);
      }
    }

  } // loop over events

    // saving efficiency histos
  // === as a function of eta ===
  char filename[500];
  sprintf(filename,"%s-MuonMisidEtaHighPt.root",outname);
  EfficiencyEvaluator MuonEffEtaHighPt(filename);
  MuonEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<(int)NoTrgMuonID.size();++icut){
    MuonEffEtaHighPt.AddNumerator(NoTrgMuonEtaHighPt[icut]);
  }

  MuonEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  MuonEffEtaHighPt.ComputeEfficiencies();
  MuonEffEtaHighPt.SetTitle("fake rate vs pT");
  MuonEffEtaHighPt.SetXaxisTitle("electron p_{T}");
  MuonEffEtaHighPt.SetYaxisTitle("Fake rate");
  MuonEffEtaHighPt.SetYaxisMin(0.0);
  MuonEffEtaHighPt.Write();

  sprintf(filename,"%s-MuonMisidEtaLowPt.root",outname);
  EfficiencyEvaluator MuonEffEtaLowPt(filename);
  MuonEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<(int)NoTrgMuonID.size();++icut){
    MuonEffEtaLowPt.AddNumerator(NoTrgMuonEtaLowPt[icut]);
  }

  MuonEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  MuonEffEtaLowPt.ComputeEfficiencies();
  MuonEffEtaLowPt.SetTitle("fake rate vs eta");
  MuonEffEtaLowPt.SetXaxisTitle("electron #eta");
  MuonEffEtaLowPt.SetYaxisTitle("Fake rate");
  MuonEffEtaLowPt.SetYaxisMin(0.0);
  MuonEffEtaLowPt.Write();

  // === as a function of pt ===
  sprintf(filename,"%s-MuonMisidPtBarrel.root",outname);
  EfficiencyEvaluator MuonEffPtBarrel(filename);
  MuonEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<(int)NoTrgMuonID.size();++icut){
    MuonEffPtBarrel.AddNumerator(NoTrgMuonPtBarrel[icut]);
  }

  MuonEffPtBarrel.SetDenominator(RecoPtBarrel);
  MuonEffPtBarrel.ComputeEfficiencies();
  MuonEffPtBarrel.SetTitle("fake rate vs pT");
  MuonEffPtBarrel.SetXaxisTitle("electron pT");
  MuonEffPtBarrel.SetYaxisTitle("Fake rate");
  MuonEffPtBarrel.SetYaxisMin(0.0);
  MuonEffPtBarrel.Write();

  sprintf(filename,"%s-MuonMisidPtEndcap.root",outname);
  EfficiencyEvaluator MuonEffPtEndcap(filename);
  MuonEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<(int)NoTrgMuonID.size();++icut){
    MuonEffPtEndcap.AddNumerator(NoTrgMuonPtEndcap[icut]);
  }

  MuonEffPtEndcap.SetDenominator(RecoPtEndcap);
  MuonEffPtEndcap.ComputeEfficiencies();
  MuonEffPtEndcap.SetTitle("fake rate vs pT");
  MuonEffPtEndcap.SetXaxisTitle("electron pT");
  MuonEffPtEndcap.SetYaxisTitle("Fake rate");
  MuonEffPtEndcap.SetYaxisMin(0.0);
  MuonEffPtEndcap.Write();

  // === as a function of PU ===
  sprintf(filename,"%s-MuonMisidPUBarrel.root",outname);
  EfficiencyEvaluator MuonEffPUBarrel(filename);
  MuonEffPUBarrel.AddNumerator(RecoPUBarrel);
  for (int icut=0;icut<(int)NoTrgMuonID.size();++icut){
    MuonEffPUBarrel.AddNumerator(NoTrgMuonPUBarrel[icut]);
  }

  MuonEffPUBarrel.SetDenominator(RecoPUBarrel);
  MuonEffPUBarrel.ComputeEfficiencies();
  MuonEffPUBarrel.SetTitle("fake rate vs vertices");
  MuonEffPUBarrel.SetXaxisTitle("# vertices");
  MuonEffPUBarrel.SetYaxisTitle("Fake rate");
  MuonEffPUBarrel.SetYaxisMin(0.0);
  MuonEffPUBarrel.Write();

  sprintf(filename,"%s-MuonMisidPUEndcap.root",outname);
  EfficiencyEvaluator MuonEffPUEndcap(filename);
  MuonEffPUEndcap.AddNumerator(RecoPUEndcap);
  for (int icut=0;icut<(int)NoTrgMuonID.size();++icut){
    MuonEffPUEndcap.AddNumerator(NoTrgMuonPUEndcap[icut]);
  }

  MuonEffPUEndcap.SetDenominator(RecoPUEndcap);
  MuonEffPUEndcap.ComputeEfficiencies();
  MuonEffPUEndcap.SetTitle("fake rate vs vertices");
  MuonEffPUEndcap.SetXaxisTitle("# vertices");
  MuonEffPUEndcap.SetYaxisTitle("Fake rate");
  MuonEffPUEndcap.SetYaxisMin(0.0);
  MuonEffPUEndcap.Write();

}
