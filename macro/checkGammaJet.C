#define checkGammaJet_cxx
#include "checkGammaJet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void checkGammaJet::Loop()
{
  TH1F *H_isoTracker0 = new TH1F("H_isoTracker0","H_isoTracker0",50, 0.,20.);
  TH1F *H_isoTracker1 = new TH1F("H_isoTracker1","H_isoTracker1",50, 0.,20.);
  TH1F *H_isoEcal0    = new TH1F("H_isoEcal0",   "H_isoEcal0",   50, 0.,20.);
  TH1F *H_isoEcal1    = new TH1F("H_isoEcal1",   "H_isoEcal1",   50, 0.,20.);
  TH1F *H_isoHcal0    = new TH1F("H_isoHcal0",   "H_isoHcal0",   50, 0.,20.);
  TH1F *H_isoHcal1    = new TH1F("H_isoHcal1",   "H_isoHcal1",   50, 0.,20.);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // histos for probes matching MC gamma
    if (isGamma==1){ 
      H_isoTracker1 -> Fill(absTrackerIsolGammaCand, weight);
      H_isoEcal1    -> Fill(absEcalIsolGammaCand,    weight);
      H_isoHcal1    -> Fill(absHcalIsolGammaCand,    weight);
    }
    
    // histos for probes not matching MC gamma
    if (isGamma==0){ 
      H_isoTracker0 -> Fill(absTrackerIsolGammaCand, weight);
      H_isoEcal0    -> Fill(absEcalIsolGammaCand,    weight);
      H_isoHcal0    -> Fill(absHcalIsolGammaCand,    weight);
    }
  }

  // rescaling
  H_isoTracker1->Scale(1./H_isoTracker1->Integral());
  H_isoEcal1   ->Scale(1./H_isoEcal1->Integral());
  H_isoHcal1   ->Scale(1./H_isoHcal1->Integral());
  H_isoTracker0->Scale(1./H_isoTracker0->Integral());
  H_isoEcal0   ->Scale(1./H_isoEcal0->Integral());
  H_isoHcal0   ->Scale(1./H_isoHcal0->Integral());

  // cosmetics
  H_isoTracker1->SetLineColor(2);
  H_isoEcal1   ->SetLineColor(2);
  H_isoHcal1   ->SetLineColor(2);
  H_isoTracker0->SetLineColor(1);
  H_isoEcal0   ->SetLineColor(1);
  H_isoHcal0   ->SetLineColor(1);

  // saving
  TFile file("file.root","RECREATE");
  H_isoTracker1->Write();
  H_isoEcal1   ->Write();
  H_isoHcal1   ->Write();
  H_isoTracker0->Write();
  H_isoEcal0   ->Write();
  H_isoHcal0   ->Write();  
}
