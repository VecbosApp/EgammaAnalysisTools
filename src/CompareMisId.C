// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>

// Offline analysis includes
#include "CommonTools/include/EfficiencyEvaluator.hh"

using namespace std;

void drawEta(const char* filename);
void drawPt(const char* filename);

int main(int argc, char* argv[]) {

  char inputFileNameEta[150];
  char inputFileNamePt[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: compareMisid inputFileEta.root inputFilePt.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameEta,argv[1]);
  strcpy(inputFileNamePt,argv[2]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  drawEta(inputFileNameEta);
  drawPt(inputFileNamePt);


}


void drawEta(const char* filename) {

  TFile *efficiencyFileEta = TFile::Open(filename);

  TH1F *FakeableJetsEta = (TH1F*)efficiencyFileEta->Get("FakeableJetsEta");
  TH1F *RecoEta = (TH1F*)efficiencyFileEta->Get("RecoEta");
  TH1F *CutIdEta = (TH1F*)efficiencyFileEta->Get("CutIdEta");
  TH1F *LHIdLooseEta = (TH1F*)efficiencyFileEta->Get("LHIdLooseEta");
  TH1F *LHIdTightEta = (TH1F*)efficiencyFileEta->Get("LHIdTightEta");

  EfficiencyEvaluator ElectronFakeRateEta("elemisid-eta.root");
  ElectronFakeRateEta.AddNumerator(FakeableJetsEta);
  ElectronFakeRateEta.AddNumerator(RecoEta);
  ElectronFakeRateEta.AddNumerator(CutIdEta);
  ElectronFakeRateEta.AddNumerator(LHIdLooseEta);
  ElectronFakeRateEta.AddNumerator(LHIdTightEta);
  ElectronFakeRateEta.SetDenominator(FakeableJetsEta);
  ElectronFakeRateEta.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronFakeRateEta.GetCumulativeEfficiencies();

  TCanvas c1;
  c1.SetLogy();
  efficiency[1]->SetLineColor(2);
  efficiency[2]->SetLineColor(4);
  efficiency[3]->SetLineColor(6);
  efficiency[4]->SetLineColor(7);

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(efficiency[1],"reconstrunction");
  leg->AddEntry(efficiency[2],"Loose Category standard");
  leg->AddEntry(efficiency[3],"Loose Likelihood");
  leg->AddEntry(efficiency[4],"Tight Likelihood");

  efficiency[1]->SetMinimum(0.0001);
  efficiency[1]->SetMaximum(1.0);
  efficiency[1]->SetTitle("");
  efficiency[1]->GetXaxis()->SetTitle("#eta of closest jet");
  efficiency[1]->GetYaxis()->SetTitle("Electron mis-id from Jets");
  efficiency[1]->Draw("hist");
  for(int i=2; i<5; i++) {
    efficiency[i]->Draw("same hist");
  }
  leg->Draw();

  c1.SaveAs("elemisid-eta.eps");
  c1.SaveAs("elemisid-eta.root");

}


void drawPt(const char* filename) {

  TFile *efficiencyFilePt = TFile::Open(filename);

  TH1F *FakeableJetsPt = (TH1F*)efficiencyFilePt->Get("FakeableJetsPt");
  TH1F *RecoPt = (TH1F*)efficiencyFilePt->Get("RecoPt");
  TH1F *CutIdPt = (TH1F*)efficiencyFilePt->Get("CutIdPt");
  TH1F *LHIdLoosePt = (TH1F*)efficiencyFilePt->Get("LHIdLoosePt");
  TH1F *LHIdTightPt = (TH1F*)efficiencyFilePt->Get("LHIdTightPt");

  EfficiencyEvaluator ElectronFakeRatePt("elemisid-pt.root");
  ElectronFakeRatePt.AddNumerator(FakeableJetsPt);
  ElectronFakeRatePt.AddNumerator(RecoPt);
  ElectronFakeRatePt.AddNumerator(CutIdPt);
  ElectronFakeRatePt.AddNumerator(LHIdLoosePt);
  ElectronFakeRatePt.AddNumerator(LHIdTightPt);
  ElectronFakeRatePt.SetDenominator(FakeableJetsPt);
  ElectronFakeRatePt.ComputeEfficiencies();

  vector<TH1F*> efficiency = ElectronFakeRatePt.GetCumulativeEfficiencies();

  TCanvas c1;
  c1.SetLogy();
  efficiency[1]->SetLineColor(2);
  efficiency[2]->SetLineColor(4);
  efficiency[3]->SetLineColor(6);
  efficiency[4]->SetLineColor(7);

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(efficiency[1],"reconstrunction");
  leg->AddEntry(efficiency[2],"Loose Category standard");
  leg->AddEntry(efficiency[3],"Loose Likelihood");
  leg->AddEntry(efficiency[4],"Tight Likelihood");

  efficiency[1]->SetMinimum(0.0001);
  efficiency[1]->SetMaximum(1.0);
  efficiency[1]->SetTitle("");
  efficiency[1]->GetXaxis()->SetTitle("p_{T} of closest jet [GeV]");
  efficiency[1]->GetYaxis()->SetTitle("Electron mis-id from Jets");
  efficiency[1]->Draw("hist");
  for(int i=2; i<5; i++) {
    efficiency[i]->Draw("same hist");
  }
  leg->Draw();

  c1.SaveAs("elemisid-pt.eps");
  c1.SaveAs("elemisid-pt.root");

}
