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
    std::cout << "usage: compareEff inputFileEta.root inputFilePt.root" << std::endl;
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

  TH1F *GenEta = (TH1F*)efficiencyFileEta->Get("GenEta");
  TH1F *RecoEta = (TH1F*)efficiencyFileEta->Get("RecoEta");
  TH1F *CutIdWP80Eta = (TH1F*)efficiencyFileEta->Get("CutIdWP80Eta");
  TH1F *CutIdTightCICEta = (TH1F*)efficiencyFileEta->Get("CutIdTightCICEta");
  TH1F *CutIdSuperTightCICEta = (TH1F*)efficiencyFileEta->Get("CutIdSuperTightCICEta");
  TH1F *LHIdLooseEta = (TH1F*)efficiencyFileEta->Get("LHIdLooseEta");
  TH1F *LHIdTightEta = (TH1F*)efficiencyFileEta->Get("LHIdTightEta");

  EfficiencyEvaluator ElectronEfficiencyEta("eleeff-eta.root");
  ElectronEfficiencyEta.AddNumerator(GenEta);
  ElectronEfficiencyEta.AddNumerator(RecoEta);
  ElectronEfficiencyEta.AddNumerator(CutIdWP80Eta);
  ElectronEfficiencyEta.AddNumerator(CutIdTightCICEta);
  ElectronEfficiencyEta.AddNumerator(CutIdSuperTightCICEta);
  ElectronEfficiencyEta.AddNumerator(LHIdLooseEta);
  ElectronEfficiencyEta.AddNumerator(LHIdTightEta);
  ElectronEfficiencyEta.SetDenominator(GenEta);
  ElectronEfficiencyEta.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronEfficiencyEta.GetCumulativeEfficiencies();

  TCanvas c1;
  efficiency[1]->SetLineColor(2);
  efficiency[2]->SetLineColor(4);
  efficiency[3]->SetLineColor(6);
  efficiency[6]->SetLineColor(8);

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(efficiency[1],"reconstrunction");
  leg->AddEntry(efficiency[2],"WP80 Cuts");
  //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
  leg->AddEntry(efficiency[6],"Tight Likelihood");

  efficiency[1]->SetMinimum(0.0);
  efficiency[1]->SetMaximum(1.0);
  efficiency[1]->SetTitle("");
  efficiency[1]->GetXaxis()->SetTitle("electron #eta");
  efficiency[1]->GetYaxis()->SetTitle("Efficiency");
  efficiency[1]->Draw("hist");
  efficiency[2]->Draw("same hist");
  //  efficiency[3]->Draw("same hist");
  efficiency[6]->Draw("same hist");
  leg->Draw();

  c1.SaveAs("eleeff-eta.eps");
  c1.SaveAs("eleeff-eta.root");

}


void drawPt(const char* filename) {

  TFile *efficiencyFilePt = TFile::Open(filename);
 
  TH1F *GenPt = (TH1F*)efficiencyFilePt->Get("GenPt");
  TH1F *RecoPt = (TH1F*)efficiencyFilePt->Get("RecoPt");
  TH1F *CutIdWP80Pt = (TH1F*)efficiencyFilePt->Get("CutIdWP80Pt");
  TH1F *CutIdTightCICPt = (TH1F*)efficiencyFilePt->Get("CutIdTightCICPt");
  TH1F *CutIdSuperTightCICPt = (TH1F*)efficiencyFilePt->Get("CutIdSuperTightCICPt");
  TH1F *LHIdLoosePt = (TH1F*)efficiencyFilePt->Get("LHIdLoosePt");
  TH1F *LHIdTightPt = (TH1F*)efficiencyFilePt->Get("LHIdTightPt");

  EfficiencyEvaluator ElectronEfficiencyPt("eleeff-pt.root");
  ElectronEfficiencyPt.AddNumerator(GenPt);
  ElectronEfficiencyPt.AddNumerator(RecoPt);
  ElectronEfficiencyPt.AddNumerator(CutIdWP80Pt);
  ElectronEfficiencyPt.AddNumerator(CutIdTightCICPt);
  ElectronEfficiencyPt.AddNumerator(CutIdSuperTightCICPt);
  ElectronEfficiencyPt.AddNumerator(LHIdLoosePt);
  ElectronEfficiencyPt.AddNumerator(LHIdTightPt);
  ElectronEfficiencyPt.SetDenominator(GenPt);
  ElectronEfficiencyPt.ComputeEfficiencies();

  vector<TH1F*> efficiency = ElectronEfficiencyPt.GetCumulativeEfficiencies();

  TCanvas c1;
  efficiency[1]->SetLineColor(2);
  efficiency[2]->SetLineColor(4);
  //  efficiency[3]->SetLineColor(6);
  efficiency[6]->SetLineColor(8);

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(efficiency[1],"reconstrunction");
  leg->AddEntry(efficiency[2],"WP80 Cuts");
  //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
  leg->AddEntry(efficiency[6],"Tight Likelihood");

  efficiency[1]->SetMinimum(0.0);
  efficiency[1]->SetMaximum(1.0);
  efficiency[1]->SetTitle("");
  efficiency[1]->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  efficiency[1]->GetYaxis()->SetTitle("Efficiency");
  efficiency[1]->Draw("hist");
  efficiency[2]->Draw("same hist");
  //  efficiency[3]->Draw("same hist");
  efficiency[6]->Draw("same hist");
  leg->Draw();

  c1.SaveAs("eleeff-pt.eps");
  c1.SaveAs("eleeff-pt.root");

}
