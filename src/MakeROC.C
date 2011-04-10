// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // MC efficiencies on W+jets Spring 11
  float eff_CutId_HighPt[5] = { 0.95, 0.88, 0.81, 0.76, 0.70 };
  float eff_LHId_HighPt[5] = { 0.94, 0.91, 0.83, 0.74, 0.67 };

  float eff_CutId_LowPt[5] = { 0.81, 0.68, 0.59, 0.54, 0.48 };
  float eff_LHId_LowPt[5] = { 0.84, 0.77, 0.65, 0.52, 0.44 };

  // MC fake rates on W+jets Spring 11
  float fr_CutId_HighPt[5] = { 0.128, 0.090, 0.061, 0.050, 0.043 };
  float fr_LHId_HighPt[5] = { 0.106, 0.087, 0.062, 0.045, 0.039 };

  float fr_CutId_LowPt[5] = { 0.093, 0.060, 0.041, 0.032, 0.028 };
  float fr_LHId_LowPt[5] = { 0.094, 0.070, 0.047, 0.030, 0.023 };

  TGraph *ROC_CutId_HighPt = new TGraph(5, eff_CutId_HighPt, fr_CutId_HighPt);
  TGraph *ROC_LHId_HighPt = new TGraph(5, eff_LHId_HighPt, fr_LHId_HighPt);

  TGraph *ROC_CutId_LowPt = new TGraph(5, eff_CutId_HighPt, fr_CutId_LowPt);
  TGraph *ROC_LHId_LowPt = new TGraph(5, eff_LHId_HighPt, fr_LHId_LowPt);

  // high pT ROC
  TCanvas *highPt = new TCanvas("highPt","highPt",600,400);
  ROC_CutId_HighPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_CutId_HighPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_CutId_HighPt->SetMarkerStyle(20);
  ROC_CutId_HighPt->SetMarkerSize(1.05);
  ROC_LHId_HighPt->SetMarkerStyle(21);
  ROC_LHId_HighPt->SetMarkerSize(1.05);
  ROC_CutId_HighPt->SetLineColor(1);
  ROC_CutId_HighPt->SetMarkerColor(1);
  ROC_LHId_HighPt->SetLineColor(2);
  ROC_LHId_HighPt->SetMarkerColor(2);

  TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(ROC_CutId_HighPt, "Cut Id","pl"); 
  leg->AddEntry(ROC_LHId_HighPt, "LH Id","pl"); 

  ROC_CutId_HighPt->Draw("apl");
  ROC_LHId_HighPt->Draw("pl");
  leg->Draw();
  highPt->SaveAs("roc_HighPt.eps");
  highPt->SaveAs("roc_HighPt.root");
  highPt->SaveAs("roc_HighPt.png");

  // low pT ROC
  TCanvas *lowPt = new TCanvas("lowPt","lowPt",600,400);
  ROC_CutId_LowPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_CutId_LowPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_CutId_LowPt->SetMarkerStyle(20);
  ROC_CutId_LowPt->SetMarkerSize(1.05);
  ROC_LHId_LowPt->SetMarkerStyle(21);
  ROC_LHId_LowPt->SetMarkerSize(1.05);
  ROC_CutId_LowPt->SetLineColor(1);
  ROC_CutId_LowPt->SetMarkerColor(1);
  ROC_LHId_LowPt->SetLineColor(2);
  ROC_LHId_LowPt->SetMarkerColor(2);

  ROC_CutId_LowPt->Draw("apl");
  ROC_LHId_LowPt->Draw("pl");
  leg->Draw();
  lowPt->SaveAs("roc_LowPt.eps");
  lowPt->SaveAs("roc_LowPt.root");
  lowPt->SaveAs("roc_LowPt.png");

}
