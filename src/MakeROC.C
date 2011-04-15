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

  // MC efficiencies on W+jets Spring 11 (EB)
  float eff_EB_CutId_HighPt[5] = { 0.97, 0.91, 0.85, 0.82, 0.77 };
  float eff_EB_LHId_HighPt[5] = { 0.96, 0.93, 0.87, 0.80, 0.73 };

  float eff_EB_CutId_LowPt[5] = { 0.87, 0.74, 0.65, 0.61, 0.55 };
  float eff_EB_LHId_LowPt[5] = { 0.90, 0.83, 0.73, 0.62, 0.50 };

  // MC efficiencies on W+jets Spring 11 (EE)
  float eff_EE_CutId_HighPt[5] = { 0.94, 0.83, 0.75, 0.68, 0.60 };
  float eff_EE_LHId_HighPt[5] = { 0.95, 0.90, 0.78, 0.68, 0.54 };

  float eff_EE_CutId_LowPt[5] = { 0.72, 0.59, 0.47, 0.39, 0.34 };
  float eff_EE_LHId_LowPt[5] = { 0.84, 0.73, 0.54, 0.37, 0.21 };

  // MC fake rates on W+jets Spring 11 (EB)
  float fr_EB_CutId_HighPt[5] = { 0.112, 0.084, 0.057, 0.050, 0.044 };
  float fr_EB_LHId_HighPt[5] = { 0.102, 0.075, 0.056, 0.046, 0.039 };

  float fr_EB_CutId_LowPt[5] = { 0.074, 0.046, 0.035, 0.030, 0.027 };
  float fr_EB_LHId_LowPt[5] = { 0.088, 0.055, 0.040, 0.030, 0.022 };

  // MC fake rates on W+jets Spring 11 (EE)
  float fr_EE_CutId_HighPt[5] = { 0.156, 0.096, 0.066, 0.046, 0.038 };
  float fr_EE_LHId_HighPt[5] = { 0.163, 0.119, 0.075, 0.047, 0.031 };

  float fr_EE_CutId_LowPt[5] = { 0.119, 0.078, 0.049, 0.033, 0.026 };
  float fr_EE_LHId_LowPt[5] = { 0.171, 0.109, 0.062, 0.033, 0.015 };

  

  TGraph *ROC_EB_CutId_HighPt = new TGraph(5, eff_EB_CutId_HighPt, fr_EB_CutId_HighPt);
  TGraph *ROC_EB_LHId_HighPt = new TGraph(5, eff_EB_LHId_HighPt, fr_EB_LHId_HighPt);

  TGraph *ROC_EE_CutId_HighPt = new TGraph(5, eff_EE_CutId_HighPt, fr_EE_CutId_HighPt);
  TGraph *ROC_EE_LHId_HighPt = new TGraph(5, eff_EE_LHId_HighPt, fr_EE_LHId_HighPt);

  TGraph *ROC_EB_CutId_LowPt = new TGraph(5, eff_EB_CutId_LowPt, fr_EB_CutId_LowPt);
  TGraph *ROC_EB_LHId_LowPt = new TGraph(5, eff_EB_LHId_LowPt, fr_EB_LHId_LowPt);

  TGraph *ROC_EE_CutId_LowPt = new TGraph(5, eff_EE_CutId_LowPt, fr_EE_CutId_LowPt);
  TGraph *ROC_EE_LHId_LowPt = new TGraph(5, eff_EE_LHId_LowPt, fr_EE_LHId_LowPt);


  // high pT ROC EB
  TCanvas *highPtEB = new TCanvas("highPtEB","highPtEB",600,400);
  ROC_EB_CutId_HighPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EB_CutId_HighPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EB_CutId_HighPt->SetMarkerStyle(20);
  ROC_EB_CutId_HighPt->SetMarkerSize(1.05);
  ROC_EB_LHId_HighPt->SetMarkerStyle(21);
  ROC_EB_LHId_HighPt->SetMarkerSize(1.05);
  ROC_EB_CutId_HighPt->SetLineColor(1);
  ROC_EB_CutId_HighPt->SetMarkerColor(1);
  ROC_EB_LHId_HighPt->SetLineColor(2);
  ROC_EB_LHId_HighPt->SetMarkerColor(2);

  TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(ROC_EB_CutId_HighPt, "Cut Id","pl"); 
  leg->AddEntry(ROC_EB_LHId_HighPt, "LH Id","pl"); 

  ROC_EB_CutId_HighPt->Draw("apl");
  ROC_EB_LHId_HighPt->Draw("pl");
  leg->Draw();
  highPtEB->SaveAs("roc_EB_HighPt.eps");
  highPtEB->SaveAs("roc_EB_HighPt.root");
  highPtEB->SaveAs("roc_EB_HighPt.png");

  // high pT ROC EE
  TCanvas *highPtEE = new TCanvas("highPtEE","highPtEE",600,400);
  ROC_EE_CutId_HighPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EE_CutId_HighPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EE_CutId_HighPt->SetMarkerStyle(20);
  ROC_EE_CutId_HighPt->SetMarkerSize(1.05);
  ROC_EE_LHId_HighPt->SetMarkerStyle(21);
  ROC_EE_LHId_HighPt->SetMarkerSize(1.05);
  ROC_EE_CutId_HighPt->SetLineColor(1);
  ROC_EE_CutId_HighPt->SetMarkerColor(1);
  ROC_EE_LHId_HighPt->SetLineColor(2);
  ROC_EE_LHId_HighPt->SetMarkerColor(2);

  ROC_EE_CutId_HighPt->Draw("apl");
  ROC_EE_LHId_HighPt->Draw("pl");
  leg->Draw();
  highPtEE->SaveAs("roc_EE_HighPt.eps");
  highPtEE->SaveAs("roc_EE_HighPt.root");
  highPtEE->SaveAs("roc_EE_HighPt.png");

  // low pT ROC EB
  TCanvas *lowPtEB = new TCanvas("lowPtEB","lowPtEB",600,400);
  ROC_EB_CutId_LowPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EB_CutId_LowPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EB_CutId_LowPt->SetMarkerStyle(20);
  ROC_EB_CutId_LowPt->SetMarkerSize(1.05);
  ROC_EB_LHId_LowPt->SetMarkerStyle(21);
  ROC_EB_LHId_LowPt->SetMarkerSize(1.05);
  ROC_EB_CutId_LowPt->SetLineColor(1);
  ROC_EB_CutId_LowPt->SetMarkerColor(1);
  ROC_EB_LHId_LowPt->SetLineColor(2);
  ROC_EB_LHId_LowPt->SetMarkerColor(2);

  ROC_EB_CutId_LowPt->Draw("apl");
  ROC_EB_LHId_LowPt->Draw("pl");
  leg->Draw();
  lowPtEB->SaveAs("roc_EB_LowPt.eps");
  lowPtEB->SaveAs("roc_EB_LowPt.root");
  lowPtEB->SaveAs("roc_EB_LowPt.png");

  // low pT ROC EE
  TCanvas *lowPtEE = new TCanvas("lowPtEE","lowPtEE",600,400);
  ROC_EE_CutId_LowPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EE_CutId_LowPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EE_CutId_LowPt->SetMarkerStyle(20);
  ROC_EE_CutId_LowPt->SetMarkerSize(1.05);
  ROC_EE_LHId_LowPt->SetMarkerStyle(21);
  ROC_EE_LHId_LowPt->SetMarkerSize(1.05);
  ROC_EE_CutId_LowPt->SetLineColor(1);
  ROC_EE_CutId_LowPt->SetMarkerColor(1);
  ROC_EE_LHId_LowPt->SetLineColor(2);
  ROC_EE_LHId_LowPt->SetMarkerColor(2);

  ROC_EE_CutId_LowPt->Draw("apl");
  ROC_EE_LHId_LowPt->Draw("pl");
  leg->Draw();
  lowPtEE->SaveAs("roc_EE_LowPt.eps");
  lowPtEE->SaveAs("roc_EE_LowPt.root");
  lowPtEE->SaveAs("roc_EE_LowPt.png");


}
