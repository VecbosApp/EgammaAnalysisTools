#include "TCanvas.h"
#include "TFile.h"
#include "RooHist.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TString.h"

void compare(const char* file1, const char* file2, const char* variable="vertices", const char* effPoint="80", bool writeEBEE=false);

void doAll() {

  // eta
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Eta_WP95.root","NewCatOpt_1oEm1oP/LowPt/effData_Eta_LHVeryLoose.root","eta","95",true);
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Eta_WP90.root","NewCatOpt_1oEm1oP/LowPt/effData_Eta_LHLoose.root","eta","90",true);
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Eta_WP85.root","NewCatOpt_1oEm1oP/LowPt/effData_Eta_LHMedium.root","eta","85",true);
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Eta_WP80.root","NewCatOpt_1oEm1oP/LowPt/effData_Eta_LHTight.root","eta","80",true);
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Eta_WP70.root","NewCatOpt_1oEm1oP/LowPt/effData_Eta_LHHyperTight.root","eta","70",true);

  // pT
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Pt_WP95.root","NewCatOpt_1oEm1oP/LowPt/effData_Pt_LHVeryLoose.root","pt","95");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Pt_WP90.root","NewCatOpt_1oEm1oP/LowPt/effData_Pt_LHLoose.root","pt","90");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Pt_WP85.root","NewCatOpt_1oEm1oP/LowPt/effData_Pt_LHMedium.root","pt","85");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Pt_WP80.root","NewCatOpt_1oEm1oP/LowPt/effData_Pt_LHTight.root","pt","80");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Pt_WP70.root","NewCatOpt_1oEm1oP/LowPt/effData_Pt_LHHyperTight.root","pt","70");

  // vertices
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Vertices_WP95.root","NewCatOpt_1oEm1oP/LowPt/effData_Vertices_LHVeryLoose.root","vertices","95");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Vertices_WP90.root","NewCatOpt_1oEm1oP/LowPt/effData_Vertices_LHLoose.root","vertices","90");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Vertices_WP85.root","NewCatOpt_1oEm1oP/LowPt/effData_Vertices_LHMedium.root","vertices","85");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Vertices_WP80.root","NewCatOpt_1oEm1oP/LowPt/effData_Vertices_LHTight.root","vertices","80");
  compare("NewCatOpt_1oEm1oP/LowPt/effData_Vertices_WP70.root","NewCatOpt_1oEm1oP/LowPt/effData_Vertices_LHHyperTight.root","vertices","70");

}

void compare(const char* file1, const char* file2, const char* variable, const char* effPoint, bool writeEBEE) {

  gROOT->SetStyle("Plain");

  TString dir = TString("eleIDdir/")+TString(variable)+TString("/fit_eff_plots/")+TString(variable)+TString("_PLOT");

  TFile *TFile1 = TFile::Open(file1);
  TCanvas *c1 = (TCanvas*)TFile1->Get(dir.Data());
  RooHist* h1 = (RooHist*)c1->GetPrimitive("hxy_fit_eff");
  
  TFile *TFile2 = TFile::Open(file2);
  TCanvas *c2 = (TCanvas*)TFile2->Get(dir.Data());
  RooHist* h2 = (RooHist*)c2->GetPrimitive("hxy_fit_eff");

  TCanvas *c3 = new TCanvas("c3","c3",600,400);
  h1->GetXaxis()->SetTitle(variable);
  h2->SetMarkerColor(kRed);
  h2->SetLineColor(kRed);
  h1->SetMinimum(0.0);
  h1->SetMaximum(1.0);
  h1->Draw("ap");
  h2->Draw("p");

  // compute the weighted average
  Double_t* x1 = h1->GetX();
  Double_t* y1 = h1->GetY();
  Double_t* y2 = h2->GetY();
  Double_t* ey1 = h1->GetEYlow();
  Double_t* ey2 = h2->GetEYlow();
  Int_t nbins = h1->GetN();

  Double_t num1, denom1, num2, denom2;
  num1 = denom1 = num2 = denom2 = 0.0;

  Double_t num1_eta[2], denom1_eta[2], num2_eta[2], denom2_eta[2];
  num1_eta[0] = denom1_eta[0] = num2_eta[0] = denom2_eta[0] = 0.0;
  num1_eta[1] = denom1_eta[1] = num2_eta[1] = denom2_eta[1] = 0.0;

  for(int i=0; i<nbins; i++) {
    std::cout << "i = " << i << " WP eff = " << ey1[i] << " LH eff = " << ey2[i] << std::endl;
    if(fabs(ey1[i])>0) {
      num1 += y1[i]/ey1[i]/ey1[i];
      denom1 += 1/ey1[i]/ey1[i];

      if(writeEBEE) {
	// now compute the EB/EE separately
	int subdet;
	if(fabs(x1[i])<1.479) subdet = 0;
	else subdet = 1;
	num1_eta[subdet] += y1[i]/ey1[i]/ey1[i];
	denom1_eta[subdet] += 1/ey1[i]/ey1[i];
      }
    }
    if(fabs(ey2[i])>0) {
      num2 += y2[i]/ey2[i]/ey2[i];
      denom2 += 1/ey2[i]/ey2[i];

      if(writeEBEE) {
	// now compute the EB/EE separately
	int subdet;
	if(fabs(x1[i])<1.479) subdet = 0;
	else subdet = 1;
	num2_eta[subdet] += y2[i]/ey2[i]/ey2[i];
	denom2_eta[subdet] += 1/ey2[i]/ey2[i];
      }
    }
  }
  std::cout << "EFFICIENCY FOR " << effPoint << " (WP) = " << num1/denom1 
	    << " (LH) = " << num2/denom2 << std::endl;

  char CutIdEff[50], LHIdEff[50];
  char CutIdEff_EB[50], LHIdEff_EB[50], CutIdEff_EE[50], LHIdEff_EE[50];
  
  sprintf(CutIdEff,"Cut Id WP%s (#epsilon = %.2f)", effPoint, num1/denom1);
  sprintf(LHIdEff,"LH Id WP%s (#epsilon = %.2f)", effPoint, num2/denom2);
  
  if(writeEBEE) {
    sprintf(CutIdEff_EB,"Cut Id WP%s (#epsilon EB = %.2f)", effPoint, num1_eta[0]/denom1_eta[0]);
    sprintf(CutIdEff_EE,"Cut Id WP%s (#epsilon EE = %.2f)", effPoint, num1_eta[1]/denom1_eta[1]);

    sprintf(LHIdEff_EB,"LH Id WP%s (#epsilon EB = %.2f)", effPoint, num2_eta[0]/denom2_eta[0]);
    sprintf(LHIdEff_EE,"LH Id WP%s (#epsilon EE = %.2f)", effPoint, num2_eta[1]/denom2_eta[1]);
  }

  TLegend* leg = new TLegend(0.15, 0.20, 0.85, 0.4);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.025);
  leg->SetFillColor(0);
  if(writeEBEE) {
    leg->AddEntry(h1,CutIdEff_EB,"p");
    leg->AddEntry(h1,CutIdEff_EE,"p");
    leg->AddEntry(h2,LHIdEff_EB,"p");
    leg->AddEntry(h2,LHIdEff_EE,"p");
  } else {
    leg->AddEntry(h1,CutIdEff,"p");
    leg->AddEntry(h2,LHIdEff,"p");
  }
  leg->Draw();

  TString fileOut = TString("eff_vs_")+TString(variable)+TString("_")+TString(effPoint)+TString("_lowPt.png");
  c3->SaveAs(fileOut.Data());

  fileOut = TString("eff_vs_")+TString(variable)+TString("_")+TString(effPoint)+TString("_lowPt.eps");
  c3->SaveAs(fileOut.Data());

}
