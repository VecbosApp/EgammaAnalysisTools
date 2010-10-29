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

#include <string>

using namespace std;

std::vector<std::string> EgammaCutBasedIDWPs;
std::vector<std::string> EgammaCiCBasedIDWPs;
std::vector<std::string> EgammaLHBasedIDWPs;

void drawEta(const char* filename);
void drawPt(const char* filename);




int main(int argc, char* argv[]) {

  EgammaCutBasedIDWPs.push_back("WP95");
  EgammaCutBasedIDWPs.push_back("WP90");
  EgammaCutBasedIDWPs.push_back("WP80");
  //  EgammaCutBasedIDWPs.push_back("WP70");

  //  EgammaCiCBasedIDWPs.push_back("CiCVeryLoose");
  EgammaCiCBasedIDWPs.push_back("CiCLoose");
  EgammaCiCBasedIDWPs.push_back("CiCMedium");
  EgammaCiCBasedIDWPs.push_back("CiCSuperTight");
  //  EgammaCiCBasedIDWPs.push_back("CiCHyperTight");

  EgammaLHBasedIDWPs.push_back("LHVeryLoose");
  EgammaLHBasedIDWPs.push_back("LHLoose");
  EgammaLHBasedIDWPs.push_back("LHMedium");
  //EgammaLHBasedIDWPs.push_back("LHTight");
  //  EgammaLHBasedIDWPs.push_back("LHHyperTight");

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

  //  TH1F *GenEta = (TH1F*)efficiencyFileEta->Get("GenEta");
  TH1F *RecoEta = (TH1F*)efficiencyFileEta->Get("RecoEta");

  std::vector<TH1F*> CutIdEta;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta");
      CutIdEta.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdEta;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta");
      LHIdEta.push_back(aHisto);
    }
  std::vector<TH1F*> CiCIdEta;
  for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta");
      CiCIdEta.push_back(aHisto);
    }




  EfficiencyEvaluator ElectronEfficiencyEta("eleeff-eta.root");
  //  ElectronEfficiencyEta.AddNumerator(GenEta);
  //  ElectronEfficiencyEta.AddNumerator(RecoEta);
  for (int i=0;i<CutIdEta.size();++i)
    ElectronEfficiencyEta.AddNumerator(CutIdEta[i]);
  for (int i=0;i<LHIdEta.size();++i)
    ElectronEfficiencyEta.AddNumerator(LHIdEta[i]);
  for (int i=0;i<CiCIdEta.size();++i)
    ElectronEfficiencyEta.AddNumerator(CiCIdEta[i]);
  ElectronEfficiencyEta.SetDenominator(RecoEta);
  ElectronEfficiencyEta.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronEfficiencyEta.GetCumulativeEfficiencies();

  for (int tightness=0;tightness<CutIdEta.size();++tightness)
    { 
      TCanvas c1;
      efficiency[0+tightness]->SetLineColor(1);
      efficiency[0+tightness]->SetMarkerColor(1);
      efficiency[CutIdEta.size()+tightness]->SetLineColor(2);
      efficiency[CutIdEta.size()+tightness]->SetMarkerColor(2);
      efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetLineColor(4);
      efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetMarkerColor(4);
      //  efficiency[6]->SetLineColor(8);
      
      TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
      leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
      leg->SetFillColor(0);
      leg->AddEntry(efficiency[0+tightness],EgammaCutBasedIDWPs[0+tightness].c_str());
      leg->AddEntry(efficiency[CutIdEta.size()+tightness],EgammaLHBasedIDWPs[0+tightness].c_str());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      leg->AddEntry(efficiency[CutIdEta.size()+LHIdEta.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].c_str());
      
      efficiency[0+tightness]->SetMinimum(0.3);
      efficiency[0+tightness]->SetMaximum(1.0);
      efficiency[0+tightness]->SetMarkerStyle(20);
      efficiency[0+tightness]->SetMarkerSize(1.05);
      efficiency[0+tightness]->SetTitle("");
      efficiency[0+tightness]->GetXaxis()->SetTitle("electron #eta");
      efficiency[0+tightness]->GetYaxis()->SetTitle("Efficiency");
      efficiency[0+tightness]->Draw("PEhist");
      efficiency[CutIdEta.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdEta.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdEta.size()+tightness]->SetMarkerSize(1.05);
      //  efficiency[3+tightness]->Draw("same hist");
      efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetMarkerSize(1.05);
      leg->Draw();
      
      char tightLevel[10];
      sprintf(tightLevel,"%d",tightness);
      c1.SaveAs("eleeff-eta-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-eta-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-eta-tightness"+TString(tightLevel)+".png");
    }
}


void drawPt(const char* filename) {
  TFile *efficiencyFilePt = TFile::Open(filename);

  //  TH1F *GenPt = (TH1F*)efficiencyFilePt->Get("GenPt");
  TH1F *RecoPt = (TH1F*)efficiencyFilePt->Get("RecoPt");

  std::vector<TH1F*> CutIdPt;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt");
      CutIdPt.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdPt;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt");
      LHIdPt.push_back(aHisto);
    }
  std::vector<TH1F*> CiCIdPt;
  for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt");
      CiCIdPt.push_back(aHisto);
    }




  EfficiencyEvaluator ElectronEfficiencyPt("eleeff-pt.root");
  //  ElectronEfficiencyPt.AddNumerator(GenPt);
  //  ElectronEfficiencyPt.AddNumerator(RecoPt);
  for (int i=0;i<CutIdPt.size();++i)
    ElectronEfficiencyPt.AddNumerator(CutIdPt[i]);
  for (int i=0;i<LHIdPt.size();++i)
    ElectronEfficiencyPt.AddNumerator(LHIdPt[i]);
  for (int i=0;i<CiCIdPt.size();++i)
    ElectronEfficiencyPt.AddNumerator(CiCIdPt[i]);
  ElectronEfficiencyPt.SetDenominator(RecoPt);
  ElectronEfficiencyPt.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronEfficiencyPt.GetCumulativeEfficiencies();

  for (int tightness=0;tightness<CutIdPt.size();++tightness)
    { 
      TCanvas c1;
      efficiency[0+tightness]->SetLineColor(1);
      efficiency[0+tightness]->SetMarkerColor(1);
      efficiency[CutIdPt.size()+tightness]->SetLineColor(2);
      efficiency[CutIdPt.size()+tightness]->SetMarkerColor(2);
      efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetLineColor(4);
      efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetMarkerColor(4);
      //  efficiency[6]->SetLineColor(8);
      
      TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
      leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
      leg->SetFillColor(0);
      leg->AddEntry(efficiency[0+tightness],EgammaCutBasedIDWPs[0+tightness].c_str());
      leg->AddEntry(efficiency[CutIdPt.size()+tightness],EgammaLHBasedIDWPs[0+tightness].c_str());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      leg->AddEntry(efficiency[CutIdPt.size()+LHIdPt.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].c_str());
      
      efficiency[0+tightness]->SetMinimum(0.3);
      efficiency[0+tightness]->SetMaximum(1.0);
      efficiency[0+tightness]->SetMarkerStyle(20);
      efficiency[0+tightness]->SetMarkerSize(1.05);
      efficiency[0+tightness]->SetTitle("");
      efficiency[0+tightness]->GetXaxis()->SetTitle("electron #pt");
      efficiency[0+tightness]->GetYaxis()->SetTitle("Efficiency");
      efficiency[0+tightness]->Draw("PEhist");
      efficiency[CutIdPt.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdPt.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdPt.size()+tightness]->SetMarkerSize(1.05);
      //  efficiency[3+tightness]->Draw("same hist");
      efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetMarkerSize(1.05);
      leg->Draw();
      
      char tightLevel[10];
      sprintf(tightLevel,"%d",tightness);
      c1.SaveAs("eleeff-pt-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-pt-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-pt-tightness"+TString(tightLevel)+".png");
    }


}
