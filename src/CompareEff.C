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
#include <TLegend.h>
#include <TCanvas.h>

// Offline analysis includes
#include "CommonTools/include/EfficiencyEvaluator.hh"

#include <string>

using namespace std;

std::vector<TString> EgammaCutBasedIDWPs;
std::vector<TString> EgammaCiCBasedIDWPs;
std::vector<TString> EgammaLHBasedIDWPs;

void drawEta(const char* filename, TString type);
void drawPt(const char* filename, TString type);

int main(int argc, char* argv[]) {

  TString IDpart;
   if (argv[3])
     {
       char idpar[100];
       strcpy(idpar,argv[3]);
       IDpart=TString(idpar);
     }
   else
     IDpart="";

   EgammaCutBasedIDWPs.push_back("WP95"+TString(IDpart));;
   EgammaCutBasedIDWPs.push_back("WP90"+TString(IDpart));;
  EgammaCutBasedIDWPs.push_back("WP80"+TString(IDpart));;
  EgammaCutBasedIDWPs.push_back("WP70"+TString(IDpart));;    
  //  EgammaCiCBasedIDWPs.push_back(TString("CiCVeryLoose")+TString(IDpart));;
  //  EgammaCiCBasedIDWPs.push_back(TString("CiCLoose")+TString(IDpart));;
  EgammaCiCBasedIDWPs.push_back(TString("CiCLoose")+TString(IDpart));;
  EgammaCiCBasedIDWPs.push_back(TString("CiCTight")+TString(IDpart));;
  EgammaCiCBasedIDWPs.push_back(TString("CiCSuperTight")+TString(IDpart));;
  EgammaCiCBasedIDWPs.push_back(TString("CiCHyperTight")+TString(IDpart));;
  //  EgammaCiCBasedIDWPs.push_back(TString("CiCHyperTight2")+TString(IDpart));;

  EgammaLHBasedIDWPs.push_back(TString("LHVeryLoose")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHLoose")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHMedium")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHTight")+TString(IDpart));;

  char inputFileNameEta[150];
  char inputFileNamePt[150];
  char type[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: compareEff inputFileEtaSuffix inputFilePtSuffix <optional id part>" << std::endl;
    return 1;
  }
  strcpy(inputFileNameEta,argv[1]);
  strcpy(inputFileNamePt,argv[2]);


  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  char input[300];
//   sprintf(input,"%s.root",inputFileNameEta);
//   drawEta(input,"");
  sprintf(input,"%sHighPt.root",inputFileNameEta);
  drawEta(input,"HighPt");
  sprintf(input,"%sLowPt.root",inputFileNameEta);
  drawEta(input,"LowPt");
//   sprintf(input,"%s.root",inputFileNamePt);
//   drawPt(input, "");
  sprintf(input,"%sBarrel.root",inputFileNamePt);
  drawPt(input, "Barrel");
  sprintf(input,"%sEndcap.root",inputFileNamePt);
  drawPt(input, "Endcap");

}


void drawEta(const char* filename,TString type) {

  TFile *efficiencyFileEta = TFile::Open(filename);

  //  TH1F *GenEta = (TH1F*)efficiencyFileEta->Get("GenEta"+type);
  TH1F *RecoEta = (TH1F*)efficiencyFileEta->Get("RecoEta"+type);

  std::vector<TH1F*> CutIdEta;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta"+type);
      CutIdEta.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFileEta->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta"+type);
//       CutIdEta.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdEta;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta"+type);
      LHIdEta.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFileEta->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta"+type);
//       LHIdEta.push_back(aHisto);
    }
  std::vector<TH1F*> CiCIdEta;
  for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta"+type);
      CiCIdEta.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFileEta->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEta"+type);
//       CiCIdEta.push_back(aHisto);
    }




  EfficiencyEvaluator ElectronEfficiencyEta("eleeff-eta_"+type+".root");
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

  for (int tightness=0;tightness<CutIdEta.size();tightness++)
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
      leg->AddEntry(efficiency[0+tightness],EgammaCutBasedIDWPs[0+tightness].Data());
      leg->AddEntry(efficiency[CutIdEta.size()+tightness],EgammaLHBasedIDWPs[0+tightness].Data());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      leg->AddEntry(efficiency[CutIdEta.size()+LHIdEta.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].Data());
      
      efficiency[0+tightness]->SetMinimum(0.1);
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
      c1.SaveAs("eleeff-eta-"+type+"-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-eta-"+type+"-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-eta-"+type+"-tightness"+TString(tightLevel)+".png");
    }

}

void drawPt(const char* filename, TString type) {
  TFile *efficiencyFilePt = TFile::Open(filename);

  //  TH1F *GenPt = (TH1F*)efficiencyFilePt->Get("GenPt"+type);
  TH1F *RecoPt = (TH1F*)efficiencyFilePt->Get("RecoPt"+type);

  std::vector<TH1F*> CutIdPt;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt"+type);
      CutIdPt.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdPt;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt"+type);
      LHIdPt.push_back(aHisto);
    }
  std::vector<TH1F*> CiCIdPt;
  for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt"+type);
      CiCIdPt.push_back(aHisto);
    }




  EfficiencyEvaluator ElectronEfficiencyPt("eleeff-pt"+type+".root");
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
      leg->AddEntry(efficiency[0+tightness],EgammaCutBasedIDWPs[0+tightness].Data());
      leg->AddEntry(efficiency[CutIdPt.size()+tightness],EgammaLHBasedIDWPs[0+tightness].Data());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      leg->AddEntry(efficiency[CutIdPt.size()+LHIdPt.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].Data());
      
      efficiency[0+tightness]->SetMinimum(0.1);
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
      c1.SaveAs("eleeff-pt-"+type+"-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-pt-"+type+"-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-pt-"+type+"-tightness"+TString(tightLevel)+".png");
    }


}
