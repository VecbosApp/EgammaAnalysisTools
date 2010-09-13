#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

#include <iostream>

#define NSPECIES 6
#define NVARIABLES 3
#define NCUTS 1

void makeDataMCPlots()
{
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x tdrstyle.C");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetOptTitle(0); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetMarkerColor(1);

  TString species[NSPECIES];
  species[0]="Data";
  species[1]="QCD";
  species[2]="#gamma jet";
  species[3]="W";
  species[4]="Z";
  species[5]="top";

  Color_t colors[NSPECIES];
  colors[0]=kBlack;
  colors[1]=kAzure+8;
  colors[2]=kYellow-7;
  colors[3]=kPink+6;
  colors[4]=kGreen; // chiara
  colors[5]=kRed;   // chiara

  Color_t lineColors[NSPECIES];
  lineColors[0]=kBlack;
  lineColors[1]=kAzure+4;
  lineColors[2]=kYellow+4;
  lineColors[3]=kMagenta+3;
  lineColors[4]=kGreen; // chiara
  lineColors[5]=kRed;   // chiara

  int legendOrder[NSPECIES];
  legendOrder[0]=0;
  legendOrder[1]=3;
  legendOrder[2]=2;
  legendOrder[3]=1;
  legendOrder[4]=4;
  legendOrder[5]=5;

  TString files[NSPECIES];
  
  files[0]="results_data/dataset_jetmettau_mergedTree.root";
  files[1]="results/trees_fullStat/QCD_tree.root";
  files[2]="results/trees_fullStat/GammaJets_tree.root";
  files[3]="results/trees_fullStat/WJetsMADGRAPH_tree.root";
  files[4]="results/trees_fullStat/ZJetsMADGRAPH_tree.root";
  files[5]="results/trees_fullStat/TOP_tree.root";

  TString plotsDir="./egammaQCD/";

  TFile* fOut=new TFile("egamma.root","RECREATE");
  
  char icut[NCUTS][100];
  TH1F* histos[NSPECIES][NCUTS][NVARIABLES]; //6 species, 10 cut levels, 3 variables
  
  TString variables[NVARIABLES];
  variables[0]="qcdInvmass";
  variables[1]="qcdDeltaphi";
  variables[2]="qcdMet";

  TString units[NVARIABLES];
  units[0]="GeV/c^{2}";
  units[1]="";
  units[2]="GeV";

  int nbins[NVARIABLES];
  nbins[0]=30;
  nbins[1]=30;
  nbins[2]=30;

  float range[NVARIABLES][2]; // 4 variables, min, max
  // invariant mass
  range[0][0]=0.;
  range[0][1]=500.;
  // deltaphi
  range[1][0]=0.;
  range[1][1]=3.14;
  // met
  range[2][0]=0.;
  range[2][1]=500.;

  TString xaxisLabel[NVARIABLES];
  xaxisLabel[0]="m_{ee}";
  xaxisLabel[1]="#Delta #phi";
  xaxisLabel[2]="met";

  TString binSize[NVARIABLES];

  for (int z=0;z<NVARIABLES;++z) {
    for (int j=0;j<NCUTS;++j){
      sprintf(icut[j],"icut%d",j);
      for (int i=0;i<NSPECIES;++i){
	histos[i][j][z]=new TH1F(variables[z]+"_W_"+species[i]+"_"+TString(icut[j]),variables[z]+"_W_"+species[i]+"_"+TString(icut[j]),nbins[z],range[z][0],range[z][1]);
	char binsiz[10];
	sprintf(binsiz,"%2.0f",(range[z][1]-range[z][0])/nbins[z]);
	binSize[z]=TString(binsiz);
      }
    }
  }

  TString cut[NCUTS];
  cut[0]="";

  TString intLumi="1.97";
  TFile *_file[NSPECIES];
  TTree *T1[NSPECIES];

  TCanvas* c1= new TCanvas("test","test",800,800);

  for (int i=0;i<NSPECIES;++i) {
    _file[i]=TFile::Open(files[i]);
    T1[i] = (TTree*)_file[i]->Get("T1");
  }
  
  int nspeciesToRun=NSPECIES;
  
  for (int z=0;z<NVARIABLES;++z) {
    for (int j=0;j<NCUTS;++j) {
      for (int i=0;i<nspeciesToRun;++i) {
	fOut->cd();
	TString histoName=variables[z]+"_W_"+species[i]+"_"+TString(icut[j]);
	std::cout << "Producing " << histoName << std::endl;
	if (T1[i]==0) {
	  std::cout << "Tree not found" << std::endl;
	  return;
	}
	if (i!=0)
	  T1[i]->Project(histoName,variables[z],cut[j]+"weight*"+intLumi);
	else
	  T1[i]->Project(histoName,variables[z],cut[j]+"1.");
	std::cout << "Done " << histoName << std::endl;
      }
      
      THStack histo_MC(variables[z]+"_MC",variables[z]+"_MC");
      for (int i=1;i<nspeciesToRun;++i) {
	histos[i][j][z]->SetFillColor(colors[i]);
	histos[i][j][z]->SetLineColor(lineColors[i]);
	histo_MC.Add(histos[i][j][z]);
      }
      
      float maximum=TMath::Max(histo_MC.GetMaximum(),histos[0][j][z]->GetMaximum());
      histo_MC.SetMinimum(0.01);
      histo_MC.SetMaximum(maximum*2.);
      histo_MC.Draw("");
      histo_MC.GetXaxis()->SetTitle(xaxisLabel[z]+" ["+units[z]+"]");
      histo_MC.GetYaxis()->SetTitle("Events/"+binSize[z]+" "+units[z]);
      
      histos[0][j][z]->SetMarkerStyle(20);
      histos[0][j][z]->SetMarkerSize(1.3);
      histos[0][j][z]->Draw("EPSAME");
      TPaveText pt1(0.6,0.83,0.8,0.9,"NDC");
      //	  pt1.SetTextFont(72);
      pt1.SetTextSize(0.028);
      pt1.SetTextAlign(12);
      pt1.SetFillColor(0);
      pt1.SetBorderSize(0);
      // pt1.AddText("CMS Preliminary 2010");
      // pt1.AddText("");
      // pt1.AddText("");
      // pt1.AddText("");
      // pt1.AddText("#sqrt{s}=7 TeV L_{int}="+intLumi+" pb^{-1}");
      // pt1.Draw();
      
      c1->Update();
      TLegendEntry *legge;
      TLegend *leg;
      leg = new TLegend(0.6,0.65,0.93,0.8,cut[j]);
      leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.025);
      leg->SetFillColor(0);
      for (int i=0;i<nspeciesToRun;++i) {
	if (i == 0)
	  legge = leg->AddEntry(histos[legendOrder[i]][j][z],"Data", "lpe");
	else
	  legge = leg->AddEntry(histos[legendOrder[i]][j][z],species[legendOrder[i]],"f");
      }
      leg->Draw();
      c1->SetLogy(0);
      c1->Update();
      c1->SaveAs(plotsDir+variables[z]+"DataMC_"+TString(icut[j])+".png");
      c1->SaveAs(plotsDir+variables[z]+"DataMC_"+TString(icut[j])+".root");
      c1->SaveAs(plotsDir+variables[z]+"DataMC_"+TString(icut[j])+".C");
      histo_MC.SetMaximum(maximum*100);
      c1->SetLogy(1);
      c1->Update();
      c1->SaveAs(plotsDir+variables[z]+"DataMC_"+TString(icut[j])+"_log.png");
    }
  }
  
  fOut->Write();
  fOut->Close();
  
}
