#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <RooDataSet.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>

void drawSPlot(const char *varname, const char *axistitle, const char *extracut, int nbins, float min, float max, int logy=1);
void makeSignalSPlots();

void makeSignalSPlots() {

  gSystem->Load("libRooFit");

  drawSPlot("EoPout","E_{seed}/P_{out}","",100,0.0,8.0);
  drawSPlot("HoE","H/E","",100,0.0,0.1);
  drawSPlot("deltaEta","#Delta #eta_{in}","",100,-0.02,0.02);
  drawSPlot("deltaPhi","#Delta #phi_{in}","",100,-0.1,0.1);
  drawSPlot("s9s25","s_{9}/s_{25}","",100,0.5,1.0);
  drawSPlot("sigmaEtaEta","#sigma_{#eta#eta}","",100,0.0,0.05);

}

void drawSPlot(const char *varname, const char *axistitle, const char *extracut, int nbins, float min, float max, int logy) {

  // Load results file...
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0) ;
  gStyle->SetOptTitle(0) ;

  TFile *ZjetsTagProbeMADGRAPH = TFile::Open("/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/ZJetsMADGRAPHPdfsDataset.root");
  RooDataSet *zjetsDataset = (RooDataSet*)ZjetsTagProbeMADGRAPH->Get("T1");
  if(strcmp(extracut,"")!=0) zjetsDataset = (RooDataSet*)zjetsDataset->reduce(extracut);

  TH1F *signalMC = new TH1F(varname, "", nbins, min, max);
  zjetsDataset->tree().Project(varname,varname);

  TFile *DataTagProbe = TFile::Open("MLFit/fitoutput/sPlotsTagAndProbe.root");
  RooDataSet *dataset = (RooDataSet*)DataTagProbe->Get("dataset");
  if(strcmp(extracut,"")!=0) dataset = (RooDataSet*)dataset->reduce(extracut);

  char splotname[200];
  sprintf(splotname,"%s_sPlotSig",varname);
  TH1F *signalsPlot = new TH1F(splotname,"", nbins, min, max);
  dataset->tree().Project(splotname,varname,"N_sig_sw");
  signalsPlot->Sumw2();

  char bkgsplotname[200];
  sprintf(bkgsplotname,"%s_sPplotBkg",varname);
  TH1F *backgroundsPlot = new TH1F(bkgsplotname,"", nbins, min, max);
  dataset->tree().Project(bkgsplotname,varname,"N_bkg_sw");
  backgroundsPlot->Sumw2();

  // sPlot has the correct normalization in 300 pb-1. Normalize the MC to this area
  float integral = signalsPlot->Integral();
  signalMC->Scale(integral/signalMC->Integral());
  backgroundsPlot->Scale(integral/backgroundsPlot->Integral());

  float max = TMath::Max(signalsPlot->GetMaximum(),backgroundsPlot->GetMaximum());
  max = TMath::Max(max,signalMC->GetMaximum());

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(signalMC,"Z+jets MC","l");
  leg->AddEntry(signalsPlot,"signal weighted data","p");
  //  leg->AddEntry(backgroundsPlot,"bkg weighted data","l");

  TCanvas c1;
  c1.SetLogy(1);
  signalsPlot->SetMarkerStyle(8);
  signalsPlot->SetMarkerSize(1);
  signalsPlot->GetXaxis()->SetTitle(axistitle);
  signalsPlot->GetYaxis()->SetTitle("events in 300 pb^{-1}");
  signalsPlot->SetMaximum(max+sqrt(max));
//   backgroundsPlot->SetLineStyle(kDashed);
  signalsPlot->Draw("pe1");
  signalMC->Draw("same hist");
//   backgroundsPlot->Draw("same hist");
  leg->Draw();

  char epsfilename[200];
  sprintf(epsfilename,"%s_sPlot.eps",varname);
  char rootfilename[200];
  sprintf(rootfilename,"%s_sPlot.root",varname);

  c1.SaveAs(epsfilename);
  c1.SaveAs(rootfilename);

}
  
