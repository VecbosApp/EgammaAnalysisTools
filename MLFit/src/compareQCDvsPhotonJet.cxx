#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <RooDataSet.h>

using namespace std;

void makePlot(const char* name, const char* title, float min, float max, int iecal, int nbins=50);

void makeAllPlots() {

  // barrel
  makePlot("EoPout","E_{seed}/P_{out}",0,20,0);
  makePlot("HoE","H/E",0,1,0);
  makePlot("deltaEta","#DELTA #eta",-0.02,0.02,0);
  makePlot("deltaPhi","#DELTA #phi",-0.1,0.1,0);
  makePlot("s9s25","s9 / s25",0,1,0);
  makePlot("sigmaIEtaIEta","#sigma_{i#eta i#eta}",0,0.1,0);

  // endcap
  makePlot("EoPout","E_{seed}/P_{out}",0,20,1);
  makePlot("HoE","H/E",0,1,1);
  makePlot("deltaEta","#DELTA #eta",-0.02,0.02,1);
  makePlot("deltaPhi","#DELTA #phi",-0.1,0.1,1);
  makePlot("s9s25","s9 / s25",0,1,1);
  makePlot("sigmaIEtaIEta","#sigma_{i#eta i#eta}",0,0.1,1);

}

void makePlot(const char* name, const char* title, float min, float max, int iecal, int nbins) {

  TFile *fileQCD = TFile::Open("datasets_QCDTaP/qcdSignal.root");
  RooDataSet *dataQCD = (RooDataSet*) fileQCD->Get("T1");
  fileQCD->Close();

  TFile *filePhotonJet = TFile::Open("datasets_QCDTaP/photonJet.root");
  RooDataSet *dataPhotonJet = (RooDataSet*)filePhotonJet->Get("T1");
  filePhotonJet->Close();

  char nameQCD[200];
  sprintf(nameQCD, "%s_qcd", name);

  char namePhotonJet[200];
  sprintf(namePhotonJet, "%s_photonjet", name);

  TH1F *hQCD = new TH1F(nameQCD, title, nbins, min, max);
  TH1F *hPhotonJet = new TH1F(namePhotonJet, title, nbins, min, max);

  TTree *treeQCD = (TTree*)&dataQCD->tree();
  TTree *treePhotonJet = (TTree*)&dataPhotonJet->tree();

  char cut[200];
  sprintf(cut,"iecal==%d",iecal);

  treeQCD->Project(nameQCD,name,cut);
  treePhotonJet->Project(namePhotonJet,name,cut);

  hPhotonJet->Sumw2();
  hQCD->Scale( 1./float(hQCD->Integral()) );
  hPhotonJet->Scale( 1./float(hPhotonJet->Integral()) );
  
  hQCD->SetLineColor(kBlue+3);
  hQCD->SetLineWidth(2);
  hPhotonJet->SetLineWidth(2);
  
  hQCD->GetXaxis()->SetTitle(title);
  hQCD->GetYaxis()->SetTitle("a.u.");

  TCanvas c("c","c",400,400);
  hQCD->Draw();
  hPhotonJet->Draw("same pe1");

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(hQCD,"fake from QCD di-jets","l");
  leg->AddEntry(hPhotonJet,"fake from #gamma + jets,","pl");
  leg->Draw();

  char nameFig[200];
  if(iecal==0) sprintf(nameFig,"%s_qcd_vs_photonjet_EB.eps", name);
  else if(iecal==1) sprintf(nameFig,"%s_qcd_vs_photonjet_EE.eps", name);
  else sprintf(nameFig,"iecalnotset.dummy");
  c.SaveAs(nameFig);

}
