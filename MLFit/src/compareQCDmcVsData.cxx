#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

void makePlot(const char* name, const char* title, float min, float max, int iecal, int nbins=50);

void makeAllPlots() {

  // barrel
  makePlot("EoPout","E_{seed}/P_{out}",0,5,0);
  makePlot("EoP",   "E/P",             0,5,0);
  makePlot("HoE",   "H/E",             0,1,0);
  makePlot("deta",  "#Delta #eta",-0.02,0.02,0);
  makePlot("dphi",  "#Delta #phi",-0.1,0.1,0);
  makePlot("see",   "#sigma_{i#eta i#eta}",0,0.05,0);
  makePlot("pt",    "p_{T}",0,100.,0);
  makePlot("eta",   "#eta",-2.5,2.5,0);
  makePlot("fbrem", "fbrem",-0.2,1.,0);

  // endcap
  makePlot("EoPout","E_{seed}/P_{out}",0,5,1);
  makePlot("EoP",   "E/P",             0,5,1);
  makePlot("HoE",   "H/E",             0,1,1);
  makePlot("deta",  "#Delta #eta",-0.02,0.02,1);
  makePlot("dphi",  "#Delta #phi",-0.1,0.1,1);
  makePlot("see",   "#sigma_{i#eta i#eta}",0,0.05,1);
  makePlot("pt",    "p_{T}",0,100.,1);
  makePlot("eta",   "#eta",-2.5,2.5,1);
  makePlot("fbrem", "fbrem",-0.2,1.,1);

  // all
  makePlot("pt",  "p_{T}",0,100.,2);
  makePlot("eta", "#eta",-2.5,2.5,2);
}

void makePlot(const char* name, const char* title, float min, float max, int iecal, int nbins) {
  
  // TFile *fileQCD = TFile::Open("results/trees/QCD_tree.root");
  TFile *fileQCD = TFile::Open("results/trees/QCD_Pt-20_TuneD6T_tree.root");
  TTree *treeQCD = (TTree*) fileQCD->Get("T1");
  
  TFile *filePhotonJet = TFile::Open("results_data/dataset_jetmettau_mergedTree.root");
  TTree *treePhotonJet = (TTree*)filePhotonJet->Get("T1");

    
  char nameQCD[200];
  sprintf(nameQCD, "%s_qcd", name);
  
  char namePhotonJet[200];
  sprintf(namePhotonJet, "%s_photonjet", name);

  TH1F *hQCD       = new TH1F(nameQCD,       title, nbins, min, max);
  TH1F *hPhotonJet = new TH1F(namePhotonJet, title, nbins, min, max);

  if(iecal==0)  treeQCD->Project(nameQCD,name,"weight*(abs(eta)<1.476)*(pt>20)");
  if(iecal==1)  treeQCD->Project(nameQCD,name,"weight*(abs(eta)>1.476)*(pt>20)");
  if(iecal==2)  treeQCD->Project(nameQCD,name,"weight*(pt>20)");

  if (iecal==0) treePhotonJet->Project(namePhotonJet,name,"(abs(eta)<1.476)*(pt>20)");
  if (iecal==1) treePhotonJet->Project(namePhotonJet,name,"(abs(eta)>1.476)*(pt>20)");
  if (iecal==2) treePhotonJet->Project(namePhotonJet,name,"(pt>20)");

  float intQCD = hQCD->Integral();
  float intPJ  = hPhotonJet->Integral();
  hQCD->Sumw2();
  hQCD->Scale( 1./intQCD );
  hPhotonJet->Sumw2();
  hPhotonJet->Scale( 1./intPJ );
  cout << "normalized = " << hQCD->Integral() << " " << hPhotonJet->Integral() << endl;

  hQCD->SetFillColor(3);
  hQCD->SetLineWidth(2);
  hPhotonJet->SetLineWidth(2);
  
  hQCD->GetXaxis()->SetTitle(title);
  hQCD->GetYaxis()->SetTitle("a.u.");

  TCanvas c = new TCanvas("c","c",400,400);
  hQCD->Draw("hist");
  hPhotonJet->Draw("same pe1");

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); 
  leg->SetBorderSize(0); 
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(hQCD,"MC","f");
  leg->AddEntry(hPhotonJet,"data","pl");
  leg->Draw();

  char nameFig[200];
  if(iecal==0) sprintf(nameFig,"%s_EB.eps", name);
  if(iecal==1) sprintf(nameFig,"%s_EE.eps", name);
  if(iecal==2) sprintf(nameFig,"%s_all.eps", name);
  c.SaveAs(nameFig);
  if(iecal==0) sprintf(nameFig,"%s_EB.root", name);
  if(iecal==1) sprintf(nameFig,"%s_EE.root", name);
  if(iecal==2) sprintf(nameFig,"%s_all.root", name);
  c.SaveAs(nameFig);
}
