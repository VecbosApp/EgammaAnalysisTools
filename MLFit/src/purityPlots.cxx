#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <RooDataSet.h>

void drawKinematics() {

  // datasets
  std::vector<TFile*> datasets;
  TFile *qcdSignal = TFile::Open("/u1/crovelli/data/Like3.2.X/datasets_QCDTaP/qcdSignal.root");
  TFile *ttJets    = TFile::Open("/u1/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root");
  TFile *zee       = TFile::Open("/u1/crovelli/data/Like3.2.X/datasets_QCDTaP/zee.root");
  TFile *wenu      = TFile::Open("/u1/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu.root");
  TFile *gammaJets = TFile::Open("/u1/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet.root");
  datasets.push_back(qcdSignal); 
  datasets.push_back(ttJets); 
  datasets.push_back(zee); 
  datasets.push_back(wenu); 
  datasets.push_back(gammaJets); 
 
  // histograms
  std::vector<TH1F*> met;
  std::vector<TH1F*> mll;
  std::vector<TH1F*> deltaPhi;
  TH1F* metH      = new TH1F("metH",      "metH",     50,0.0,100.0);
  TH1F* mllH      = new TH1F("mllH",      "mllH",     50,0.0,600.0);
  TH1F* deltaPhiH = new TH1F("deltaPhiH", "deltaPhiH",50,0.5,3.1415);

  // projecting trees into histos
  for(int i=0; i<(int)datasets.size(); i++) {
    
    RooDataSet *theRoodataset = (RooDataSet*)datasets[i]->Get("T1");
    TTree *tree = (TTree*)&theRoodataset->tree();      
    std::cout << "dataset " << i << " has " << tree->GetEntries() << " entries" << endl;

    // met
    char buf[50];
    sprintf(buf,"met_%d",i);
    TH1F* metProcessX = (TH1F*) metH->Clone(buf);
    met.push_back(metProcessX);
    tree->Project(buf,"qcdMet","weight");

    // mll
    sprintf(buf,"mll_%d",i);
    TH1F* mllProcessX = (TH1F*) mllH->Clone(buf);
    mll.push_back(mllProcessX);
    tree->Project(buf,"qcdInvmass","weight");

    // deltaPhi
    sprintf(buf,"deltaPhi_%d",i);
    TH1F* deltaPhiProcessX = (TH1F*) deltaPhiH->Clone(buf);
    deltaPhi.push_back(deltaPhiProcessX);
    tree->Project(buf,"qcdDeltaphi","weight");
  }


  // summary
  cout << "Numbers: "  << endl;
  cout << "signal: "   << met[0]->Integral() << endl;
  cout << "ttbar: "    << met[1]->Integral() << endl;
  cout << "zee: "      << met[2]->Integral() << endl;
  cout << "wenu: "     << met[3]->Integral() << endl;
  cout << "gammaJet: " << met[4]->Integral() << endl;


  // plots
  TLegend *leg = new TLegend(0.11,0.65,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(met[0],"QCD","f");
  leg->AddEntry(met[1],"ttbar","f");
  leg->AddEntry(met[2],"Z->ee","f");
  leg->AddEntry(met[3],"W->enu","f");
  leg->AddEntry(met[4],"gamma+jets","f");

  gStyle->SetOptStat(0);

  // --------------
  // draw met
  TCanvas cmet("cmet","cmet",600,600);
  cmet.SetLogy();
  met[0]->SetTitle("");
  met[0]->GetXaxis()->SetTitle("Missing E_{T} [GeV]");
  met[0]->Draw();
  met[0]->SetFillColor(1);      met[0]->Draw("same");
  met[4]->SetFillColor(5);      met[4]->Draw("same");
  met[3]->SetFillColor(4);      met[3]->Draw("same");  
  met[1]->SetFillColor(2);      met[1]->Draw("same");
  met[2]->SetFillColor(3);      met[2]->Draw("same");
  leg->Draw();
  cmet.SaveAs("cmet.root");


  // --------------  
  // draw mll
  TCanvas cmll("cmll","cmll",600,600);
  cmll.SetLogy();
  mll[0]->SetTitle("");
  mll[0]->GetXaxis()->SetTitle("Invariant mass [GeV]");
  mll[0]->SetFillColor(1);      mll[0]->Draw();  
  mll[4]->SetFillColor(5);      mll[4]->Draw("same");
  mll[3]->SetFillColor(4);      mll[3]->Draw("same");
  mll[1]->SetFillColor(2);      mll[1]->Draw("same");
  mll[2]->SetFillColor(3);      mll[2]->Draw("same");  
  leg->Draw();
  cmll.SaveAs("cmll.root");

  // --------------
  // draw deltaPhi
  TCanvas cdeltaPhi("cdeltaPhi","cdeltaPhi",600,600);
  cdeltaPhi.SetLogy();
  deltaPhi[0]->SetTitle("");
  deltaPhi[0]->GetXaxis()->SetTitle("#Delta #phi");
  deltaPhi[0]->SetFillColor(1);      deltaPhi[0]->Draw();  
  deltaPhi[4]->SetFillColor(5);      deltaPhi[4]->Draw("same");
  deltaPhi[3]->SetFillColor(4);      deltaPhi[3]->Draw("same");
  deltaPhi[1]->SetFillColor(2);      deltaPhi[1]->Draw("same");
  deltaPhi[2]->SetFillColor(3);      deltaPhi[2]->Draw("same");  
  leg->Draw();
  cdeltaPhi.SaveAs("cdeltaPhi.root");

}



