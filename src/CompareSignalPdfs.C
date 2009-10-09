// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

void makePlots(TH1F *signalHistos1[2][2][2], TH1F *signalHistos2[2][2][2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal1[150];
  char inputFileNameSignal2[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignal1.root fileSignal2.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal1,argv[1]);
  strcpy(inputFileNameSignal2,argv[2]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile1 = TFile::Open(inputFileNameSignal1);
  signalFile1->cd();

  //[iecal:0=EB;1=EE][iptbin:0=<15GeV;1=>15GeV][iclass:0=nonshowering;1=showering]
  TH1F *dPhiClassEle1[2][2][2];
  TH1F *dEtaClassEle1[2][2][2];
  TH1F *EoPoutClassEle1[2][2][2];
  TH1F *EoPClassEle1[2][2][2];
  TH1F *HoEClassEle1[2][2][2];
  TH1F *sigmaIEtaIEtaClassEle1[2][2][2];
  TH1F *s9s25ClassEle1[2][2][2];
  TH1F *s1s9ClassEle1[2][2][2];

  int iptbin=1; // draw only >15 GeV

  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {
    
    for(int iclass=0; iclass<2; iclass++) {
      
      sprintf(histo,"dPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      dPhiClassEle1[iecal][iptbin][iclass]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      dEtaClassEle1[iecal][iptbin][iclass]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      EoPoutClassEle1[iecal][iptbin][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"EoPClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      EoPClassEle1[iecal][iptbin][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      HoEClassEle1[iecal][iptbin][iclass]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      sigmaIEtaIEtaClassEle1[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
      s9s25ClassEle1[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"s1s9Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
      s1s9ClassEle1[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo);

    }

  }


  TFile *signalFile2 = TFile::Open(inputFileNameSignal2);
  signalFile2->cd();

  TH1F *dPhiClassEle2[2][2][2];
  TH1F *dEtaClassEle2[2][2][2];
  TH1F *EoPoutClassEle2[2][2][2];
  TH1F *EoPClassEle2[2][2][2];
  TH1F *HoEClassEle2[2][2][2];
  TH1F *sigmaIEtaIEtaClassEle2[2][2][2];
  TH1F *s9s25ClassEle2[2][2][2];
  TH1F *s1s9ClassEle2[2][2][2];

  for(int iecal=0; iecal<2; iecal++) {

    for(int iclass=0; iclass<2; iclass++) {

      sprintf(histo,"dPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      dPhiClassEle2[iecal][iptbin][iclass]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      dEtaClassEle2[iecal][iptbin][iclass]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      EoPoutClassEle2[iecal][iptbin][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"EoPClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      EoPClassEle2[iecal][iptbin][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      HoEClassEle2[iecal][iptbin][iclass]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      sigmaIEtaIEtaClassEle2[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
      s9s25ClassEle2[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"s1s9Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
      s1s9ClassEle2[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo);
      
      dPhiClassEle2[iecal][iptbin][iclass]->Sumw2();
      dEtaClassEle2[iecal][iptbin][iclass]->Sumw2();
      EoPoutClassEle2[iecal][iptbin][iclass]->Sumw2();
      EoPClassEle2[iecal][iptbin][iclass]->Sumw2();
      HoEClassEle2[iecal][iptbin][iclass]->Sumw2();
      sigmaIEtaIEtaClassEle2[iecal][iptbin][iclass]->Sumw2();
      s9s25ClassEle2[iecal][iptbin][iclass]->Sumw2();
      s1s9ClassEle2[iecal][iptbin][iclass]->Sumw2();
      
    }

  }

  makePlots(dPhiClassEle1,dPhiClassEle2,"dPhiIn","#Delta #phi_{in} (rad)");
  makePlots(dEtaClassEle1,dEtaClassEle2,"dEtaIn","#Delta #eta_{in} (rad)");
  makePlots(EoPoutClassEle1,EoPoutClassEle2,"EoPout","E_{seed}/p_{out}");
  makePlots(EoPClassEle1,EoPClassEle2,"EoP","E_{SC}/p_{in}");
  makePlots(HoEClassEle1,HoEClassEle2,"HoE","H/E");
  makePlots(sigmaIEtaIEtaClassEle1,sigmaIEtaIEtaClassEle2,"sigmaIEtaIEta","#sigma_{#eta#eta} (rad^{2})");
  makePlots(s9s25ClassEle1,s9s25ClassEle2,"s9s25","s_{9}/s_{25}");
  makePlots(s1s9ClassEle1,s1s9ClassEle2,"s1s9","s_{1}/s_{9}");

}

void makePlots(TH1F *signalHistos1[2][2][2], TH1F *signalHistos2[2][2][2], const char *namevar, const char *axistitle) {

  int iptbin=1;

  for(int iecal=0; iecal<2; iecal++) {

    TCanvas c1;
    c1.SetLogy();

    if(signalHistos1[iecal][iptbin][0]) {
      signalHistos1[iecal][iptbin][0]->SetLineWidth(2);
      signalHistos1[iecal][iptbin][0]->SetLineColor(kRed+1);
      signalHistos1[iecal][iptbin][0]->Scale(1.0/signalHistos1[iecal][iptbin][0]->Integral());
    }
    if(signalHistos1[iecal][iptbin][1]) {
      signalHistos1[iecal][iptbin][1]->SetLineWidth(2);
      signalHistos1[iecal][iptbin][1]->SetLineColor(kBlue+3);
      signalHistos1[iecal][iptbin][1]->Scale(1.0/signalHistos1[iecal][iptbin][1]->Integral());
    }
    if(signalHistos2[iecal][iptbin][0]) {
      signalHistos2[iecal][iptbin][0]->SetLineWidth(2);
      signalHistos2[iecal][iptbin][0]->SetLineColor(kRed+1);
      signalHistos2[iecal][iptbin][0]->SetLineStyle(kDashed);
      signalHistos2[iecal][iptbin][0]->Scale(1.0/signalHistos2[iecal][iptbin][0]->Integral());
    }
    if(signalHistos2[iecal][iptbin][1]) {
      signalHistos2[iecal][iptbin][1]->SetLineWidth(2);
      signalHistos2[iecal][iptbin][1]->SetLineColor(kBlue+3);
      signalHistos2[iecal][iptbin][1]->SetLineStyle(kDashed);
      signalHistos2[iecal][iptbin][1]->Scale(1.0/signalHistos2[iecal][iptbin][1]->Integral());
    }

//     double  max=TMath::Max(signalHistos1[iecal][iptbin][0]->GetMaximum(),signalHistos1[iecal][iptbin][1]->GetMaximum());
//     max=TMath::Max(max,backgroundHistos[iecal][iptbin]->GetMaximum());
//     signalHistos1[iecal][iptbin][0]->SetMaximum(max+sqrt(max));
//     signalHistos1[iecal][iptbin][1]->SetMaximum(max+sqrt(max));
//     backgroundHistos[iecal][iptbin]->SetMaximum(max+sqrt(max));

    signalHistos1[iecal][iptbin][0]->SetTitle("");
    signalHistos1[iecal][iptbin][0]->GetXaxis()->SetTitle(axistitle);
    signalHistos1[iecal][iptbin][0]->GetYaxis()->SetTitle("a.u.");
    signalHistos1[iecal][iptbin][0]->Draw();
    signalHistos1[iecal][iptbin][1]->Draw("same");
    signalHistos2[iecal][iptbin][0]->Draw("same pe1");
    signalHistos2[iecal][iptbin][1]->Draw("same pe1");

    TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(signalHistos1[iecal][iptbin][0],"MC non-showering");
    leg->AddEntry(signalHistos1[iecal][iptbin][1],"MC showering");
    leg->AddEntry(signalHistos2[iecal][iptbin][0],"data non-showering");
    leg->AddEntry(signalHistos2[iecal][iptbin][1],"data showering");
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    if(iecal==0) {
      sprintf(epsfilename,"%s_EB.eps",namevar);
      sprintf(rootfilename,"%s_EB.root",namevar);
    }
    else if(iecal==1) {
      sprintf(epsfilename,"%s_EE.eps",namevar); 
      sprintf(rootfilename,"%s_EE.root",namevar);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
