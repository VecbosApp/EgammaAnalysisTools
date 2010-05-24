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
char inputFactor[150];

void makePlots(TH1F *signalHistos1[2], TH1F *signalHistos2[2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal1[150];
  char inputFileNameSignal2[150];
  if ( argc < 4 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignalMC.root fileSignalDATA.root factor" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal1,argv[1]);
  strcpy(inputFileNameSignal2,argv[2]);
  strcpy(inputFactor,argv[3]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile1 = TFile::Open(inputFileNameSignal1);
  signalFile1->cd();

  TH1F *dPhiClassEle1[2];
  TH1F *dEtaClassEle1[2];
  TH1F *EoPClassEle1[2];
  TH1F *HoEClassEle1[2];
  TH1F *sigmaIEtaIEtaClassEle1[2];

  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {
      
    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiClassEle1[iecal]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle1[iecal]        = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"EoPClass_electrons_%d",iecal);
    EoPClassEle1[iecal]      = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEClassEle1[iecal]         = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
    sigmaIEtaIEtaClassEle1[iecal] = (TH1F*)gDirectory->Get(histo); 
    
    dPhiClassEle1[iecal]   -> SetMinimum(0.001);
    dEtaClassEle1[iecal]   -> SetMinimum(0.001);
    EoPClassEle1[iecal]    -> SetMinimum(0.001);
    HoEClassEle1[iecal]   -> SetMinimum(0.001);
    sigmaIEtaIEtaClassEle1[iecal] -> SetMinimum(0.1);
    
  }


  TFile *signalFile2 = TFile::Open(inputFileNameSignal2);
  signalFile2->cd();

  TH1F *dPhiClassEle2[2];
  TH1F *dEtaClassEle2[2];
  TH1F *EoPClassEle2[2];
  TH1F *HoEClassEle2[2];
  TH1F *sigmaIEtaIEtaClassEle2[2];

  for(int iecal=0; iecal<2; iecal++) {
      
    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiClassEle2[iecal]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle2[iecal]        = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"EoPClass_electrons_%d",iecal);
    EoPClassEle2[iecal]      = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEClassEle2[iecal]         = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
    sigmaIEtaIEtaClassEle2[iecal] = (TH1F*)gDirectory->Get(histo); 
    
    dPhiClassEle2[iecal]   -> SetMinimum(0.001);
    dEtaClassEle2[iecal]   -> SetMinimum(0.001);
    EoPClassEle2[iecal]    -> SetMinimum(0.001);
    HoEClassEle2[iecal]   -> SetMinimum(0.001);
    sigmaIEtaIEtaClassEle2[iecal] -> SetMinimum(0.1);
    
  }


  makePlots(dPhiClassEle1,dPhiClassEle2,"dPhiIn","#Delta #phi_{in} (rad)");
  makePlots(dEtaClassEle1,dEtaClassEle2,"dEtaIn","#Delta #eta_{in} (rad)");
  makePlots(EoPClassEle1,EoPClassEle2,"EoP","E_{SC}/p_{in}");
  makePlots(HoEClassEle1,HoEClassEle2,"HoE","H/E");
  makePlots(sigmaIEtaIEtaClassEle1,sigmaIEtaIEtaClassEle2,"sigmaIEtaIEta","#sigma_{#eta#eta} (rad^{2})");

}

void makePlots(TH1F *signalHistos1[2], TH1F *signalHistos2[2], const char *namevar, const char *axistitle) {

  int iptbin=1;

  for(int iecal=0; iecal<2; iecal++) {

    TCanvas c1;
    // c1.SetLogy();

    // MC
    if(signalHistos1[iecal]) {
      signalHistos1[iecal]->SetLineWidth(2);
      signalHistos1[iecal]->SetLineColor(kRed+1);
      signalHistos1[iecal]->SetFillColor(kYellow+1);
      signalHistos1[iecal]->Scale((float)signalHistos2[iecal]->GetSum()/(float)signalHistos1[iecal]->GetSum());
      std::cout << "Normalized to " << (float)signalHistos2[iecal]->GetSum() << std::endl;
    }
    // data
    if(signalHistos2[iecal]) {
      signalHistos2[iecal]->SetLineWidth(2);
      signalHistos2[iecal]->SetLineColor(kBlack);
      signalHistos2[iecal]->SetMarkerColor(kBlack);
      signalHistos2[iecal]->SetMarkerStyle(8);
      // signalHistos2[iecal]->Scale(1.0/signalHistos2[iecal]->Integral());
    }

    double  max=TMath::Max(signalHistos1[iecal]->GetMaximum(),signalHistos2[iecal]->GetMaximum());
    signalHistos1[iecal]->SetMaximum(max+sqrt(max));

    double  min=TMath::Min(signalHistos1[iecal]->GetMinimum(),signalHistos2[iecal]->GetMinimum());
    signalHistos1[iecal]->SetMinimum(-10.);

    signalHistos1[iecal]->SetTitle("");
     signalHistos1[iecal]->GetXaxis()->SetTitle(axistitle);
    signalHistos1[iecal]->GetYaxis()->SetTitle("electrons in 200 nb^{-1}");
    signalHistos1[iecal]->Draw("hist");
    signalHistos2[iecal]->Draw("pE1 same");

    TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(signalHistos1[iecal],"MC");
    leg->AddEntry(signalHistos2[iecal],"data");
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    if(iecal==0) {
      sprintf(epsfilename,"%s_EB_x%s.eps",namevar,inputFactor);
      sprintf(rootfilename,"%s_EB_x%s.root",namevar,inputFactor);
    }
    else if(iecal==1) {
      sprintf(epsfilename,"%s_EE_x%s.eps",namevar,inputFactor); 
      sprintf(rootfilename,"%s_EE_x%s.root",namevar,inputFactor);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
