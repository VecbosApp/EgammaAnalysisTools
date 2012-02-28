#include "TCanvas.h"
#include "TFile.h"
#include "RooHist.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TLine.h"
#include "TString.h"
#include "TGraphErrors.h"

#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <math.h>

using namespace std;

void compareEff() {

  vector<TString> filesHWW, filesHZZ;
  filesHWW.push_back(TString("results/efficiency_hww2011_Data2011_Pt10To20.root"));
  filesHWW.push_back(TString("results/efficiency_hww2011_Data2011_Pt20To200.root"));
  filesHZZ.push_back(TString("results/efficiency_hzz2011_Data2011_Pt10To20.root"));
  filesHZZ.push_back(TString("results/efficiency_hzz2011_Data2011_Pt20To200.root"));

  vector<TString> etabins;
  char etabin[100];
  for(int i=0;i<5;++i) {
    sprintf(etabin,"pt_PLOT_abseta_bin%d",i);
    etabins.push_back(TString(etabin));
  }

  TString dir = TString("eleIDdir/etapt/fit_eff_plots/");

  vector<double> pts, effww, effzz, erreffww, erreffzz;
  vector<int> etas;

  for(int ifile=0; ifile<2; ++ifile) { // two separate fits are done for pT</> 20 GeV
    TFile *fileHWW = TFile::Open(filesHWW[ifile]);
    TFile *fileHZZ = TFile::Open(filesHZZ[ifile]);
    for(int ieta=0;ieta<5;++ieta) { // loop over eta bins
      TString plot = dir+etabins[ieta];
      TCanvas *cHWW = (TCanvas*)fileHWW->Get(plot);
      TCanvas *cHZZ = (TCanvas*)fileHZZ->Get(plot);
      RooHist* hWW = (RooHist*)cHWW->GetPrimitive("hxy_fit_eff");
      RooHist* hZZ = (RooHist*)cHZZ->GetPrimitive("hxy_fit_eff");

      // get the points
      Int_t nbins = hWW->GetN();
      Double_t* x = hWW->GetX();
      Double_t* yWW = hWW->GetY();
      Double_t* yZZ = hZZ->GetY();
      Double_t* eyWW = hWW->GetEYlow();
      Double_t* eyZZ = hZZ->GetEYlow();
      
      for(int ipt=0;ipt<nbins;++ipt) { // loop over pt bins
	pts.push_back(x[ipt]);
	etas.push_back(ieta);
	effww.push_back(yWW[ipt]);
	erreffww.push_back(eyWW[ipt]);
	effzz.push_back(yZZ[ipt]);
	erreffzz.push_back(eyZZ[ipt]);
      }
    }
  }

  vector<TGraphErrors*> graphs;
  for(int ieta=0;ieta<5;++ieta) {
    char name[200];
    sprintf(name,"effratio_abseta_bin%d",ieta);
    TGraphErrors *graph =  new TGraphErrors(5);
    graph->SetName(name);
    graph->SetTitle("");
    int ip=0;
    for(int i=0;i<(int)pts.size();++i) {
      if(etas[i]!=ieta) continue; // bins are in one 1D vector. Keep only the actual eta bin
      float ratio = effzz[i]/effww[i];
      float errratio = ratio * sqrt(pow(erreffzz[i]/effzz[i],2)+pow(erreffww[i]/effww[i],2));
      cout << "eta bin = " << ieta << "ipt = " << ip << "\tpt = " << pts[i] << "\tratio zz/ww = " << ratio << " +/- " << errratio << endl;
      graph->SetPoint(ip,float(pts[i]),ratio);
      graph->SetPointError(ip,0,errratio);
      ip++;

      // cosmetics
      graph->SetMinimum(0.7);
      graph->SetMaximum(2.0);
      graph->SetMarkerSize(2);
      graph->SetMarkerStyle(20);
      graph->SetMarkerColor(kAzure-6);
      graph->SetLineColor(kAzure-6);
      graph->GetXaxis()->SetTitle("p_{T} [GeV]");
      graph->GetYaxis()->SetTitle("efficiency ratio");
    }
    graphs.push_back(graph);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  for(int ieta=0;ieta<5;++ieta) {
    TString name = (TString)graphs[ieta]->GetName();
    graphs[ieta]->Draw("AP");
    TLine *l = new TLine(c1->GetUxmin(),1,c1->GetUxmax(),1);
    l->SetLineWidth(2);
    l->SetLineColor(kRed+1);
    l->Draw("same");
    c1->SaveAs(name+TString(".pdf"));
    c1->SaveAs(name+TString(".png"));
  }

}
