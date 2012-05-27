#include "TCanvas.h"
#include "TFile.h"
#include "RooHist.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TLine.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <math.h>

using namespace std;

void plotTPeff() {

  vector<TString> filesData, filesMC;
  filesData.push_back("results/hzzref/newhzzWP_eff_Data7TeV_Kine_Pt10To1000.root");
  filesMC.push_back("results/hzzref/newhzzWP_eff_MC7TeV_Kine_Pt10To1000.root");

  map<int,TString> etabins;
  etabins.insert(make_pair(0,"0<|#eta|<0.8"));
  etabins.insert(make_pair(1,"0.8<|#eta|<1.4442"));
  etabins.insert(make_pair(2,"1.4442<|#eta|<1.566"));
  etabins.insert(make_pair(3,"1.566<|#eta|<1.8"));
  etabins.insert(make_pair(4,"1.8<|#eta|<2.0"));
  etabins.insert(make_pair(5,"2.0<|#eta|<2.5"));


  for(unsigned int e=0;e<etabins.size();e++) {

    char filename[100];
    sprintf(filename,"eff_etabin%d",e);
    
    vector<TString> files;
    files.push_back(filename);
    
    TString dirdata = TString("eleIDdir/etapt/fit_eff_plots/");
    TString dirmc = TString("eleIDdir/etapt/cnt_eff_plots/");
    char plotname[100];
    sprintf(plotname,"pt_PLOT_abseta_bin%d_&_vertices_bin0",e);
    TString plot(plotname);
    
    TCanvas c1("c1","c1",600,600);
    c1.SetGridx();
    c1.SetGridy();
    TLegend* legend = new TLegend(0.60, 0.16, 0.83, 0.32);
    legend->SetBorderSize(   0);
    legend->SetFillColor (   0);
    legend->SetTextAlign (  12);
    legend->SetTextFont  (  42);
    legend->SetTextSize  (0.05);
    
    vector<RooHist*> plotdata, plotmc;
    for(int i=0;i<(int)filesData.size();++i) {
      TFile *fileData = TFile::Open(filesData[i]);
      TFile *fileMC = TFile::Open(filesMC[i]);

      TCanvas *cdata = (TCanvas*)fileData->Get(dirdata+plot);
      TCanvas *cmc = (TCanvas*)fileMC->Get(dirmc+plot);

      plotdata.push_back((RooHist*)cdata->GetPrimitive("hxy_fit_eff"));
      plotmc.push_back((RooHist*)cmc->GetPrimitive("hxy_cnt_eff"));

      c1.cd();
    
      plotdata[i]->SetMinimum(0.4);
      plotdata[i]->SetMaximum(1.2);
      plotdata[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      plotdata[i]->GetYaxis()->SetTitle("efficiency");
      plotdata[i]->GetXaxis()->SetRangeUser(0,80);

      // cosmetics

      // this is to cure the features of the error bars of the saved RooPlot
      int np = plotdata[i]->GetN();
      for(int p=0;p<np;p++) {
	plotdata[i]->SetPointEXlow(p,0);
	plotdata[i]->SetPointEXhigh(p,0);
	plotmc[i]->SetPointEXlow(p,0);
	plotmc[i]->SetPointEXhigh(p,0);
	plotdata[i]->SetPointEYhigh(p,plotdata[i]->GetErrorYlow(p));
	plotmc[i]->SetPointEYhigh(p,plotmc[i]->GetErrorYlow(p));
      }
    

      plotdata[i]->SetMarkerSize(1.5);
      plotdata[i]->SetMarkerStyle(20);
      plotmc[i]->SetMarkerSize(1.5);
      plotmc[i]->SetMarkerStyle(20);

      plotdata[i]->SetLineColor(kRed+1);
      plotmc[i]->SetLineColor(kAzure-6);

      plotdata[i]->SetMarkerColor(kRed+1);
      plotmc[i]->SetMarkerColor(kAzure-6);

      if(i==0) {
	legend->AddEntry(plotdata[i],"Data 2011","pl");
	legend->AddEntry(plotmc[i],"MC","pl");
      }

      plotdata[i]->Draw("ape");
      plotmc[i]->Draw("pe");
      legend->Draw();

      TPaveText* text  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
      text->AddText("#sqrt{s} = 7 TeV 2011, L = 4.9 fb^{-1}");
      text->AddText(etabins[e]);
      text->SetBorderSize(0);
      text->SetFillStyle(0);
      text->SetTextAlign(12);
      text->SetTextFont(32);
      text->SetTextSize(0.06);
      text->Draw();

      c1.SaveAs(files[i]+TString(".pdf"));
      c1.SaveAs(files[i]+TString(".png"));
    }
  }

}
