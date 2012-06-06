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

#include "plotUtil.cxx"
#include "collageMuonID_ZZ4L.cxx"

using namespace std;

void plotTPeff() {

  vector<TString> filesData, filesMC;
  filesData.push_back("results/hzzref/2012/newhzzWP_eff_Data8TeV_Kine_Pt10To1000.root");
  filesMC.push_back("results/hzzref/2012/newhzzWP_eff_MC8TeV_Kine_Pt10To1000.root");

  map<int,TString> etabins;
  etabins.insert(make_pair(0,"0<|#eta|<0.8"));
  etabins.insert(make_pair(1,"0.8<|#eta|<1.4442"));
  etabins.insert(make_pair(2,"1.4442<|#eta|<1.566"));
  etabins.insert(make_pair(3,"1.566<|#eta|<2.0"));
  etabins.insert(make_pair(4,"2.0<|#eta|<2.5"));


  for(unsigned int e=0;e<etabins.size();e++) {

    char filename[100];
    sprintf(filename,"eff_etabin%d",e);
    
    vector<TString> files;
    files.push_back(filename);
    
    TString dirdata = TString("eleIDdir/etapt/fit_eff_plots/");
    TString dirmc = TString("eleIDdir/etapt/cnt_eff_plots/");
    char plotname[100];
    sprintf(plotname,"pt_PLOT_abseta_bin%d",e);
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
    
    vector<TGraphAsymmErrors*> plotdata, plotmc;
    for(int i=0;i<(int)filesData.size();++i) {
      TFile *fileData = TFile::Open(filesData[i]);
      TFile *fileMC = TFile::Open(filesMC[i]);

      TCanvas *cdata = (TCanvas*)fileData->Get(dirdata+plot);
      TCanvas *cmc = (TCanvas*)fileMC->Get(dirmc+plot);
      plotdata.push_back((TGraphAsymmErrors*)cdata->FindObject("hxy_fit_eff"));
      plotmc.push_back((TGraphAsymmErrors*)cmc->FindObject("hxy_cnt_eff"));

      c1.cd();
    
      plotdata[i]->SetMinimum(0.0);
      plotdata[i]->SetMaximum(1.2);
      plotdata[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      plotdata[i]->GetYaxis()->SetTitle("efficiency");
      plotdata[i]->GetXaxis()->SetRangeUser(0,80);

      // cosmetics
      for(int b=0;b<plotdata[i]->GetN();b++) {
	// minos give some large upper error when close to 1. Put the symmetric one from lower error
	float yerrlo=plotdata[i]->GetErrorYlow(b);
	plotdata[i]->SetPointEYhigh(b,yerrlo);
      }

      plotdata[i]->SetMarkerSize(1.5);
      plotdata[i]->SetMarkerStyle(20);
      plotmc[i]->SetMarkerSize(1.5);
      plotmc[i]->SetMarkerStyle(20);

      plotdata[i]->SetLineColor(kAzure-6);
      plotmc[i]->SetLineColor(kRed+1);

      plotdata[i]->SetMarkerColor(kAzure-6);
      plotmc[i]->SetMarkerColor(kRed+1);

      if(i==0) {
	legend->AddEntry(plotdata[i],"Data 2012","pl");
	legend->AddEntry(plotmc[i],"MC","pl");
      }

      plotdata[i]->Draw("ape");
      plotmc[i]->Draw("pe");
      legend->Draw();

      TPaveText* text  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
      text->AddText("#sqrt{s} = 8 TeV 2012, L = 1.6 fb^{-1}");
      text->AddText(etabins[e]);
      text->SetBorderSize(0);
      text->SetFillStyle(0);
      text->SetTextAlign(12);
      text->SetTextFont(32);
      text->SetTextSize(0.06);
      text->Draw();

      c1.SaveAs(files[i]+TString(".pdf"));
      c1.SaveAs(files[i]+TString(".png"));

      doRatio(plotdata[i],plotmc[i],files[i],"p_{T} [GeV]",text);

    }
  }

}



void plotTPeffHLT() {

  vector<TString> filesData, filesMC;
  filesData.push_back("Electrons_CBVoigtianPlusChebychevBackground_TnP_Z_DATA_2012_passHLTEle17Leg_Latinos_PtEtaBins.root");

  map<int,TString> etabins;
  etabins.insert(make_pair(0,"0<|#eta|<1.479"));
  etabins.insert(make_pair(1,"1.479<|#eta|<2.5"));

  for(unsigned int e=0;e<etabins.size();e++) {

    char filename[100];
    sprintf(filename,"eff_etabin%d",e);
    
    vector<TString> files;
    files.push_back(filename);
    
    TString dirdata = TString("eleIDdir/etapt/fit_eff_plots/");
    TString dirmc = TString("eleIDdir/etapt/cnt_eff_plots/");
    char plotname[100];
    sprintf(plotname,"pt_PLOT_abseta_bin%d_&_passConvR_pass_&_passID_pass_&_passIP_pass_&_passIso_pass",e);
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
    
    vector<TGraphAsymmErrors*> plotdata, plotmc;
    for(int i=0;i<(int)filesData.size();++i) {
      TFile *fileData = TFile::Open(filesData[i]);

      TCanvas *cdata = (TCanvas*)fileData->Get(dirdata+plot);
      plotdata.push_back((TGraphAsymmErrors*)cdata->FindObject("hxy_fit_eff"));

      c1.cd();
    
      plotdata[i]->SetMinimum(0.0);
      plotdata[i]->SetMaximum(1.2);
      plotdata[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      plotdata[i]->GetYaxis()->SetTitle("efficiency");
      plotdata[i]->GetXaxis()->SetRangeUser(0,80);

      TF1 *fitFunc = paramTurnOnHLT(plotdata[i], 7, 200, 200);
      fitFunc->SetLineColor(100); fitFunc->SetLineWidth(2);

      // cosmetics

//       for(int b=0;b<plotdata[i]->GetN();b++) {
// 	// minos give some large upper error when close to 1. Put the symmetric one from lower error
// 	float yerrlo=plotdata[i]->GetErrorYlow(b);
// 	plotdata[i]->SetPointEYhigh(b,yerrlo);
//       }

      plotdata[i]->SetMarkerSize(1.5);
      plotdata[i]->SetMarkerStyle(20);

      plotdata[i]->SetLineColor(kAzure-6);

      plotdata[i]->SetMarkerColor(kAzure-6);

      if(i==0) {
	legend->AddEntry(plotdata[i],"Data 2012","pl");
      }

      plotdata[i]->Draw("ape");
      fitFunc->Draw("L SAME");
      legend->Draw();

      TPaveText* text  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
      text->AddText("#sqrt{s} = 8 TeV 2012, L = 1.6 fb^{-1}");
      text->AddText(etabins[e]);
      text->SetBorderSize(0);
      text->SetFillStyle(0);
      text->SetTextAlign(12);
      text->SetTextFont(32);
      text->SetTextSize(0.06);
      text->Draw();

      c1.SaveAs(files[i]+TString(".pdf"));
      c1.SaveAs(files[i]+TString(".png"));

      //      doRatio(plotdata[i],plotmc[i],files[i],"p_{T} [GeV]",text);

    }
  }

}
