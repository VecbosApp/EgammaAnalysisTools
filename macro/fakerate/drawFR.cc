#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TLegend.h>

using namespace std;

void drawOneComparison(vector<TH1F*> histos, vector<TString> descr, TString xaxislabel, const char *filename) {

  if(histos.size()>3) {
    cout << "more than 3 histos not implemented." << endl;
    return;
  }

  if(histos.size()!=descr.size()) {
    cout << "description not complete!" << endl;
    return;
  }

  vector<int> colors;
  colors.push_back(kRed+1);
  colors.push_back(kAzure-6);
  colors.push_back(kTeal+3);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetGridx();
  c1->SetGridy();  
  TLegend* legend = new TLegend(0.20, 0.70, 0.43, 0.86);
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);

  for(int i=0;i<(int)histos.size();++i) {
    
    histos[i]->SetMinimum(0);
    histos[i]->SetMaximum(0.2);
    histos[i]->SetMarkerSize(2);
    histos[i]->SetMarkerStyle(20);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetTitle("");
    if(TString(histos[i]->GetName()).Contains("PU")) {
      histos[i]->Fit("pol1","","same",4,35);
      histos[i]->GetFunction("pol1")->SetLineColor(colors[i]);
    }
    histos[i]->GetXaxis()->SetTitle(xaxislabel);
    histos[i]->GetYaxis()->SetTitle("efficiency");

    legend->AddEntry(histos[i],descr[i]);

    if(i==0) { 
      histos[i]->Draw("pe1");
      for(int bin=1;bin<=histos[i]->GetNbinsX();bin++) {
	//cout << "bin i = " << bin << ": " <<  histos[i]->GetBinContent(bin) << " +/- " << histos[i]->GetBinError(bin) << endl;
	cout << "m35_fakeRate[" << bin << "] = " << histos[i]->GetBinContent(bin) << ";" << endl;
      }
      for(int bin=1;bin<=histos[i]->GetNbinsX();bin++) {
	cout << "m35_fakeRate_err[" << bin << "] = " << histos[i]->GetBinError(bin) << ";" << endl;
      }
    }
    else histos[i]->Draw("same pe1");
  }
  legend->Draw();

  TString basename(filename);
  basename.ReplaceAll("Eff","FR");
  c1->SaveAs(basename+TString(".png"));
  c1->SaveAs(basename+TString(".pdf"));

}

void drawOneToOne(vector<TH1F*> set1, vector<TH1F*> set2, const char* desc1, const char* desc2, const char* xaxislabel) {
  if(set1.size()!=set2.size()) {
    cout << "first set and second set of histos have different sizes! ERROR! " << endl;
    return;
  }
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    histos.push_back(set2[i]);
    desc.push_back(TString(desc1));
    desc.push_back(TString(desc2));
    if(set1[i]==0 || set2[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}

void drawOneToTwo(vector<TH1F*> set1, vector<TH1F*> set2, vector<TH1F*> set3, const char* desc1, const char* desc2, const char* desc3, const char* xaxislabel) {
  if(set1.size()!=set2.size() || set1.size()!=set3.size()) {
    cout << "first set and second/third set of histos have different sizes! ERROR! " << endl;
    return;
  }
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    histos.push_back(set2[i]);
    histos.push_back(set3[i]);
    desc.push_back(TString(desc1));
    desc.push_back(TString(desc2));
    desc.push_back(TString(desc3));
    if(set1[i]==0 || set2[i]==0 || set3[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}

void drawIsolations() {
  
  // here isolation only, WP80. Corrected and not
  TFile *file = TFile::Open("fakerates.root");

  // eta
  TH1F *TrgEleIsoWP80EtaHighPt = (TH1F*)file->Get("TrgEleIsoWP80EtaHighPt_Eff");
  TH1F *TrgEleIsoWP80EtaLowPt = (TH1F*)file->Get("TrgEleIsoWP80EtaLowPt_Eff");
  // ---> EA corrected
  TH1F *TrgEleIsoWP80EAEtaHighPt = (TH1F*)file->Get("TrgEleIsoWP80EAEtaHighPt_Eff");
  TH1F *TrgEleIsoWP80EAEtaLowPt = (TH1F*)file->Get("TrgEleIsoWP80EAEtaLowPt_Eff");
  // ---> not EA corrected
  TH1F *TrgEleIsoWP80ZZNoEAEtaHighPt = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAEtaHighPt_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAEtaLowPt = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2, etaSet3;
  etaSet1.push_back(TrgEleIsoWP80EtaHighPt);
  etaSet1.push_back(TrgEleIsoWP80EtaLowPt);
  etaSet2.push_back(TrgEleIsoWP80EAEtaHighPt);
  etaSet2.push_back(TrgEleIsoWP80EAEtaLowPt);
  etaSet3.push_back(TrgEleIsoWP80ZZNoEAEtaHighPt);
  etaSet3.push_back(TrgEleIsoWP80ZZNoEAEtaLowPt);
  
  drawOneToTwo(etaSet1,etaSet3,etaSet2,"PF iso","PF iso, new vetoes","PF iso, new vetoes, EA","#eta");

  // pt
  TH1F *TrgEleIsoWP80PtBarrel1 = (TH1F*)file->Get("TrgEleIsoWP80PtBarrel1_Eff");
  TH1F *TrgEleIsoWP80PtBarrel2 = (TH1F*)file->Get("TrgEleIsoWP80PtBarrel2_Eff");
  TH1F *TrgEleIsoWP80PtEndcap1 = (TH1F*)file->Get("TrgEleIsoWP80PtEndcap1_Eff");
  TH1F *TrgEleIsoWP80PtEndcap2 = (TH1F*)file->Get("TrgEleIsoWP80PtEndcap2_Eff");
  // ---> EA corrected
  TH1F *TrgEleIsoWP80EAPtBarrel1 = (TH1F*)file->Get("TrgEleIsoWP80EAPtBarrel1_Eff");
  TH1F *TrgEleIsoWP80EAPtBarrel2 = (TH1F*)file->Get("TrgEleIsoWP80EAPtBarrel2_Eff");
  TH1F *TrgEleIsoWP80EAPtEndcap1 = (TH1F*)file->Get("TrgEleIsoWP80EAPtEndcap1_Eff");
  TH1F *TrgEleIsoWP80EAPtEndcap2 = (TH1F*)file->Get("TrgEleIsoWP80EAPtEndcap2_Eff");
  // ---> not EA corrected
  TH1F *TrgEleIsoWP80ZZNoEAPtBarrel1 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPtBarrel1_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAPtBarrel2 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPtBarrel2_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAPtEndcap1 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPtEndcap1_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAPtEndcap2 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPtEndcap2_Eff");

  vector<TH1F*> ptSet1, ptSet2, ptSet3;
  ptSet1.push_back(TrgEleIsoWP80PtBarrel1);
  ptSet1.push_back(TrgEleIsoWP80PtBarrel2);
  ptSet1.push_back(TrgEleIsoWP80PtEndcap1);
  ptSet1.push_back(TrgEleIsoWP80PtEndcap2);
  ptSet2.push_back(TrgEleIsoWP80EAPtBarrel1);
  ptSet2.push_back(TrgEleIsoWP80EAPtBarrel2);
  ptSet2.push_back(TrgEleIsoWP80EAPtEndcap1);
  ptSet2.push_back(TrgEleIsoWP80EAPtEndcap2);
  ptSet3.push_back(TrgEleIsoWP80ZZNoEAPtBarrel1);
  ptSet3.push_back(TrgEleIsoWP80ZZNoEAPtBarrel2);
  ptSet3.push_back(TrgEleIsoWP80ZZNoEAPtEndcap1);
  ptSet3.push_back(TrgEleIsoWP80ZZNoEAPtEndcap2);

  drawOneToTwo(ptSet1,ptSet3,ptSet2,"PF iso","PF iso, new vetoes","PF iso, new vetoes, EA","p_{T} [GeV]");

  // PU
  TH1F *TrgEleIsoWP80PUBarrel1 = (TH1F*)file->Get("TrgEleIsoWP80PUBarrel1_Eff");
  TH1F *TrgEleIsoWP80PUBarrel2 = (TH1F*)file->Get("TrgEleIsoWP80PUBarrel2_Eff");
  TH1F *TrgEleIsoWP80PUEndcap1 = (TH1F*)file->Get("TrgEleIsoWP80PUEndcap1_Eff");
  TH1F *TrgEleIsoWP80PUEndcap2 = (TH1F*)file->Get("TrgEleIsoWP80PUEndcap2_Eff");
  // ---> EA corrected
  TH1F *TrgEleIsoWP80EAPUBarrel1 = (TH1F*)file->Get("TrgEleIsoWP80EAPUBarrel1_Eff");
  TH1F *TrgEleIsoWP80EAPUBarrel2 = (TH1F*)file->Get("TrgEleIsoWP80EAPUBarrel2_Eff");
  TH1F *TrgEleIsoWP80EAPUEndcap1 = (TH1F*)file->Get("TrgEleIsoWP80EAPUEndcap1_Eff");
  TH1F *TrgEleIsoWP80EAPUEndcap2 = (TH1F*)file->Get("TrgEleIsoWP80EAPUEndcap2_Eff");
  // ---> not EA corrected
  TH1F *TrgEleIsoWP80ZZNoEAPUBarrel1 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPUBarrel1_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAPUBarrel2 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPUBarrel2_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAPUEndcap1 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPUEndcap1_Eff");
  TH1F *TrgEleIsoWP80ZZNoEAPUEndcap2 = (TH1F*)file->Get("TrgEleIsoWP80ZZNoEAPUEndcap2_Eff");

  vector<TH1F*> puSet1, puSet2, puSet3;
  puSet1.push_back(TrgEleIsoWP80PUBarrel1);
  puSet1.push_back(TrgEleIsoWP80PUBarrel2);
  puSet1.push_back(TrgEleIsoWP80PUEndcap1);
  puSet1.push_back(TrgEleIsoWP80PUEndcap2);
  puSet2.push_back(TrgEleIsoWP80EAPUBarrel1);
  puSet2.push_back(TrgEleIsoWP80EAPUBarrel2);
  puSet2.push_back(TrgEleIsoWP80EAPUEndcap1);
  puSet2.push_back(TrgEleIsoWP80EAPUEndcap2);
  puSet3.push_back(TrgEleIsoWP80ZZNoEAPUBarrel1);
  puSet3.push_back(TrgEleIsoWP80ZZNoEAPUBarrel2);
  puSet3.push_back(TrgEleIsoWP80ZZNoEAPUEndcap1);
  puSet3.push_back(TrgEleIsoWP80ZZNoEAPUEndcap2);

  drawOneToTwo(puSet1,puSet3,puSet2,"PF iso","PF iso, new vetoes","PF iso, new vetoes, EA","# vertices");

}

void drawIdsBiased() {

  TFile *file = TFile::Open("fakerates_trigger.root");

  // eta
  TH1F *TrgEleWP80EtaHighPt = (TH1F*)file->Get("TrgEleWP80EtaHighPt_Eff");
  TH1F *TrgEleWP80EtaLowPt = (TH1F*)file->Get("TrgEleWP80EtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *TrgElenewWPHWWEtaHighPt = (TH1F*)file->Get("TrgElenewWPHWWEtaHighPt_Eff");
  TH1F *TrgElenewWPHWWEtaLowPt = (TH1F*)file->Get("TrgElenewWPHWWEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(TrgEleWP80EtaHighPt);
  etaSet1.push_back(TrgEleWP80EtaLowPt);
  etaSet2.push_back(TrgElenewWPHWWEtaHighPt);
  etaSet2.push_back(TrgElenewWPHWWEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H #rightarrow WW 2011","H #rightarrow WW 2012 (same-eff)","#eta");

  // pt
  TH1F *TrgEleWP80PtBarrel1 = (TH1F*)file->Get("TrgEleWP80PtBarrel1_Eff");
  TH1F *TrgEleWP80PtBarrel2 = (TH1F*)file->Get("TrgEleWP80PtBarrel2_Eff");
  TH1F *TrgEleWP80PtEndcap1 = (TH1F*)file->Get("TrgEleWP80PtEndcap1_Eff");
  TH1F *TrgEleWP80PtEndcap2 = (TH1F*)file->Get("TrgEleWP80PtEndcap2_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *TrgElenewWPHWWPtBarrel1 = (TH1F*)file->Get("TrgElenewWPHWWPtBarrel1_Eff");
  TH1F *TrgElenewWPHWWPtBarrel2 = (TH1F*)file->Get("TrgElenewWPHWWPtBarrel2_Eff");
  TH1F *TrgElenewWPHWWPtEndcap1 = (TH1F*)file->Get("TrgElenewWPHWWPtEndcap1_Eff");
  TH1F *TrgElenewWPHWWPtEndcap2 = (TH1F*)file->Get("TrgElenewWPHWWPtEndcap2_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(TrgEleWP80PtBarrel1);
  ptSet1.push_back(TrgEleWP80PtBarrel2);
  ptSet1.push_back(TrgEleWP80PtEndcap1);
  ptSet1.push_back(TrgEleWP80PtEndcap2);
  ptSet2.push_back(TrgElenewWPHWWPtBarrel1);
  ptSet2.push_back(TrgElenewWPHWWPtBarrel2);
  ptSet2.push_back(TrgElenewWPHWWPtEndcap1);
  ptSet2.push_back(TrgElenewWPHWWPtEndcap2);

  drawOneToOne(ptSet1,ptSet2,"H #rightarrow WW 2011","H #rightarrow WW 2012 (same-eff)","p_{T} [GeV]");

  // PU
  TH1F *TrgEleWP80PUBarrel1 = (TH1F*)file->Get("TrgEleWP80PUBarrel1_Eff");
  TH1F *TrgEleWP80PUBarrel2 = (TH1F*)file->Get("TrgEleWP80PUBarrel2_Eff");
  TH1F *TrgEleWP80PUEndcap1 = (TH1F*)file->Get("TrgEleWP80PUEndcap1_Eff");
  TH1F *TrgEleWP80PUEndcap2 = (TH1F*)file->Get("TrgEleWP80PUEndcap2_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *TrgElenewWPHWWPUBarrel1 = (TH1F*)file->Get("TrgElenewWPHWWPUBarrel1_Eff");
  TH1F *TrgElenewWPHWWPUBarrel2 = (TH1F*)file->Get("TrgElenewWPHWWPUBarrel2_Eff");
  TH1F *TrgElenewWPHWWPUEndcap1 = (TH1F*)file->Get("TrgElenewWPHWWPUEndcap1_Eff");
  TH1F *TrgElenewWPHWWPUEndcap2 = (TH1F*)file->Get("TrgElenewWPHWWPUEndcap2_Eff");

  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(TrgEleWP80PUBarrel1);
  puSet1.push_back(TrgEleWP80PUBarrel2);
  puSet1.push_back(TrgEleWP80PUEndcap1);
  puSet1.push_back(TrgEleWP80PUEndcap2);
  puSet2.push_back(TrgElenewWPHWWPUBarrel1);
  puSet2.push_back(TrgElenewWPHWWPUBarrel2);
  puSet2.push_back(TrgElenewWPHWWPUEndcap1);
  puSet2.push_back(TrgElenewWPHWWPUEndcap2);

  drawOneToOne(puSet1,puSet2,"H #rightarrow WW 2011","H #rightarrow WW 2012 (same-eff)","# vertices");


}


void drawIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zee1fake.root");

  // eta
  TH1F *NoTrgEleCiCMediumEtaHighPt = (TH1F*)file->Get("NoTrgEleCiCMediumEtaHighPt_Eff");
  TH1F *NoTrgEleCiCMediumEtaLowPt = (TH1F*)file->Get("NoTrgEleCiCMediumEtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElenewWPHZZEtaHighPt = (TH1F*)file->Get("NoTrgElenewWPHZZEtaHighPt_Eff");
  TH1F *NoTrgElenewWPHZZEtaLowPt = (TH1F*)file->Get("NoTrgElenewWPHZZEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(NoTrgEleCiCMediumEtaHighPt);
  etaSet1.push_back(NoTrgEleCiCMediumEtaLowPt);
  etaSet2.push_back(NoTrgElenewWPHZZEtaHighPt);
  etaSet2.push_back(NoTrgElenewWPHZZEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H#rightarrow ZZ CiC","H#rightarrow ZZ 2012","#eta");

  // pt
  TH1F *NoTrgEleCiCMediumPtBarrel = (TH1F*)file->Get("NoTrgEleCiCMediumPtBarrel_Eff");
  TH1F *NoTrgEleCiCMediumPtEndcap = (TH1F*)file->Get("NoTrgEleCiCMediumPtEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElenewWPHZZPtBarrel = (TH1F*)file->Get("NoTrgElenewWPHZZPtBarrel_Eff");
  TH1F *NoTrgElenewWPHZZPtEndcap = (TH1F*)file->Get("NoTrgElenewWPHZZPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(NoTrgEleCiCMediumPtBarrel);
  ptSet1.push_back(NoTrgEleCiCMediumPtEndcap);
  ptSet2.push_back(NoTrgElenewWPHZZPtBarrel);
  ptSet2.push_back(NoTrgElenewWPHZZPtEndcap);

  drawOneToOne(ptSet1,ptSet2,"H#rightarrow ZZ CiC","H#rightarrow ZZ 2012","p_{T} [GeV]");

  // PU
  TH1F *NoTrgEleCiCMediumPUBarrel = (TH1F*)file->Get("NoTrgEleCiCMediumPUBarrel_Eff");
  TH1F *NoTrgEleCiCMediumPUEndcap = (TH1F*)file->Get("NoTrgEleCiCMediumPUEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElenewWPHZZPUBarrel = (TH1F*)file->Get("NoTrgElenewWPHZZPUBarrel_Eff");
  TH1F *NoTrgElenewWPHZZPUEndcap = (TH1F*)file->Get("NoTrgElenewWPHZZPUEndcap_Eff");

  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(NoTrgEleCiCMediumPUBarrel);
  puSet1.push_back(NoTrgEleCiCMediumPUEndcap);
  puSet2.push_back(NoTrgElenewWPHZZPUBarrel);
  puSet2.push_back(NoTrgElenewWPHZZPUEndcap);

  drawOneToOne(puSet1,puSet2,"H#rightarrow ZZ CiC","H#rightarrow ZZ 2012","# vertices");


}

void drawIdsUnbiasedNewWPs() {

  TFile *file = TFile::Open("fakerates_zee1fake.root");

  // eta
  TH1F *NoTrgElehzzMvaLooseEtaHighPt = (TH1F*)file->Get("NoTrgElehzzMvaLooseEtaHighPt_Eff");
  TH1F *NoTrgElehzzMvaLooseEtaLowPt = (TH1F*)file->Get("NoTrgElehzzMvaLooseEtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElehzzMvaTightEtaHighPt = (TH1F*)file->Get("NoTrgElehzzMvaTightEtaHighPt_Eff");
  TH1F *NoTrgElehzzMvaTightEtaLowPt = (TH1F*)file->Get("NoTrgElehzzMvaTightEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(NoTrgElehzzMvaLooseEtaHighPt);
  etaSet1.push_back(NoTrgElehzzMvaLooseEtaLowPt);
  etaSet2.push_back(NoTrgElehzzMvaTightEtaHighPt);
  etaSet2.push_back(NoTrgElehzzMvaTightEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","#eta");

  // pt
  TH1F *NoTrgElehzzMvaLoosePtBarrel = (TH1F*)file->Get("NoTrgElehzzMvaLoosePtBarrel_Eff");
  TH1F *NoTrgElehzzMvaLoosePtEndcap = (TH1F*)file->Get("NoTrgElehzzMvaLoosePtEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElehzzMvaTightPtBarrel = (TH1F*)file->Get("NoTrgElehzzMvaTightPtBarrel_Eff");
  TH1F *NoTrgElehzzMvaTightPtEndcap = (TH1F*)file->Get("NoTrgElehzzMvaTightPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(NoTrgElehzzMvaLoosePtBarrel);
  ptSet1.push_back(NoTrgElehzzMvaLoosePtEndcap);
  ptSet2.push_back(NoTrgElehzzMvaTightPtBarrel);
  ptSet2.push_back(NoTrgElehzzMvaTightPtEndcap);

  drawOneToOne(ptSet1,ptSet2,"H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","p_{T} [GeV]");

  // PU
  TH1F *NoTrgElehzzMvaLoosePUBarrel = (TH1F*)file->Get("NoTrgElehzzMvaLoosePUBarrel_Eff");
  TH1F *NoTrgElehzzMvaLoosePUEndcap = (TH1F*)file->Get("NoTrgElehzzMvaLoosePUEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElehzzMvaTightPUBarrel = (TH1F*)file->Get("NoTrgElehzzMvaTightPUBarrel_Eff");
  TH1F *NoTrgElehzzMvaTightPUEndcap = (TH1F*)file->Get("NoTrgElehzzMvaTightPUEndcap_Eff");

  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(NoTrgElehzzMvaLoosePUBarrel);
  puSet1.push_back(NoTrgElehzzMvaLoosePUEndcap);
  puSet2.push_back(NoTrgElehzzMvaTightPUBarrel);
  puSet2.push_back(NoTrgElehzzMvaTightPUEndcap);

  drawOneToOne(puSet1,puSet2,"H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","# vertices");


}


void drawIdsUnbiasedIsoBins() {

  TFile *file = TFile::Open("fakerates_zee1fake.root");

  // pt
  TH1F *NoTrgEleCiCMediumIso1PtBarrel = (TH1F*)file->Get("NoTrgEleCiCMediumIso1PtBarrel_Eff");
  TH1F *NoTrgEleCiCMediumIso2PtBarrel = (TH1F*)file->Get("NoTrgEleCiCMediumIso2PtBarrel_Eff");
  TH1F *NoTrgEleCiCMediumIso1PtEndcap = (TH1F*)file->Get("NoTrgEleCiCMediumIso1PtEndcap_Eff");
  TH1F *NoTrgEleCiCMediumIso2PtEndcap = (TH1F*)file->Get("NoTrgEleCiCMediumIso2PtEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgElenewWPHZZIso1PtBarrel = (TH1F*)file->Get("NoTrgElenewWPHZZIso1PtBarrel_Eff");
  TH1F *NoTrgElenewWPHZZIso2PtBarrel = (TH1F*)file->Get("NoTrgElenewWPHZZIso2PtBarrel_Eff");
  TH1F *NoTrgElenewWPHZZIso1PtEndcap = (TH1F*)file->Get("NoTrgElenewWPHZZIso1PtEndcap_Eff");
  TH1F *NoTrgElenewWPHZZIso2PtEndcap = (TH1F*)file->Get("NoTrgElenewWPHZZIso2PtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(NoTrgEleCiCMediumIso1PtBarrel);
  ptSet1.push_back(NoTrgEleCiCMediumIso2PtBarrel);
  ptSet1.push_back(NoTrgEleCiCMediumIso1PtEndcap);
  ptSet1.push_back(NoTrgEleCiCMediumIso2PtEndcap);
  ptSet2.push_back(NoTrgElenewWPHZZIso1PtBarrel);
  ptSet2.push_back(NoTrgElenewWPHZZIso2PtBarrel);
  ptSet2.push_back(NoTrgElenewWPHZZIso1PtEndcap);
  ptSet2.push_back(NoTrgElenewWPHZZIso2PtEndcap);

  drawOneToOne(ptSet1,ptSet2,"H#rightarrow ZZ CiC","H#rightarrow ZZ 2012","p_{T} [GeV]");

}
