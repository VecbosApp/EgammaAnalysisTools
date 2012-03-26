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
    histos[i]->SetMaximum(0.4);
    histos[i]->SetMarkerSize(2);
    histos[i]->SetMarkerStyle(20);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetTitle("");
    if(TString(histos[i]->GetName()).Contains("PU")) {
      histos[i]->Fit("pol1","","same",4,25);
      histos[i]->GetFunction("pol1")->SetLineColor(colors[i]);
    }
    histos[i]->GetXaxis()->SetTitle(xaxislabel);
    histos[i]->GetYaxis()->SetTitle("efficiency");

    legend->AddEntry(histos[i],descr[i]);

    if(i==0) histos[i]->Draw("pe1");
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
  TH1F *BdtHWWIsoWP80EtaHighPt = (TH1F*)file->Get("BdtHWWIsoWP80EtaHighPt_Eff");
  TH1F *BdtHWWIsoWP80EtaLowPt = (TH1F*)file->Get("BdtHWWIsoWP80EtaLowPt_Eff");
  // ---> EA corrected
  TH1F *BdtHWWIsoWP80EAEtaHighPt = (TH1F*)file->Get("BdtHWWIsoWP80EAEtaHighPt_Eff");
  TH1F *BdtHWWIsoWP80EAEtaLowPt = (TH1F*)file->Get("BdtHWWIsoWP80EAEtaLowPt_Eff");
  // ---> not EA corrected
  TH1F *BdtHWWIsoWP80ZZNoEAEtaHighPt = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAEtaHighPt_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAEtaLowPt = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2, etaSet3;
  etaSet1.push_back(BdtHWWIsoWP80EtaHighPt);
  etaSet1.push_back(BdtHWWIsoWP80EtaLowPt);
  etaSet2.push_back(BdtHWWIsoWP80EAEtaHighPt);
  etaSet2.push_back(BdtHWWIsoWP80EAEtaLowPt);
  etaSet3.push_back(BdtHWWIsoWP80ZZNoEAEtaHighPt);
  etaSet3.push_back(BdtHWWIsoWP80ZZNoEAEtaLowPt);
  
  drawOneToTwo(etaSet1,etaSet3,etaSet2,"PF iso","PF iso, new vetoes","PF iso, new vetoes, EA","#eta");

  // pt
  TH1F *BdtHWWIsoWP80PtBarrel1 = (TH1F*)file->Get("BdtHWWIsoWP80PtBarrel1_Eff");
  TH1F *BdtHWWIsoWP80PtBarrel2 = (TH1F*)file->Get("BdtHWWIsoWP80PtBarrel2_Eff");
  TH1F *BdtHWWIsoWP80PtEndcap1 = (TH1F*)file->Get("BdtHWWIsoWP80PtEndcap1_Eff");
  TH1F *BdtHWWIsoWP80PtEndcap2 = (TH1F*)file->Get("BdtHWWIsoWP80PtEndcap2_Eff");
  // ---> EA corrected
  TH1F *BdtHWWIsoWP80EAPtBarrel1 = (TH1F*)file->Get("BdtHWWIsoWP80EAPtBarrel1_Eff");
  TH1F *BdtHWWIsoWP80EAPtBarrel2 = (TH1F*)file->Get("BdtHWWIsoWP80EAPtBarrel2_Eff");
  TH1F *BdtHWWIsoWP80EAPtEndcap1 = (TH1F*)file->Get("BdtHWWIsoWP80EAPtEndcap1_Eff");
  TH1F *BdtHWWIsoWP80EAPtEndcap2 = (TH1F*)file->Get("BdtHWWIsoWP80EAPtEndcap2_Eff");
  // ---> not EA corrected
  TH1F *BdtHWWIsoWP80ZZNoEAPtBarrel1 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPtBarrel1_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAPtBarrel2 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPtBarrel2_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAPtEndcap1 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPtEndcap1_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAPtEndcap2 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPtEndcap2_Eff");

  vector<TH1F*> ptSet1, ptSet2, ptSet3;
  ptSet1.push_back(BdtHWWIsoWP80PtBarrel1);
  ptSet1.push_back(BdtHWWIsoWP80PtBarrel2);
  ptSet1.push_back(BdtHWWIsoWP80PtEndcap1);
  ptSet1.push_back(BdtHWWIsoWP80PtEndcap2);
  ptSet2.push_back(BdtHWWIsoWP80EAPtBarrel1);
  ptSet2.push_back(BdtHWWIsoWP80EAPtBarrel2);
  ptSet2.push_back(BdtHWWIsoWP80EAPtEndcap1);
  ptSet2.push_back(BdtHWWIsoWP80EAPtEndcap2);
  ptSet3.push_back(BdtHWWIsoWP80ZZNoEAPtBarrel1);
  ptSet3.push_back(BdtHWWIsoWP80ZZNoEAPtBarrel2);
  ptSet3.push_back(BdtHWWIsoWP80ZZNoEAPtEndcap1);
  ptSet3.push_back(BdtHWWIsoWP80ZZNoEAPtEndcap2);

  drawOneToTwo(ptSet1,ptSet3,ptSet2,"PF iso","PF iso, new vetoes","PF iso, new vetoes, EA","p_{T} [GeV]");

  // PU
  TH1F *BdtHWWIsoWP80PUBarrel1 = (TH1F*)file->Get("BdtHWWIsoWP80PUBarrel1_Eff");
  TH1F *BdtHWWIsoWP80PUBarrel2 = (TH1F*)file->Get("BdtHWWIsoWP80PUBarrel2_Eff");
  TH1F *BdtHWWIsoWP80PUEndcap1 = (TH1F*)file->Get("BdtHWWIsoWP80PUEndcap1_Eff");
  TH1F *BdtHWWIsoWP80PUEndcap2 = (TH1F*)file->Get("BdtHWWIsoWP80PUEndcap2_Eff");
  // ---> EA corrected
  TH1F *BdtHWWIsoWP80EAPUBarrel1 = (TH1F*)file->Get("BdtHWWIsoWP80EAPUBarrel1_Eff");
  TH1F *BdtHWWIsoWP80EAPUBarrel2 = (TH1F*)file->Get("BdtHWWIsoWP80EAPUBarrel2_Eff");
  TH1F *BdtHWWIsoWP80EAPUEndcap1 = (TH1F*)file->Get("BdtHWWIsoWP80EAPUEndcap1_Eff");
  TH1F *BdtHWWIsoWP80EAPUEndcap2 = (TH1F*)file->Get("BdtHWWIsoWP80EAPUEndcap2_Eff");
  // ---> not EA corrected
  TH1F *BdtHWWIsoWP80ZZNoEAPUBarrel1 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPUBarrel1_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAPUBarrel2 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPUBarrel2_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAPUEndcap1 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPUEndcap1_Eff");
  TH1F *BdtHWWIsoWP80ZZNoEAPUEndcap2 = (TH1F*)file->Get("BdtHWWIsoWP80ZZNoEAPUEndcap2_Eff");

  vector<TH1F*> puSet1, puSet2, puSet3;
  puSet1.push_back(BdtHWWIsoWP80PUBarrel1);
  puSet1.push_back(BdtHWWIsoWP80PUBarrel2);
  puSet1.push_back(BdtHWWIsoWP80PUEndcap1);
  puSet1.push_back(BdtHWWIsoWP80PUEndcap2);
  puSet2.push_back(BdtHWWIsoWP80EAPUBarrel1);
  puSet2.push_back(BdtHWWIsoWP80EAPUBarrel2);
  puSet2.push_back(BdtHWWIsoWP80EAPUEndcap1);
  puSet2.push_back(BdtHWWIsoWP80EAPUEndcap2);
  puSet3.push_back(BdtHWWIsoWP80ZZNoEAPUBarrel1);
  puSet3.push_back(BdtHWWIsoWP80ZZNoEAPUBarrel2);
  puSet3.push_back(BdtHWWIsoWP80ZZNoEAPUEndcap1);
  puSet3.push_back(BdtHWWIsoWP80ZZNoEAPUEndcap2);

  drawOneToTwo(puSet1,puSet3,puSet2,"PF iso","PF iso, new vetoes","PF iso, new vetoes, EA","# vertices");

}

void drawIdsBiased() {

  TFile *file = TFile::Open("fakerates_trigger.root");

  // eta
  TH1F *BdtHWWWP80EtaHighPt = (TH1F*)file->Get("BdtHWWWP80EtaHighPt_Eff");
  TH1F *BdtHWWWP80EtaLowPt = (TH1F*)file->Get("BdtHWWWP80EtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHWWnewWP70x80EtaHighPt = (TH1F*)file->Get("BdtHWWnewWP70x80EtaHighPt_Eff");
  TH1F *BdtHWWnewWP70x80EtaLowPt = (TH1F*)file->Get("BdtHWWnewWP70x80EtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(BdtHWWWP80EtaHighPt);
  etaSet1.push_back(BdtHWWWP80EtaLowPt);
  etaSet2.push_back(BdtHWWnewWP70x80EtaHighPt);
  etaSet2.push_back(BdtHWWnewWP70x80EtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H #rightarrow WW 2011","H #rightarrow WW 2012 (same-eff)","#eta");

  // pt
  TH1F *BdtHWWWP80PtBarrel1 = (TH1F*)file->Get("BdtHWWWP80PtBarrel1_Eff");
  TH1F *BdtHWWWP80PtBarrel2 = (TH1F*)file->Get("BdtHWWWP80PtBarrel2_Eff");
  TH1F *BdtHWWWP80PtEndcap1 = (TH1F*)file->Get("BdtHWWWP80PtEndcap1_Eff");
  TH1F *BdtHWWWP80PtEndcap2 = (TH1F*)file->Get("BdtHWWWP80PtEndcap2_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHWWnewWP70x80PtBarrel1 = (TH1F*)file->Get("BdtHWWnewWP70x80PtBarrel1_Eff");
  TH1F *BdtHWWnewWP70x80PtBarrel2 = (TH1F*)file->Get("BdtHWWnewWP70x80PtBarrel2_Eff");
  TH1F *BdtHWWnewWP70x80PtEndcap1 = (TH1F*)file->Get("BdtHWWnewWP70x80PtEndcap1_Eff");
  TH1F *BdtHWWnewWP70x80PtEndcap2 = (TH1F*)file->Get("BdtHWWnewWP70x80PtEndcap2_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(BdtHWWWP80PtBarrel1);
  ptSet1.push_back(BdtHWWWP80PtBarrel2);
  ptSet1.push_back(BdtHWWWP80PtEndcap1);
  ptSet1.push_back(BdtHWWWP80PtEndcap2);
  ptSet2.push_back(BdtHWWnewWP70x80PtBarrel1);
  ptSet2.push_back(BdtHWWnewWP70x80PtBarrel2);
  ptSet2.push_back(BdtHWWnewWP70x80PtEndcap1);
  ptSet2.push_back(BdtHWWnewWP70x80PtEndcap2);

  drawOneToOne(ptSet1,ptSet2,"H #rightarrow WW 2011","H #rightarrow WW 2012 (same-eff)","p_{T} [GeV]");

  // PU
  TH1F *BdtHWWWP80PUBarrel1 = (TH1F*)file->Get("BdtHWWWP80PUBarrel1_Eff");
  TH1F *BdtHWWWP80PUBarrel2 = (TH1F*)file->Get("BdtHWWWP80PUBarrel2_Eff");
  TH1F *BdtHWWWP80PUEndcap1 = (TH1F*)file->Get("BdtHWWWP80PUEndcap1_Eff");
  TH1F *BdtHWWWP80PUEndcap2 = (TH1F*)file->Get("BdtHWWWP80PUEndcap2_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHWWnewWP70x80PUBarrel1 = (TH1F*)file->Get("BdtHWWnewWP70x80PUBarrel1_Eff");
  TH1F *BdtHWWnewWP70x80PUBarrel2 = (TH1F*)file->Get("BdtHWWnewWP70x80PUBarrel2_Eff");
  TH1F *BdtHWWnewWP70x80PUEndcap1 = (TH1F*)file->Get("BdtHWWnewWP70x80PUEndcap1_Eff");
  TH1F *BdtHWWnewWP70x80PUEndcap2 = (TH1F*)file->Get("BdtHWWnewWP70x80PUEndcap2_Eff");

  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(BdtHWWWP80PUBarrel1);
  puSet1.push_back(BdtHWWWP80PUBarrel2);
  puSet1.push_back(BdtHWWWP80PUEndcap1);
  puSet1.push_back(BdtHWWWP80PUEndcap2);
  puSet2.push_back(BdtHWWnewWP70x80PUBarrel1);
  puSet2.push_back(BdtHWWnewWP70x80PUBarrel2);
  puSet2.push_back(BdtHWWnewWP70x80PUEndcap1);
  puSet2.push_back(BdtHWWnewWP70x80PUEndcap2);

  drawOneToOne(puSet1,puSet2,"H #rightarrow WW 2011","H #rightarrow WW 2012 (same-eff)","# vertices");


}


void drawIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zee1fake.root");

  // eta
  TH1F *BdtHZZWP80EtaHighPt = (TH1F*)file->Get("BdtHZZWP80EtaHighPt_Eff");
  TH1F *BdtHZZWP80EtaLowPt = (TH1F*)file->Get("BdtHZZWP80EtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHZZWP80chEtaHighPt = (TH1F*)file->Get("BdtHZZWP80chEtaHighPt_Eff");
  TH1F *BdtHZZWP80chEtaLowPt = (TH1F*)file->Get("BdtHZZWP80chEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(BdtHZZWP80EtaHighPt);
  etaSet1.push_back(BdtHZZWP80EtaLowPt);
  etaSet2.push_back(BdtHZZWP80chEtaHighPt);
  etaSet2.push_back(BdtHZZWP80chEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"non-trigger WP80 full iso","non-trigger WP80 ch.iso","#eta");

  // pt
  TH1F *BdtHZZWP80PtBarrel = (TH1F*)file->Get("BdtHZZWP80PtBarrel_Eff");
  TH1F *BdtHZZWP80PtEndcap = (TH1F*)file->Get("BdtHZZWP80PtEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHZZWP80chPtBarrel = (TH1F*)file->Get("BdtHZZWP80chPtBarrel_Eff");
  TH1F *BdtHZZWP80chPtEndcap = (TH1F*)file->Get("BdtHZZWP80chPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(BdtHZZWP80PtBarrel);
  ptSet1.push_back(BdtHZZWP80PtEndcap);
  ptSet2.push_back(BdtHZZWP80chPtBarrel);
  ptSet2.push_back(BdtHZZWP80chPtEndcap);

  drawOneToOne(ptSet1,ptSet2,"non-trigger WP80 full iso","non-trigger WP80 ch.iso","p_{T} [GeV]");

  // PU
  TH1F *BdtHZZWP80PUBarrel = (TH1F*)file->Get("BdtHZZWP80PUBarrel_Eff");
  TH1F *BdtHZZWP80PUEndcap = (TH1F*)file->Get("BdtHZZWP80PUEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHZZWP80chPUBarrel = (TH1F*)file->Get("BdtHZZWP80chPUBarrel_Eff");
  TH1F *BdtHZZWP80chPUEndcap = (TH1F*)file->Get("BdtHZZWP80chPUEndcap_Eff");

  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(BdtHZZWP80PUBarrel);
  puSet1.push_back(BdtHZZWP80PUEndcap);
  puSet2.push_back(BdtHZZWP80chPUBarrel);
  puSet2.push_back(BdtHZZWP80chPUEndcap);

  drawOneToOne(puSet1,puSet2,"non-trigger WP80 full iso","non-trigger WP80 ch.iso","# vertices");


}
