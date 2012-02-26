#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
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

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  TLegend* legend = new TLegend(0.20, 0.70, 0.43, 0.86);
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);

  for(int i=0;i<(int)histos.size();++i) {
    
    histos[i]->SetMinimum(0);
    histos[i]->SetMaximum(0.5);
    histos[i]->SetMarkerSize(2);
    histos[i]->SetMarkerStyle(20);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetTitle("");
    histos[i]->GetXaxis()->SetTitle(xaxislabel);
    histos[i]->GetYaxis()->SetTitle("efficiency");

    legend->AddEntry(histos[i],descr[i]);

    if(i==0) histos[i]->Draw("pe1");
    else histos[i]->Draw("same pe1");
  }
  legend->Draw();

  TString basename(filename);
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

void drawIsolations() {
  
  // here isolation only, WP80. Corrected and not
  TFile *file = TFile::Open("fakerates.root");

  // eta
  TH1F *BdtHWWIdIsoWP80EtaHighPt = (TH1F*)file->Get("BdtHWWIdIsoWP80EtaHighPt_Eff");
  TH1F *BdtHWWIdIsoWP80EtaLowPt = (TH1F*)file->Get("BdtHWWIdIsoWP80EtaLowPt_Eff");
  // ---> EA corrected
  TH1F *BdtHWWIdIsoWP80EAEtaHighPt = (TH1F*)file->Get("BdtHWWIdIsoWP80EAEtaHighPt_Eff");
  TH1F *BdtHWWIdIsoWP80EAEtaLowPt = (TH1F*)file->Get("BdtHWWIdIsoWP80EAEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(BdtHWWIdIsoWP80EtaHighPt);
  etaSet1.push_back(BdtHWWIdIsoWP80EtaLowPt);
  etaSet2.push_back(BdtHWWIdIsoWP80EAEtaHighPt);
  etaSet2.push_back(BdtHWWIdIsoWP80EAEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"PF isolation","PF isolation, EA corr","#eta");

  // pt
  TH1F *BdtHWWIdIsoWP80PtBarrel1 = (TH1F*)file->Get("BdtHWWIdIsoWP80PtBarrel1_Eff");
  TH1F *BdtHWWIdIsoWP80PtBarrel2 = (TH1F*)file->Get("BdtHWWIdIsoWP80PtBarrel2_Eff");
  TH1F *BdtHWWIdIsoWP80PtEndcap1 = (TH1F*)file->Get("BdtHWWIdIsoWP80PtEndcap1_Eff");
  TH1F *BdtHWWIdIsoWP80PtEndcap2 = (TH1F*)file->Get("BdtHWWIdIsoWP80PtEndcap2_Eff");
  // ---> EA corrected
  TH1F *BdtHWWIdIsoWP80EAPtBarrel1 = (TH1F*)file->Get("BdtHWWIdIsoWP80EAPtBarrel1_Eff");
  TH1F *BdtHWWIdIsoWP80EAPtBarrel2 = (TH1F*)file->Get("BdtHWWIdIsoWP80EAPtBarrel2_Eff");
  TH1F *BdtHWWIdIsoWP80EAPtEndcap1 = (TH1F*)file->Get("BdtHWWIdIsoWP80EAPtEndcap1_Eff");
  TH1F *BdtHWWIdIsoWP80EAPtEndcap2 = (TH1F*)file->Get("BdtHWWIdIsoWP80EAPtEndcap2_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(BdtHWWIdIsoWP80PtBarrel1);
  ptSet1.push_back(BdtHWWIdIsoWP80PtBarrel2);
  ptSet1.push_back(BdtHWWIdIsoWP80PtEndcap1);
  ptSet1.push_back(BdtHWWIdIsoWP80PtEndcap2);
  ptSet2.push_back(BdtHWWIdIsoWP80EAPtBarrel1);
  ptSet2.push_back(BdtHWWIdIsoWP80EAPtBarrel2);
  ptSet2.push_back(BdtHWWIdIsoWP80EAPtEndcap1);
  ptSet2.push_back(BdtHWWIdIsoWP80EAPtEndcap2);

  drawOneToOne(ptSet1,ptSet2,"PF isolation","PF isolation, EA corr","p_{T} [GeV]");

  // PU
  TH1F *BdtHWWIdIsoWP80PU = (TH1F*)file->Get("BdtHWWIdIsoWP80PU_Eff");
  // ---> EA corrected
  TH1F *BdtHWWIdIsoWP80EAPU = (TH1F*)file->Get("BdtHWWIdIsoWP80EAPU_Eff");
  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(BdtHWWIdIsoWP80PU);
  puSet2.push_back(BdtHWWIdIsoWP80EAPU);

  drawOneToOne(puSet1,puSet2,"PF isolation","PF isolation, EA corr","# vertices");

}

void drawIds() {

  // here isolation only, WP80. Corrected and not
  TFile *file = TFile::Open("fakerates.root");

  // eta
  TH1F *BdtHWWIdWP80EtaHighPt = (TH1F*)file->Get("BdtHWWIdWP80EtaHighPt_Eff");
  TH1F *BdtHWWIdWP80EtaLowPt = (TH1F*)file->Get("BdtHWWIdWP80EtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHZZIdWP80EAEtaHighPt = (TH1F*)file->Get("BdtHZZIdWP80EAEtaHighPt_Eff");
  TH1F *BdtHZZIdWP80EAEtaLowPt = (TH1F*)file->Get("BdtHZZIdWP80EAEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(BdtHWWIdWP80EtaHighPt);
  etaSet1.push_back(BdtHWWIdWP80EtaLowPt);
  etaSet2.push_back(BdtHZZIdWP80EAEtaHighPt);
  etaSet2.push_back(BdtHZZIdWP80EAEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H #rightarrow WW 2011","H #rightarrow ZZ 2012","#eta");

  // pt
  TH1F *BdtHWWIdWP80PtBarrel1 = (TH1F*)file->Get("BdtHWWIdWP80PtBarrel1_Eff");
  TH1F *BdtHWWIdWP80PtBarrel2 = (TH1F*)file->Get("BdtHWWIdWP80PtBarrel2_Eff");
  TH1F *BdtHWWIdWP80PtEndcap1 = (TH1F*)file->Get("BdtHWWIdWP80PtEndcap1_Eff");
  TH1F *BdtHWWIdWP80PtEndcap2 = (TH1F*)file->Get("BdtHWWIdWP80PtEndcap2_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *BdtHZZIdWP80EAPtBarrel1 = (TH1F*)file->Get("BdtHZZIdWP80EAPtBarrel1_Eff");
  TH1F *BdtHZZIdWP80EAPtBarrel2 = (TH1F*)file->Get("BdtHZZIdWP80EAPtBarrel2_Eff");
  TH1F *BdtHZZIdWP80EAPtEndcap1 = (TH1F*)file->Get("BdtHZZIdWP80EAPtEndcap1_Eff");
  TH1F *BdtHZZIdWP80EAPtEndcap2 = (TH1F*)file->Get("BdtHZZIdWP80EAPtEndcap2_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(BdtHWWIdWP80PtBarrel1);
  ptSet1.push_back(BdtHWWIdWP80PtBarrel2);
  ptSet1.push_back(BdtHWWIdWP80PtEndcap1);
  ptSet1.push_back(BdtHWWIdWP80PtEndcap2);
  ptSet2.push_back(BdtHZZIdWP80EAPtBarrel1);
  ptSet2.push_back(BdtHZZIdWP80EAPtBarrel2);
  ptSet2.push_back(BdtHZZIdWP80EAPtEndcap1);
  ptSet2.push_back(BdtHZZIdWP80EAPtEndcap2);

  drawOneToOne(ptSet1,ptSet2,"H #rightarrow WW 2011","H #rightarrow ZZ 2012","p_{T} [GeV]");

  // PU
  TH1F *BdtHWWIdWP80PU = (TH1F*)file->Get("BdtHWWIdWP80PU_Eff");
  // ---> EA corrected
  TH1F *BdtHZZIdWP80EAPU = (TH1F*)file->Get("BdtHZZIdWP80EAPU_Eff");
  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(BdtHWWIdWP80PU);
  puSet2.push_back(BdtHZZIdWP80EAPU);

  drawOneToOne(puSet1,puSet2,"H #rightarrow WW 2011","H #rightarrow ZZ 2012","# vertices");

}
