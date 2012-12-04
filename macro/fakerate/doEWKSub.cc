#include <iostream>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>
#include <TLegend.h>
#include <TPaveText.h>

using namespace std;

// constants: hardcoded
struct constants {
  Float_t lumiprescaled, nWlnuProc, nDYllProc, xsecWlnu, xsecDYll;
  Float_t xsecDYllMeas, xsecWlnuMeas;
};

TH1F* drawOne(TH1F* dataD, TH1F* dataN, TH1F* WD, TH1F* WN, TH1F* ZD, TH1F* ZN, const char *basename, const char* t) {

  constants c;
  c.lumiprescaled=12.9; // pb
  c.nWlnuProc=45.161288e+6;
  c.nDYllProc=16.844034e+6;
  c.xsecWlnu=37509.0;
  c.xsecDYll=3532.8149;

  gStyle->SetOptStat(0);

  float wEffL=c.nWlnuProc / c.xsecWlnu;
  float zEffL=c.nDYllProc / c.xsecDYll;

  cout << "Lumi presc = " << c.lumiprescaled << " anw W / Z lumi = " << wEffL << "   /   " << zEffL << endl;

  WD->Scale(c.lumiprescaled / wEffL);  WN->Scale(c.lumiprescaled / wEffL);
  ZD->Scale(c.lumiprescaled / zEffL);  ZN->Scale(c.lumiprescaled / zEffL);

  WN->Sumw2(); ZN->Sumw2();
  TH1F *ewkN = (TH1F*)WN->Clone();
  TH1F *ewkD = (TH1F*)WD->Clone();
  ewkN->Add(ZN);
  ewkD->Add(ZD);

  // first draw the Numerator for data and EWK and show the relative contribution of EWK
  TLegend* legend = new TLegend(0.15, 0.70, 0.30, 0.80);
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);

  dataN->SetMarkerColor(kBlack);
  dataN->SetLineColor(kBlack);
  dataN->SetMarkerStyle(8);
  ewkN->SetMarkerColor(kRed);
  ewkN->SetLineColor(kRed);
  ewkN->SetMarkerStyle(8);

  legend->AddEntry(dataN,"data","pl");
  legend->AddEntry(ewkN,"EWK contribution","pl");

  dataN->Sumw2(); dataD->Sumw2();
  dataN->GetXaxis()->SetRangeUser(10,35);
  dataN->GetXaxis()->SetTitle("p_{T} [GeV]");
  dataN->SetMinimum(0.0);

  TPaveText* text  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
  text->AddText("#sqrt{s} = 8 TeV, L = 12.1 fb^{-1}");
  text->AddText(t);
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(32);
  text->SetTextSize(0.05);

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(1,2);

  c1->cd(1);
  dataN->SetTitle("numerator counts");
  dataN->Draw("pe1");
  ewkN->Draw("same");
  legend->Draw();

  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();  
  TH1F* ewkRel = (TH1F*)ewkN->Clone("ratio");
  ewkRel->Sumw2();
  ewkRel->GetXaxis()->SetRangeUser(10,35);
  ewkRel->GetXaxis()->SetTitle("p_{T} [GeV]");
  ewkRel->SetTitle("relative EKW contribution");
  ewkRel->Divide(dataN);
  ewkRel->SetMarkerStyle(8);
  ewkRel->SetMarkerColor(kBlack);
  ewkRel->SetLineColor(kBlack);
  ewkRel->Draw("pe");
  
  stringstream fileContr;
  fileContr << basename << "_contr.png";
  c1->SaveAs(fileContr.str().c_str());

  
  // do the subtraction
  TH1F* dataNSub = (TH1F*)dataN->Clone("dataNSub");
  TH1F* dataDSub = (TH1F*)dataD->Clone("dataDSub");
  dataNSub->Sumw2(); dataDSub->Sumw2();
  dataNSub->Add(ewkN,-1.);
  dataDSub->Add(ewkD,-1.);

  dataNSub->Divide(dataDSub);
  dataN->Divide(dataD);

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  dataNSub->SetLineColor(kAzure);
  dataNSub->SetMarkerColor(kAzure);
  dataNSub->SetMarkerStyle(8);
  dataN->SetTitle("");
  dataN->SetMaximum(0.4);
  dataN->GetYaxis()->SetTitle("fake rate");

  TLegend* legend2 = new TLegend(0.15, 0.60, 0.30, 0.70);
  legend2->SetBorderSize(   0);
  legend2->SetFillColor (   0);
  legend2->SetTextAlign (  12);
  legend2->SetTextFont  (  42);
  legend2->SetTextSize  (0.05);
  legend2->AddEntry(dataN,"fake rate (plain)","pl");
  legend2->AddEntry(dataNSub,"fake rate (EWK-corr.)","pl");
  dataN->Draw("pe1");
  dataNSub->Draw("pe1 same");
  text->Draw();
  legend2->Draw();

  stringstream fileCorr;
  fileCorr << basename << ".png";
  c2->SaveAs(fileCorr.str().c_str());

  return dataNSub;

}

void doEWKSub() {

  TFile *fileData = TFile::Open("fakerates_trigger.root");
  TFile *fileW = TFile::Open("subewkw_trigger.root");
  TFile *fileZ = TFile::Open("subewkz_trigger.root");

  // --- DATA ---
  // pt DEN
  TH1F *RecoPtBarrel1 = (TH1F*)fileData->Get("RecoPtBarrel1");
  TH1F *RecoPtBarrel2 = (TH1F*)fileData->Get("RecoPtBarrel2");
  TH1F *RecoPtEndcap1 = (TH1F*)fileData->Get("RecoPtEndcap1");
  TH1F *RecoPtEndcap2 = (TH1F*)fileData->Get("RecoPtEndcap2");
  // pt NUM
  TH1F *TrgElenewWPHWWPtBarrel1 = (TH1F*)fileData->Get("TrgElenewWPHWWPtBarrel1");
  TH1F *TrgElenewWPHWWPtBarrel2 = (TH1F*)fileData->Get("TrgElenewWPHWWPtBarrel2");
  TH1F *TrgElenewWPHWWPtEndcap1 = (TH1F*)fileData->Get("TrgElenewWPHWWPtEndcap1");
  TH1F *TrgElenewWPHWWPtEndcap2 = (TH1F*)fileData->Get("TrgElenewWPHWWPtEndcap2");
  // -----------------------------------

  // --- WLNU ---
  // pt DEN
  TH1F *WRecoPtBarrel1 = (TH1F*)fileW->Get("RecoPtBarrel1");
  TH1F *WRecoPtBarrel2 = (TH1F*)fileW->Get("RecoPtBarrel2");
  TH1F *WRecoPtEndcap1 = (TH1F*)fileW->Get("RecoPtEndcap1");
  TH1F *WRecoPtEndcap2 = (TH1F*)fileW->Get("RecoPtEndcap2");
  // pt NUM
  TH1F *WTrgElenewWPHWWPtBarrel1 = (TH1F*)fileW->Get("TrgElenewWPHWWPtBarrel1");
  TH1F *WTrgElenewWPHWWPtBarrel2 = (TH1F*)fileW->Get("TrgElenewWPHWWPtBarrel2");
  TH1F *WTrgElenewWPHWWPtEndcap1 = (TH1F*)fileW->Get("TrgElenewWPHWWPtEndcap1");
  TH1F *WTrgElenewWPHWWPtEndcap2 = (TH1F*)fileW->Get("TrgElenewWPHWWPtEndcap2");
  // -----------------------------------

  // --- ZLL ---
  // pt DEN
  TH1F *ZRecoPtBarrel1 = (TH1F*)fileZ->Get("RecoPtBarrel1");
  TH1F *ZRecoPtBarrel2 = (TH1F*)fileZ->Get("RecoPtBarrel2");
  TH1F *ZRecoPtEndcap1 = (TH1F*)fileZ->Get("RecoPtEndcap1");
  TH1F *ZRecoPtEndcap2 = (TH1F*)fileZ->Get("RecoPtEndcap2");
  // pt NUM
  TH1F *ZTrgElenewWPHWWPtBarrel1 = (TH1F*)fileZ->Get("TrgElenewWPHWWPtBarrel1");
  TH1F *ZTrgElenewWPHWWPtBarrel2 = (TH1F*)fileZ->Get("TrgElenewWPHWWPtBarrel2");
  TH1F *ZTrgElenewWPHWWPtEndcap1 = (TH1F*)fileZ->Get("TrgElenewWPHWWPtEndcap1");
  TH1F *ZTrgElenewWPHWWPtEndcap2 = (TH1F*)fileZ->Get("TrgElenewWPHWWPtEndcap2");
  // -----------------------------------

  TFile *filecorr = TFile::Open("fakerates_trigger_ewksub_hcp.root","recreate");
  TH1F *frBarrel1 = drawOne(RecoPtBarrel1,TrgElenewWPHWWPtBarrel1,WRecoPtBarrel1,WTrgElenewWPHWWPtBarrel1,ZRecoPtBarrel1,ZTrgElenewWPHWWPtBarrel1,"ewksub_barrel1","|#eta|<1");
  TH1F *frBarrel2 = drawOne(RecoPtBarrel2,TrgElenewWPHWWPtBarrel2,WRecoPtBarrel2,WTrgElenewWPHWWPtBarrel2,ZRecoPtBarrel2,ZTrgElenewWPHWWPtBarrel2,"ewksub_barrel2","1<|#eta|<1.48");
  TH1F *frEndcap1 = drawOne(RecoPtEndcap1,TrgElenewWPHWWPtEndcap1,WRecoPtEndcap1,WTrgElenewWPHWWPtEndcap1,ZRecoPtEndcap1,ZTrgElenewWPHWWPtEndcap1,"ewksub_endcap1","1.48<|#eta|<2");
  TH1F *frEndcap2 = drawOne(RecoPtEndcap2,TrgElenewWPHWWPtEndcap2,WRecoPtEndcap2,WTrgElenewWPHWWPtEndcap2,ZRecoPtEndcap2,ZTrgElenewWPHWWPtEndcap2,"ewksub_endcap2","|#eta|>2");

  // use same names as the originals to make application easier
  frBarrel1->SetName("TrgElenewWPHWWPtBarrel1");
  frBarrel2->SetName("TrgElenewWPHWWPtBarrel2");
  frEndcap1->SetName("TrgElenewWPHWWPtEndcap1");
  frEndcap2->SetName("TrgElenewWPHWWPtEndcap2");

  filecorr->cd();
  frBarrel1->Write();
  frBarrel2->Write();
  frEndcap1->Write();
  frEndcap2->Write();
  filecorr->Close();

}
