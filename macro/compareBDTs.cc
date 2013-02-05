// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include "TObjArray.h"
#include "TObjString.h"

void makeBDTDistributions(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",30,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",30,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",30,-1.0,1.0);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",30,-1.0,1.0);

  treeSig->Project("bdthww_sig","bdthzz[3]",cutSig);
  treeBkg->Project("bdthww_bkg","bdthzz[3]",cutBkg);
  treeSig->Project("bdthzz_sig","newbdthww[3]",cutSig);
  treeBkg->Project("bdthzz_bkg","newbdthww[3]",cutBkg);

  bdthww_sig->Sumw2();
  bdthzz_sig->Sumw2();

  bdthww_sig->SetLineColor(kAzure-7);
  bdthzz_sig->SetLineColor(kAzure-7);
  bdthww_sig->SetMarkerColor(kAzure-7);
  bdthzz_sig->SetMarkerColor(kAzure-7);
  bdthww_sig->SetMarkerStyle(kFullDotLarge);
  bdthzz_sig->SetMarkerStyle(kFullDotLarge);
  bdthww_sig->SetFillColor(kAzure-6);
  bdthzz_sig->SetFillColor(kAzure-6);

  bdthww_bkg->SetLineColor(kRed);
  bdthzz_bkg->SetLineColor(kRed);
  bdthww_bkg->SetMarkerColor(kRed);
  bdthzz_bkg->SetMarkerColor(kRed);
  bdthww_sig->SetMarkerStyle(kFullDotLarge);
  bdthzz_sig->SetMarkerStyle(kFullDotLarge);
  bdthww_sig->SetFillStyle(3004);
  bdthzz_sig->SetFillStyle(3004);
  bdthww_bkg->SetFillColor(kRed+1);
  bdthzz_bkg->SetFillColor(kRed+1);
  bdthww_bkg->SetFillStyle(3005);
  bdthzz_bkg->SetFillStyle(3005);

  bdthww_bkg->SetMarkerStyle(8);
  bdthzz_bkg->SetMarkerStyle(8);  

  // draw the legend
  TLegend* legend = new TLegend(0.14, 0.70, 0.37, 0.85);
  
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);
  // legend->AddEntry(bdthww_sig, "Z #rightarrow ee 2011 data");
  // legend->AddEntry(bdthww_bkg, "Z #rightarrow ee 2012 data");
  //  legend->AddEntry(bdthww_bkg, "fake rate data");
  legend->AddEntry(bdthww_bkg, "Z(#rightarrow ll)+1 fake 2011");
  legend->AddEntry(bdthww_bkg, "Z(#rightarrow ll)+1 fake 2012");
  //  legend->AddEntry(bdthww_bkg, "W(#rightarrow l#nu)+1 fake");

  // legend->AddEntry(bdthww_sig, "Z(#rightarrow ee) MC Summer12");
  // legend->AddEntry(bdthww_bkg, "Z(#rightarrow ee) data 2012");

  TCanvas c1("c1","",600,600);
  c1.Divide(1,2);
  c1.cd(1);
  gPad->SetLogy();

  bdthww_sig->GetXaxis()->SetTitle("electron HWW 2011 BDT");
  bdthzz_sig->GetXaxis()->SetTitle("electron e#gamma 2012 BDT");

  float bkgnorm = bdthww_bkg->Integral();
  float signorm = bdthww_sig->Integral();
  bdthww_sig->Scale(bkgnorm/signorm);
  bdthzz_sig->Scale(bkgnorm/signorm);

  float max = TMath::Max(bdthww_sig->GetMaximum(),bdthww_bkg->GetMaximum());
  max = max + 0.2 * max;
  bdthww_sig->SetMaximum(max);
  max = TMath::Max(bdthzz_sig->GetMaximum(),bdthzz_bkg->GetMaximum());
  max = max + 0.2 * max;
  bdthzz_sig->SetMaximum(max);

  bdthww_sig->Draw("hist");
  bdthww_sig->Draw("same pe1");
  bdthww_bkg->Draw("same hist");
  bdthww_bkg->Draw("same pe1");

  c1.cd(2);
  gPad->SetLogy();
  bdthzz_sig->Draw("hist");
  bdthzz_sig->Draw("same pe1");
  bdthzz_bkg->Draw("same hist");
  bdthzz_bkg->Draw("same pe1");
  legend->Draw();
  
  TString fullname(namefile);
  TObjArray *tokens = fullname.Tokenize(".");
  const char *basename = (((TObjString*)(*tokens)[0])->GetString()).Data();

  TString pdf = TString("figs/")+TString(basename)+TString(".pdf");
  TString png = TString("figs/")+TString(basename)+TString(".png");
  TString root = TString("figs/")+TString(basename)+TString(".root");

  c1.SaveAs(pdf);
  c1.SaveAs(png);
  c1.SaveAs(root);

}

void makeInputVarDistributions(TString var, TString title, pair<float,float> range, TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *var_sig = new TH1F("var_sig","",30,range.first,range.second);
  TH1F *var_bkg = new TH1F("var_bkg","",30,range.first,range.second);

  var_sig->Sumw2();
  var_bkg->Sumw2();
  
  treeSig->Project("var_sig",var,cutSig);
  treeBkg->Project("var_bkg",var,cutBkg);

  float bkgnorm = var_bkg->Integral();
  float signorm = var_sig->Integral();
  var_sig->Scale(bkgnorm/signorm);
  var_sig->Scale(bkgnorm/signorm);

  var_sig->Scale(signorm/bkgnorm);
  float max = TMath::Max(var_sig->GetMaximum(),var_bkg->GetMaximum());
  max = max + 0.2 * max;
  var_sig->SetMaximum(max);

  var_sig->SetLineColor(kAzure+3);
  var_sig->SetLineWidth(2);
  var_sig->SetFillColor(kAzure+5);
  //  var_sig->SetMarkerColor(kAzure-7);
  var_sig->SetFillStyle(1001);
  //  var_sig->SetMarkerStyle(kFullDotLarge);

  var_sig->GetXaxis()->SetLabelFont(42);
  var_sig->GetXaxis()->SetLabelSize(0.035);
  var_sig->GetXaxis()->SetTitleSize(0.08);
  var_sig->GetXaxis()->SetTitleFont(42);
  var_sig->GetYaxis()->SetLabelFont(42);
  var_sig->GetYaxis()->SetLabelSize(0.035);
  var_sig->GetYaxis()->SetTitleSize(0.08);
  var_sig->GetYaxis()->SetTitleOffset(0.2);
  var_sig->GetYaxis()->SetTitleFont(42);
  var_sig->GetZaxis()->SetLabelFont(42);
  var_sig->GetZaxis()->SetLabelSize(0.035);
  var_sig->GetZaxis()->SetTitleSize(0.035);
  var_sig->GetZaxis()->SetTitleFont(42);


  var_bkg->SetLineColor(kOrange+2);
  var_bkg->SetLineWidth(2);
  //  var_bkg->SetMarkerColor(kRed);
  var_bkg->SetFillColor(kOrange+1);
  var_bkg->SetFillStyle(3002);

  //  var_bkg->SetMarkerStyle(8);

  // draw the legend
  TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.031);
  legend->AddEntry(var_sig, " Z(ee) MC", "f");
  // legend->AddEntry(var_sig, "Z #rightarrow ee 2012 data", "pe");
  //  legend->AddEntry(var_bkg, "fake rate data");
  //legend->AddEntry(var_bkg, "Z(#rightarrow ee/#mu#mu) + 1 electron", "pe");
  legend->AddEntry(var_bkg, " Z+1 fake data", "f");
  // legend->AddEntry(var_sig, "W(#rightarrow e/#mu #nu) + 1 electron", "pe");

  // legend->AddEntry(var_sig, "Z(#rightarrow ee) simulation");
  // legend->AddEntry(var_bkg, "Z(#rightarrow ee) data");

  // cosmetics
  TPaveText *text = new TPaveText(0.15,0.90,0.77,0.98,"brNDC");
  text->AddText("CMS Preliminary           #sqrt{s} = 8 TeV,  L = 19.6 fb^{-1}");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(132);
  text->SetTextSize(0.04);

  TCanvas c1("c1","",600,600);
  c1.Range(-1.146789,-2319.078,5.688073,12419.95);
  c1.SetFillColor(0);
  c1.SetBorderMode(0);
  c1.SetBorderSize(2);
  c1.SetLeftMargin(0.1677852);
  c1.SetFrameBorderMode(0);
  c1.SetFrameBorderMode(0);
  c1.cd();

  if(var.Contains("HoE")||var.Contains("Iso")||var.Contains("bdthzz")) {
    var_sig->SetMinimum(1);
    c1.SetLogy();
  }

  var_sig->GetXaxis()->SetTitle(title);
  var_sig->GetYaxis()->SetTitle("events");
  var_sig->GetYaxis()->SetTitleOffset(1.8);
  var_sig->GetXaxis()->SetTitleOffset(1.0);
  var_sig->GetXaxis()->SetTitleSize(0.04);
  var_sig->GetYaxis()->SetTitleSize(0.04);
  var_sig->Draw("hist");
  //  var_sig->Draw("same pe1");
  var_bkg->Draw("same hist");
  //  var_bkg->Draw("same pe1");

  legend->Draw();
  text->Draw();
  
  TString fullname(namefile);
  fullname.ReplaceAll("/","Over");
  fullname.ReplaceAll("[","_");
  fullname.ReplaceAll("]","_");
  TObjArray *tokens = fullname.Tokenize(".");
  const char *basename = (((TObjString*)(*tokens)[0])->GetString()).Data();

  TString pdf = TString("figs/")+TString(basename)+TString(".pdf");
  TString png = TString("figs/")+TString(basename)+TString(".png");
  TString macro = TString("figs/")+TString(basename)+TString(".C");

  c1.SaveAs(pdf);
  c1.SaveAs(png);
  c1.SaveAs(macro);

}

void compareBDTs(bool applydenom) {

  // sig/bkg comparison
  // if applydenom=1 => signal = Z->ee (data), bkg = WW fake rate sample (data). This is because with denom applied Z->ee is quite clean
  // if applydenom=0 => signal = Z->ee (mc), bkg = W->ln + 1jet (data)

  gStyle->SetOptStat(0);

  TFile *fileSig, *fileBkg;
  TTree *treeSig, *treeBkg;
  fileSig = fileBkg = 0;
  treeSig = treeBkg = 0;

  if(applydenom) {
    // fileSig = TFile::Open("results_data/electrons_zeemc.root");
    fileSig = TFile::Open("results_data_2011/electrons.root");
    //    fileBkg = TFile::Open("results_data_2012/electrons.root");

    // fileSig = TFile::Open("results_data_2011/fakes-zeeOneFake.root");
    // fileBkg = TFile::Open("results_data_2012/fakes-zeeOneFake.root");

    //fileBkg = TFile::Open("results_data/fakes.root");
    //fileBkg = TFile::Open("results_data/electrons.root");
  } else {
    fileSig = TFile::Open("results_data_2012/electrons_zeemc.root");
    //fileSig = TFile::Open("results_data_2012/fakes-unbiased-wlnu.root");
    // fileBkg = TFile::Open("results_data_2012/electrons.root");
    fileBkg = TFile::Open("results_data_2012/fakes-zll1e.root");
    // fileBkg = TFile::Open("results_data_2012/fakes-zeeOneFake.root");
    // fileSig = TFile::Open("results_data/fakes-zeeOneFake.root");
    //    fileBkg = TFile::Open("results_data/fakes-unbiased-wlnu.root");
  }
  if( fileSig && fileBkg ) {
    fileSig->cd();
    treeSig = (TTree*)fileSig->Get("eleIDdir/T1");
    fileBkg->cd();
    treeBkg = (TTree*)fileBkg->Get("eleIDdir/T1");
  } else {
    cout << "File " << fileSig << " or " << fileBkg << " not existing !" << endl;
    return;
  }

  if(!treeSig || !treeBkg) {
    cout << "Tree eleIDdir/T1 not existing inside signal or background files!" << endl;
    return;
  }

  treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data_2012/fakes-zll1e_hzzisoFriend.root");
  //  treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data_2012/fakes-zeeOneFake_hzzisoFriend.root");
  //treeSig->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data_2012/fakes-unbiased-wlnu_hzzisoFriend.root");
  treeSig->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data_2012/electrons_zeemc_hzzisoFriend.root");
  //treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data_2012/electrons_hzzisoFriend.root");
  //  treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data/electrons_hzzisoFriend.root");
  //  treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","results_data/fakes_hzzisoFriend.root");

  cout << "All files OK!" << endl;

  vector<TString> cutBase;
  cutBase.push_back(TString("abs(eta)<1.0                   && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt<20"));
  cutBase.push_back(TString("abs(eta)<1.0                   && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt>20"));

  vector<TString> cutSignal;
  for(int i=0;i<(int)cutBase.size();++i) {
    if(applydenom) cutSignal.push_back(cutBase[i]+TString("&& DenomFake && abs(mass-91.1876)<7.5"));
    else cutSignal.push_back(cutBase[i]+TString("&& mcmatch && pt<30"));
    // else cutSignal.push_back(cutBase[i]+TString("&& pt<30"));
  }

  vector<TString> cutBackground;
  for(int i=0;i<(int)cutBase.size();++i) {
    //    if(applydenom) cutBackground.push_back(cutBase[i]+TString("&& DenomFakeSmurf"));
    if(applydenom) cutBackground.push_back(cutBase[i]+TString("&& DenomFake && abs(mass-91.1876)<7.5"));
    else cutBackground.push_back(cutBase[i] + TString("&& pt<30"));
  }

  vector<TString> id;
  id.push_back(TString("electronsEleMVA_inEB_LowPt.pdf"));
  id.push_back(TString("electronsEleMVA_outEB_LowPt.pdf"));
  id.push_back(TString("electronsEleMVA_EE_LowPt.pdf"));
  id.push_back(TString("electronsEleMVA_inEB_HighPt.pdf"));
  id.push_back(TString("electronsEleMVA_outEB_HighPt.pdf"));
  id.push_back(TString("electronsEleMVA_EE_HighPt.pdf"));

  vector<TString> suffix;
  suffix.push_back("_inEB_LowPt");
  suffix.push_back("_outEB_LowPt");
  suffix.push_back("_EE_LowPt");
  suffix.push_back("_inEB_HighPt");
  suffix.push_back("_outEB_HighPt");
  suffix.push_back("_EE_HighPt");

  // for(int i=0;i<(int)cutBase.size();++i) 
  //   makeBDTDistributions(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);

  vector<TString> input;
  input.push_back("bdthzz[3]");
  input.push_back("pt");
  input.push_back("EoPout");
  input.push_back("eleEoPout");
  input.push_back("IoEmIoP");
  input.push_back("EoP");
  input.push_back("HoE");
  input.push_back("eledeta");
  input.push_back("deta");
  input.push_back("dphi");
  input.push_back("see");
  input.push_back("sep");
  input.push_back("spp");
  input.push_back("phiwidth");
  input.push_back("etawidth");
  input.push_back("missHits");
  input.push_back("fbrem");
  input.push_back("nbrem");
  input.push_back("dist");
  input.push_back("dcot");
  input.push_back("d0");
  input.push_back("ip3d");
  input.push_back("ip3ds");
  input.push_back("kfhits");
  input.push_back("kfchi2");
  input.push_back("e1x5e5x5");
  input.push_back("ecaldriven");
  input.push_back("(trkIso+rho*TMath::Pi()*0.3*0.3)/pt"); // rho has been subtracted when making the ntuple
  input.push_back("(ecalIso+rho*TMath::Pi()*0.3*0.3)/pt");
  input.push_back("(hcalIso+rho*TMath::Pi()*0.3*0.3)/pt");
  input.push_back("chaPFIso[3]/pt");
  input.push_back("neuPFIso[3]/pt");
  input.push_back("phoPFIso[3]/pt");
  input.push_back("combPFIsoHWW/pt");
  input.push_back("combPFIsoHZZ/pt");
  input.push_back("vertices");

  vector<TString> title;
  title.push_back("H #rightarrow ZZ electron BDT");
  title.push_back("p_{T} [GeV]");
  title.push_back("E_{seed}/p_{out}");
  title.push_back("E_{match cluster}/p_{out}");
  title.push_back("1/E-1/P");
  title.push_back("E_{SC}/p_{in}");
  title.push_back("H/E");
  title.push_back("#Delta #eta (track-cluster)");
  title.push_back("#Delta #eta");
  title.push_back("#Delta #phi");
  title.push_back("#sigma_{i#eta i#eta}");
  title.push_back("#sigma_{i#eta i#phi}");
  title.push_back("#sigma_{i#phi i#phi}");
  title.push_back("#phi width");
  title.push_back("#eta width");
  title.push_back("miss hits");
  title.push_back("f_{brem}");
  title.push_back("n_{brem}");
  title.push_back("conv. dist");
  title.push_back("conv. ctg #theta");
  title.push_back("d_{0}");
  title.push_back("IP 3D");
  title.push_back("IP 3D sign.");
  title.push_back("KF hits");
  title.push_back("KF #chi^{2}");
  title.push_back("(E5x5-E1x5)/E5x5");
  title.push_back("ECAL seed");
  title.push_back("rel. tracker iso.");
  title.push_back("rel. ECAL iso.");
  title.push_back("rel. HCAL iso.");
  title.push_back("rel. charged PF iso.");
  title.push_back("rel. nh PF iso.");
  title.push_back("rel. #gamma PF iso.");
  title.push_back("rel. comb. PF iso. old");
  title.push_back("rel. comb. PF iso. new");
  title.push_back("n. vertices");

  vector< pair<float,float> > range;
  range.push_back(std::make_pair(-1.0,1.0)); // HZZ BDT
  range.push_back(std::make_pair(0.0,35.)); // pt
  range.push_back(std::make_pair(0.0,5.0)); // EoPout
  range.push_back(std::make_pair(0.0,5.0)); // eleEoPout
  range.push_back(std::make_pair(-0.1,0.1)); // 1/E - 1/P
  range.push_back(std::make_pair(0.0,3.0)); // EoP
  range.push_back(std::make_pair(0.0,0.15)); // H/E
  range.push_back(std::make_pair(-0.03,0.03)); // ele deta
  range.push_back(std::make_pair(-0.01,0.01)); // deta in
  range.push_back(std::make_pair(-0.2,0.2)); // dphi in
  range.push_back(std::make_pair(0.0,0.04)); // see
  range.push_back(std::make_pair(0.0,0.08)); // sep
  range.push_back(std::make_pair(0.0,0.08)); // spp
  range.push_back(std::make_pair(0.0,0.04)); // eta width
  range.push_back(std::make_pair(0.0,0.08)); // phi width
  range.push_back(std::make_pair(0.0,10)); // miss hits
  range.push_back(std::make_pair(-0.2,1.0)); // fbrem
  range.push_back(std::make_pair(0.0,5.0)); // nbrem
  range.push_back(std::make_pair(-2.0,2.0)); // dist
  range.push_back(std::make_pair(-20.0,20.)); // dcot
  range.push_back(std::make_pair(-0.04,0.04)); // d0
  range.push_back(std::make_pair(-0.05,0.05)); // ip3d
  range.push_back(std::make_pair(-10.0,10.0)); // ip3ds
  range.push_back(std::make_pair(0,20)); // kfhits
  range.push_back(std::make_pair(0,10)); // kfchi2
  range.push_back(std::make_pair(0,1)); // (e5x5-e1x5)/e5x5
  range.push_back(std::make_pair(0,1)); // ECAL seed
  range.push_back(std::make_pair(0,0.3)); // rel tk iso
  range.push_back(std::make_pair(-0.5,0.3)); // rel ECAL iso
  range.push_back(std::make_pair(0.,0.3)); // rel HCAL iso
  range.push_back(std::make_pair(0,2.0)); // cha PF iso
  range.push_back(std::make_pair(0,2.0)); // nh PF iso
  range.push_back(std::make_pair(0,2.0)); // pho PF iso
  range.push_back(std::make_pair(0,4.0)); // comb PF iso HWW
  range.push_back(std::make_pair(0,4.0)); // comb PF iso HZZ
  range.push_back(std::make_pair(0,30)); // vertices

  for(int v=0;v<(int)input.size();++v) {
    for(int i=0;i<(int)cutBase.size();++i) {
      TString filen(input[v]);
      filen.Append(suffix[i]);
      makeInputVarDistributions(input[v],title[v],range[v],treeSig,treeBkg,cutSignal[i],cutBackground[i],filen);
    }
  }

}

