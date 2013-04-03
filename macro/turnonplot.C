#define turnonplot_cxx
#include "turnonplot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void turnonplot::Loop()
{
  if (fChain == 0) return;

  gStyle->SetOptStat(0);

  // histos
  TH1F *Hden30_EB = new TH1F("Hden30_EB","Hden30_EB",40,20.,60.);
  TH1F *Hden50_EB = new TH1F("Hden50_EB","Hden50_EB",40,40.,80.);
  TH1F *Hden75_EB = new TH1F("Hden75_EB","Hden75_EB",40,70.,110.);
  TH1F *Hden90_EB = new TH1F("Hden90_EB","Hden90_EB",40,85.,125.);
  //
  TH1F *Hnum30_EB = new TH1F("Hnum30_EB","Hnum30_EB",40,20.,60.);
  TH1F *Hnum50_EB = new TH1F("Hnum50_EB","Hnum50_EB",40,40.,80.);
  TH1F *Hnum75_EB = new TH1F("Hnum75_EB","Hnum75_EB",40,70.,110.);
  TH1F *Hnum90_EB = new TH1F("Hnum90_EB","Hnum90_EB",40,85.,125.);
  //
  TH1F *Hden30_EE = new TH1F("Hden30_EE","Hden30_EE",40,20.,60.);
  TH1F *Hden50_EE = new TH1F("Hden50_EE","Hden50_EE",40,40.,80.);
  TH1F *Hden75_EE = new TH1F("Hden75_EE","Hden75_EE",40,70.,110.);
  TH1F *Hden90_EE = new TH1F("Hden90_EE","Hden90_EE",40,85.,125.);
  //
  TH1F *Hnum30_EE = new TH1F("Hnum30_EE","Hnum30_EE",40,20.,60.);
  TH1F *Hnum50_EE = new TH1F("Hnum50_EE","Hnum50_EE",40,40.,80.);
  TH1F *Hnum75_EE = new TH1F("Hnum75_EE","Hnum75_EE",40,70.,110.);
  TH1F *Hnum90_EE = new TH1F("Hnum90_EE","Hnum90_EE",40,85.,125.);

  // for errors
  Hden30_EB -> Sumw2();
  Hden50_EB -> Sumw2();
  Hden75_EB -> Sumw2();
  Hden90_EB -> Sumw2();
  // 
  Hnum30_EB -> Sumw2();
  Hnum50_EB -> Sumw2();
  Hnum75_EB -> Sumw2();
  Hnum90_EB -> Sumw2();
  //
  Hden30_EE -> Sumw2();
  Hden50_EE -> Sumw2();
  Hden75_EE -> Sumw2();
  Hden90_EE -> Sumw2();
  // 
  Hnum30_EE -> Sumw2();
  Hnum50_EE -> Sumw2();
  Hnum75_EE -> Sumw2();
  Hnum90_EE -> Sumw2();

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // ask a reco photon to pass the offline selection and match the HLT candidate
    if (!hltmatchGamma) continue;
    if (!idGamma)       continue;

    // ask for Z
    // if (!isAZ)       continue;
    // if (!isAZ_tight) continue;

    // preparing the numerator
    if (fabs(etaGamma)<1.5) {
      Hden30_EB -> Fill(ptGamma);
      Hden50_EB -> Fill(ptGamma);
      Hden75_EB -> Fill(ptGamma);
      Hden90_EB -> Fill(ptGamma);
      if (passHLT30) Hnum30_EB -> Fill(ptGamma);
      if (passHLT50) Hnum50_EB -> Fill(ptGamma);
      if (passHLT75) Hnum75_EB -> Fill(ptGamma);
      if (passHLT90) Hnum90_EB -> Fill(ptGamma);
    } else {
      Hden30_EE -> Fill(ptGamma);
      Hden50_EE -> Fill(ptGamma);
      Hden75_EE -> Fill(ptGamma);
      Hden90_EE -> Fill(ptGamma);
      if (passHLT30) Hnum30_EE -> Fill(ptGamma);
      if (passHLT50) Hnum50_EE -> Fill(ptGamma);
      if (passHLT75) Hnum75_EE -> Fill(ptGamma);
      if (passHLT90) Hnum90_EE -> Fill(ptGamma);
    }
  }


  // --------------------------------------------------
  // efficiency plots, HLT pT 30
  TH1F *Heff30_EB  = (TH1F*)Hnum30_EB ->Clone("Heff30_EB");
  TH1F *Heff30_EE  = (TH1F*)Hnum30_EE ->Clone("Heff30_EE");
  Heff30_EB  -> Divide(Hden30_EB);
  Heff30_EE  -> Divide(Hden30_EE);
  Heff30_EB  -> SetLineWidth(2); Heff30_EB  -> SetLineColor(4); 
  Heff30_EE  -> SetLineWidth(2); Heff30_EE  -> SetLineColor(4); 
  cout << endl;
  cout << "Heff30_EB "  << Heff30_EB->GetEntries()  << endl;
  cout << "Heff30_EE "  << Heff30_EE->GetEntries()  << endl;

  TH2F *myH2_30 = new TH2F("myH2_30","myH2_30",100,20.,50.,100,0.,1.2);
  myH2_30 -> GetXaxis()->SetTitle("reco photon pT [GeV]");
  myH2_30 -> SetTitle("HLT_Photon30_CaloIdVL_IsoL");

  TCanvas c30("c30","c30",1);
  c30.Divide(2,1);
  c30.cd(1); myH2_30->Draw(); Heff30_EB->Draw("pEsame"); 
  c30.cd(2); myH2_30->Draw(); Heff30_EE->Draw("pEsame"); 
  c30.SaveAs("pippo30.root");


  // --------------------------------------------------
  // efficiency plots, HLT pT 50
  TH1F *Heff50_EB  = (TH1F*)Hnum50_EB ->Clone("Heff50_EB");
  TH1F *Heff50_EE  = (TH1F*)Hnum50_EE ->Clone("Heff50_EE");
  Heff50_EB  -> Divide(Hden50_EB);
  Heff50_EE  -> Divide(Hden50_EE);
  Heff50_EB  -> SetLineWidth(2); Heff50_EB  -> SetLineColor(4); 
  Heff50_EE  -> SetLineWidth(2); Heff50_EE  -> SetLineColor(4); 

  cout << endl;
  cout << "Heff50_EB "  << Heff50_EB->GetEntries()  << endl;
  cout << "Heff50_EE "  << Heff50_EE->GetEntries()  << endl;

  TH2F *myH2_50 = new TH2F("myH2_50","myH2_50",100,40.,70.,100,0.,1.2);
  myH2_50 -> GetXaxis()->SetTitle("reco photon pT [GeV]");
  myH2_50 -> SetTitle("HLT_Photon50_CaloIdVL_IsoL");

  TCanvas c50("c50","c50",1);
  c50.Divide(2,1);
  c50.cd(1); myH2_50->Draw(); Heff50_EB->Draw("samepE"); 
  c50.cd(2); myH2_50->Draw(); Heff50_EE->Draw("samepE"); 
  c50.SaveAs("pippo50.root");


  // --------------------------------------------------
  // efficiency plots, HLT pT 75
  TH1F *Heff75_EB  = (TH1F*)Hnum75_EB ->Clone("Heff75_EB");
  TH1F *Heff75_EE  = (TH1F*)Hnum75_EE ->Clone("Heff75_EE");
  Heff75_EB  -> Divide(Hden75_EB);
  Heff75_EE  -> Divide(Hden75_EE);
  Heff75_EB  -> SetLineWidth(2); Heff75_EB  -> SetLineColor(4); 
  Heff75_EE  -> SetLineWidth(2); Heff75_EE  -> SetLineColor(4); 

  cout << endl;
  cout << "Heff75_EB "  << Heff75_EB->GetEntries()  << endl;
  cout << "Heff75_EE "  << Heff75_EE->GetEntries()  << endl;

  TH2F *myH2_75 = new TH2F("myH2_75","myH2_75",100,70.,100.,100,0.,1.2);
  myH2_75 -> GetXaxis()->SetTitle("reco photon pT [GeV]");
  myH2_75 -> SetTitle("HLT_Photon75_CaloIdVL_IsoL");

  TCanvas c75("c75","c75",1);
  c75.Divide(2,1);
  c75.cd(1); myH2_75->Draw(); Heff75_EB->Draw("samepE"); 
  c75.cd(2); myH2_75->Draw(); Heff75_EE->Draw("samepE"); 
  c75.SaveAs("pippo75.root");


  // --------------------------------------------------
  // efficiency plots, HLT pT 90
  TH1F *Heff90_EB  = (TH1F*)Hnum90_EB ->Clone("Heff90_EB");
  TH1F *Heff90_EE  = (TH1F*)Hnum90_EE ->Clone("Heff90_EE");
  Heff90_EB  -> Divide(Hden90_EB);
  Heff90_EE  -> Divide(Hden90_EE);
  Heff90_EB  -> SetLineWidth(2); Heff90_EB  -> SetLineColor(4); 
  Heff90_EE  -> SetLineWidth(2); Heff90_EE  -> SetLineColor(4); 

  cout << endl;
  cout << "Heff90_EB "  << Heff90_EB->GetEntries()  << endl;
  cout << "Heff90_EE "  << Heff90_EE->GetEntries()  << endl;

  TH2F *myH2_90 = new TH2F("myH2_90","myH2_90",100,80.,115.,100,0.,1.2);
  myH2_90 -> GetXaxis()->SetTitle("reco photon pT [GeV]");
  myH2_90 -> SetTitle("HLT_Photon90_CaloIdVL_IsoL");

  TCanvas c90("c90","c90",1);
  c90.Divide(2,1);
  c90.cd(1); myH2_90->Draw(); Heff90_EB->Draw("samepE"); 
  c90.cd(2); myH2_90->Draw(); Heff90_EE->Draw("samepE"); 
  c90.SaveAs("pippo90.root");
}
