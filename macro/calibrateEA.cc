#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

void calibrateEA(const char* filename, const char *cut="1") {

  gStyle->SetOptStat(0);

  TFile *file = 0;
  TTree *tree = 0;

  file = TFile::Open(filename);
  if( file ) {
    file->cd();
    tree = (TTree*)file->Get("eleIDdir/T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  if(!tree) {
    cout << "Tree eleIDdir/T1 not existing inside file " << filename << "!" << endl;
    return;
  }

  TProfile *rho_prof = new TProfile("rho_prof","rho_prof",30,1,30);
  TProfile *cha_prof = new TProfile("cha_prof","cha_prof",30,1,30);
  TProfile *neu_prof = new TProfile("neu_prof","neu_prof",30,1,30);
  TProfile *pho_prof = new TProfile("pho_prof","pho_prof",30,1,30);

  TString fullCut("(");
  fullCut += TString(cut);
  fullCut += TString(" && CutBasedIdOlyID[3] && abs(mass-91.1876)<7.5)"); // requiring WP80 ID only

  tree->Project("rho_prof","0.4*0.4*3.14*rho:vertices",fullCut.Data());
  tree->Project("cha_prof","chaPFIso:vertices",fullCut.Data());
  //  tree->Project("cha_prof","combPFIsoHWW:vertices",fullCut.Data());
  tree->Project("neu_prof","neuPFIso:vertices",fullCut.Data());
  tree->Project("pho_prof","phoPFIso:vertices",fullCut.Data());

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

  rho_prof->GetXaxis()->SetTitle("# vertices");
  rho_prof->GetYaxis()->SetTitle("<E flow> [GeV]");
  rho_prof->SetMarkerStyle(24);
  rho_prof->SetMarkerColor(kBlack);

  cha_prof->SetMarkerStyle(20);
  cha_prof->SetMarkerColor(kAzure-6);

  neu_prof->SetMarkerStyle(20);
  neu_prof->SetMarkerColor(kRed+1);

  pho_prof->SetMarkerStyle(20);
  pho_prof->SetMarkerColor(kTeal+3);
  
  rho_prof->Draw();
  cha_prof->Draw("same pe1");
  neu_prof->Draw("same pe1");
  pho_prof->Draw("same pe1");

  // draw the legend
  TLegend* legend = new TLegend(0.24, 0.74, 0.47, 0.90);
    
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);
    
  legend->AddEntry(rho_prof, "#rho (#Delta R=0.4)");
  legend->AddEntry(cha_prof, "charged particles");
  legend->AddEntry(neu_prof, "neutral hadrons");
  legend->AddEntry(pho_prof, "photons");

  legend->Draw();

  c1->SaveAs("isopf_EAcalib.eps");
  
}
