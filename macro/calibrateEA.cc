#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

// GLOBAL VARIABLES
double rho[2]; // val[0] +/- err[1]
double neu[2];
double pho[2];

vector<double> Aeff_neu, Aeff_pho;
vector<double> Aeff_neu_err, Aeff_pho_err;

Double_t polfun(Double_t *x, Double_t *par) {
  // par[0] = linear coefficient
  // par[1] = offset
  return par[0]*x[0] + par[1];
}

void calibrateEA(const char* filename, const char* outfilename, const char *cut="1", double max=24, int rebin=1) {

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

  TProfile *rho_prof = new TProfile("rho_prof","",30,1,30);
  TProfile *cha_prof = new TProfile("cha_prof","",30,1,30);
  TProfile *neu_prof = new TProfile("neu_prof","",30,1,30);
  TProfile *pho_prof = new TProfile("pho_prof","",30,1,30);

  TString fullCut("(");
  fullCut += TString(cut);
  fullCut += TString(" && CutBasedIdOlyID[3] && abs(mass-91.1876)<7.5)"); // requiring WP80 ID only

  tree->Project("rho_prof","rho:vertices",fullCut.Data());
  tree->Project("cha_prof","chaPFIso[3]:vertices",fullCut.Data());
  // this if you want separate NH/PHO
  tree->Project("neu_prof","neuPFIso[3]:vertices",fullCut.Data());
  tree->Project("pho_prof","phoPFIso[3]:vertices",fullCut.Data());
  // this is combined NH+PHO
  // tree->Project("neu_prof","neuPFIso[3]+phoPFIso[3]:vertices",fullCut.Data());

  rho_prof->Rebin(rebin);
  cha_prof->Rebin(rebin);
  neu_prof->Rebin(rebin);
  pho_prof->Rebin(rebin);

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
  
  // draw the legend
  TLegend* legend = new TLegend(0.20, 0.70, 0.43, 0.86);
    
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);
    
  legend->AddEntry(rho_prof, "#rho");
  legend->AddEntry(cha_prof, "charged particles");
  legend->AddEntry(neu_prof, "neutral hadrons");
  legend->AddEntry(pho_prof, "photons");

  // =========================
  // And now apply the fit
  // =========================

  // fit NH
  TF1 *f2nh = new TF1("mypol2nh",polfun,2.,max,2);
  f2nh->SetParNames("a","b");
  f2nh->SetLineWidth(2);
  f2nh->SetLineColor(kRed+1);
  neu_prof->Fit("mypol2nh","RB");
  neu[0] = f2nh->GetParameter(0);
  neu[1] = f2nh->GetParError(0);
  //  float neu_ChiSqr = f2->GetChisquare();  // obtain chi^2
  //  float neu_NDF = f2->GetNDF();           // obtain ndf
  std::cout << "===> NEU HAD FIT: A = " << neu[0] << " +/- " << neu[1] << std::endl;

  // fit GAMMAs
  TF1 *f2pho = new TF1("mypol2pho",polfun,2.,max,2);
  f2pho->SetParNames("a","b");
  f2pho->SetLineWidth(2);
  f2pho->SetLineColor(kTeal+3);
  pho_prof->Fit("mypol2pho","RB");
  pho[0] = f2pho->GetParameter(0);
  pho[1] = f2pho->GetParError(0);
  std::cout << "===> GAMMAs FIT: A = " << pho[0] << " +/- " << pho[1] << std::endl;

  // fit rho
  TF1 *f1 = new TF1("mypol",polfun,2,max,2);
  f1->SetParNames("a","b");
  f1->SetLineWidth(2);
  f1->SetLineColor(kBlack);
  rho_prof->Fit("mypol","RB");
  rho[0] = f1->GetParameter(0);
  rho[1] = f1->GetParError(0);
  std::cout << "===> RHO FIT: A = " << rho[0] << " +/- " << rho[1] << std::endl;
  
  rho_prof->Draw("same");
  cha_prof->Draw("same pe1");
  neu_prof->Draw("same pe1");
  pho_prof->Draw("same pe1");
  legend->Draw();
  c1->SaveAs(outfilename);

  TString png(outfilename);
  png.ReplaceAll("pdf","png");
  c1->SaveAs(png);

}

void calcAllAeff() {

  vector<TString> cut;
  cut.push_back(TString("abs(eta)<1.0"));                    // 0
  cut.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479"));  // 1 
  cut.push_back(TString("abs(eta)>1.479 && abs(eta)<2.0"));  // 2
  cut.push_back(TString("abs(eta)>2.0 && abs(eta)<2.2"));    // 3
  cut.push_back(TString("abs(eta)>2.2 && abs(eta)<2.3"));    // 4
  cut.push_back(TString("abs(eta)>2.3 && abs(eta)<2.4"));    // 5
  cut.push_back(TString("abs(eta)>2.4"));                    // 6
  
  vector<TString> id;
  id.push_back(TString("isopf_EAcalib_EB_10.pdf"));
  id.push_back(TString("isopf_EAcalib_EB_15.pdf"));
  id.push_back(TString("isopf_EAcalib_EE_20.pdf"));
  id.push_back(TString("isopf_EAcalib_EE_22.pdf"));
  id.push_back(TString("isopf_EAcalib_EE_23.pdf"));
  id.push_back(TString("isopf_EAcalib_EE_24.pdf"));
  id.push_back(TString("isopf_EAcalib_EE_30.pdf"));

  vector<float> max;
  for(int i=0; i<(int)cut.size()-3; i++) max.push_back(24);
  max[4] = 20;
  max[5] = 15;
  max[6] = 15;
  
  vector<float> rebin;
  for(int i=0; i<(int)cut.size(); i++) rebin.push_back(1);

  Aeff_neu.reserve(id.size());
  Aeff_neu_err.reserve(id.size());
  Aeff_pho.reserve(id.size());
  Aeff_pho_err.reserve(id.size());

  for(int i=0; i<(int)cut.size(); ++i) {
    calibrateEA("results_data/electrons.root",id[i].Data(),cut[i].Data(),max[i],rebin[i]);
    Aeff_neu[i] = neu[0]/rho[0];
    Aeff_neu_err[i] = Aeff_neu[i] * sqrt(pow(neu[1]/neu[0],2)+pow(rho[1]/rho[0],2));
    Aeff_pho[i] = pho[0]/rho[0];
    Aeff_pho_err[i] = Aeff_pho[i] * sqrt(pow(pho[1]/pho[0],2)+pow(rho[1]/rho[0],2));
  }

  cout << " ==== RESULTS OF EFFECTIVE AREAS==== " << endl;
  for(int i=0; i<(int)cut.size(); ++i) {
    cout << "CUT = " << cut[i].Data();
    cout << "\tAeff(NH) = " << Aeff_neu[i] << " +/- " << Aeff_neu_err[i] << "\t\t";
    cout << "\tAeff(PHO) = " << Aeff_pho[i] << " +/- " << Aeff_pho_err[i] << endl;    
  }
}


void rhoEASubtracted(const char* filename, const char* outfilename, float aeffneu, float aeffpho, const char *cut="1", double max=24, int rebin=1) {
  
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

  TProfile *rho_prof = new TProfile("rho_prof","",30,1,30);
  TProfile *iso_prof = new TProfile("iso_prof","",30,1,30);
  TProfile *isocorr_prof = new TProfile("isocorr_prof","",30,1,30);

  TString fullCut("(");
  fullCut += TString(cut);
  fullCut += TString(" && CutBasedIdOlyID[3] && abs(mass-91.1876)<7.5)"); // requiring WP80 ID only

  tree->Project("rho_prof","rho:vertices",fullCut.Data());
  tree->Project("iso_prof","chaPFIso[3]+neuPFIso[3]+phoPFIso[3]:vertices",fullCut.Data());

  // form the formula for the corrected isolation from calibrated EA
  char formula[1000];
  sprintf(formula,"chaPFIso[3] + (neuPFIso[3] - %f * rho) + (phoPFIso[3] - %f * rho):vertices",aeffneu,aeffpho);
  tree->Project("isocorr_prof",formula,fullCut.Data());

  rho_prof->Rebin(rebin);
  iso_prof->Rebin(rebin);
  isocorr_prof->Rebin(rebin);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

  rho_prof->GetXaxis()->SetTitle("# vertices");
  rho_prof->GetYaxis()->SetTitle("<E flow> [GeV]");
  rho_prof->SetMarkerStyle(24);
  rho_prof->SetMarkerColor(kBlack);

  iso_prof->SetMarkerStyle(20);
  iso_prof->SetMarkerColor(kMagenta+3);

  isocorr_prof->SetMarkerStyle(20);
  isocorr_prof->SetMarkerColor(kTeal+3);

  // draw the legend
  TLegend* legend = new TLegend(0.20, 0.70, 0.43, 0.86);
    
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);
    
  legend->AddEntry(rho_prof, "#rho");
  legend->AddEntry(iso_prof, "PF isolation");
  legend->AddEntry(isocorr_prof, "PF isolation, EA corr");

  // ==============================
  // And now apply the fit (verify)
  // ==============================

  // fit iso corr
  TF1 *f1 = new TF1("mypol",polfun,2.,max,2);
  f1->SetParNames("a","b");
  f1->SetLineWidth(2);
  f1->SetLineColor(kTeal+3);
  isocorr_prof->Fit("mypol","RB");
  float isocorr_coeff = f1->GetParameter(0);
  float isocorr_coeff_err = f1->GetParError(0);
  std::cout << "===> ISO FIT: A = " << isocorr_coeff << " +/- " << isocorr_coeff_err << std::endl;
  
  isocorr_prof->SetMaximum(20);
  isocorr_prof->SetMinimum(-1);
  rho_prof->Draw("same");
  iso_prof->Draw("same pe1");
  isocorr_prof->Draw("same pe1");
  legend->Draw();
  c1->SaveAs(outfilename);

  TString png(outfilename);
  png.ReplaceAll("pdf","png");
  c1->SaveAs(png);

}

void doAll() {

  // first, calibrate the effective area (EA) in 4 bins
  calcAllAeff();

  vector<TString> cut;
  cut.push_back(TString("abs(eta)<1.0"));
  cut.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479"));
  cut.push_back(TString("abs(eta)>1.479 && abs(eta)<2.0"));
  cut.push_back(TString("abs(eta)>2.0 && abs(eta)<2.2"));
  cut.push_back(TString("abs(eta)>2.2 && abs(eta)<2.3"));
  cut.push_back(TString("abs(eta)>2.3 && abs(eta)<2.4"));
  cut.push_back(TString("abs(eta)>2.4"));
  
  vector<TString> id;
  id.push_back(TString("isocorrpf_EAcalib_EB_10.pdf"));
  id.push_back(TString("isocorrpf_EAcalib_EB_15.pdf"));
  id.push_back(TString("isocorrpf_EAcalib_EE_20.pdf"));
  id.push_back(TString("isocorrpf_EAcalib_EE_22.pdf"));
  id.push_back(TString("isocorrpf_EAcalib_EE_23.pdf"));
  id.push_back(TString("isocorrpf_EAcalib_EE_24.pdf"));
  id.push_back(TString("isocorrpf_EAcalib_EE_30.pdf"));

  vector<float> max;
  for(int i=0; i<(int)cut.size()-3; i++) max.push_back(24);
  max[4] = 20;
  max[5] = 20;
  max[6] = 20;

  vector<float> rebin;
  for(int i=0; i<(int)cut.size(); i++) rebin.push_back(1);
  
  // plot the raw / corrected isolation
  for(int i=0; i<(int)cut.size(); ++i) 
    rhoEASubtracted("results_data/electrons.root",id[i].Data(),Aeff_neu[i],Aeff_pho[i],cut[i].Data(),max[i],rebin[i]);

}
