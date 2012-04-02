// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>

#include "include/HZZEleIDSelector.hh"

using namespace std;


enum idType {
  kIsoHWW2011 = 0, // HWW cuts 2011
  kBDTHWW2011_withIP, // HWW cuts 2011
  kIso,
  kIsoEACorr,
  kBDTHZZ_withIP
};

float Aeff_neu_dr04[7] = { 0.045, 0.062, 0.061, 0.041, 0.050, 0.051, 0.110 };
float Aeff_pho_dr04[7] = { 0.140, 0.130, 0.080, 0.130, 0.140, 0.160, 0.180 };

bool passHWWID(Float_t eta, Float_t pt, Float_t bdthww, Float_t bdthzz, Float_t rho, Float_t combIso, Float_t combPFIsoHWW, idType type) {
  if(type == kIsoHWW2011) {
    if(fabs(eta)<1.479) return (combPFIsoHWW/pt < 0.13);
    else return (combPFIsoHWW/pt < 0.09);
  }

  if(type == kIsoEACorr) {
    if(pt>20) {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.23); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.20);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.12);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.11);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.049);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.070);
      if(fabs(eta) >=  2.4) return (combIso/pt < 0.010);
    } else {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.20); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.21);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.13);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.10);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.04);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.03);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.03);
    }
  }

  if(type == kIso) {
    if(fabs(eta) < 1.479) return (combIso/pt < 0.29);
    else return (combIso/pt < 0.21);
  }

  if(type == kBDTHWW2011_withIP) {
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.884);
  }

  if(type == kBDTHZZ_withIP) {
    // WP with same fake rate as HWW with IP
    if(pt < 20 && fabs(eta) < 1.0) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.091);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthzz > 0.064);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.071);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.067);
  }

  return false;
}

void makeFriendHZZIsolation(const char* file) {

  TFile *pF = TFile::Open(file);
  TTree *pT = (TTree*)pF->Get("eleIDdir/T1");
  
  TString nF(file);
  nF.ReplaceAll(".root","_hzzisoFriend.root");
  TFile *fF = TFile::Open(nF,"recreate");

  Float_t chaPFIso[8], neuPFIso[8], phoPFIso[8], rho, eta;
  pT->SetBranchAddress("chaPFIso", chaPFIso);
  pT->SetBranchAddress("neuPFIso", neuPFIso);
  pT->SetBranchAddress("phoPFIso", phoPFIso);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("eta", &eta);

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  Float_t combIso, combIsoNoEA;
  fT->Branch("combPFIsoHZZ",&combIso,"combPFIsoHZZ/F");
  fT->Branch("combPFIsoHZZNoEA",&combIsoNoEA,"combPFIsoHZZNoEA/F");
  
  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     combIsoNoEA = chaPFIso[3]+neuPFIso[3]+phoPFIso[3]; // [0] = cone 0.1, then steps of 0.1 up to 0.7 [7] is directional isolation with DR=0.4
     combIso = 0.0;
     if(fabs(eta) <  1.0) combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[0]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[0]*rho,Float_t(0.));
     else if(fabs(eta) < 1.479) combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[1]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[1]*rho,Float_t(0.));
     else if(fabs(eta) < 2.0) combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[2]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[2]*rho,Float_t(0.));
     else if(fabs(eta) < 2.2) combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[3]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[3]*rho,Float_t(0.));
     else if(fabs(eta) < 2.3) combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[4]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[4]*rho,Float_t(0.));
     else if(fabs(eta) < 2.4) combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[5]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[5]*rho,Float_t(0.));
     else combIso = chaPFIso[3] + TMath::Max(neuPFIso[3]-Aeff_neu_dr04[6]*rho,Float_t(0.)) + TMath::Max(phoPFIso[3]-Aeff_pho_dr04[6]*rho,Float_t(0.));
     fT->Fill();
  }

  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF.Data() << endl;

}

void makeFriendHZZIdBits(const char* file) {

  TFile *pF = TFile::Open(file);
  TTree *pT = (TTree*)pF->Get("eleIDdir/T1");

  TString isofriend = TString(file);
  isofriend.ReplaceAll(".root","_hzzisoFriend.root");
  pT->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", isofriend );
 
  TString nF(file);
  nF.ReplaceAll(".root","_hzzidbitsFriend.root");
  TFile *fF = TFile::Open(nF,"recreate");

  Float_t eta, abseta, pt, rho, vertices;
  Float_t bdthww[2], newbdthww[4], combPFIsoHWW;
  Float_t combPFIsoHZZ, bdthzz[4];
  Float_t chaPFIso[8], neuPFIso[8], phoPFIso[8];
  Float_t mass; // not dummy only for TP trees
  pT->SetBranchAddress("bdthww", bdthww);
  pT->SetBranchAddress("newbdthww", newbdthww);
  pT->SetBranchAddress("bdthzz",bdthzz);
  pT->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
  pT->SetBranchAddress("combPFIsoHZZ",&combPFIsoHZZ);
  pT->SetBranchAddress("chaPFIso", chaPFIso);
  pT->SetBranchAddress("neuPFIso", neuPFIso);
  pT->SetBranchAddress("phoPFIso", phoPFIso);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("vertices", &vertices);
  if(!TString(file).Contains("fake")) pT->SetBranchAddress("mass", &mass);
  else mass=-1.0;

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  // the new WPs with full isolation for triggering electrons
  Int_t WP95trg, WP90trg, WP85trg, WP80trg, WP70trg, WP80x70trg;
  // the new WPs with charged only isolation for non triggering electrons
  Int_t WP95notrg, WP90notrg, WP85notrg, WP80notrg, WP70notrg;
  // the hww2011 WP and the one with the same efficiency
  // hzz WP is using the MVA for the unbiased electrons
  Int_t hwwWP, newhwwWP;
  // first 4 variables needed for TP
  fT->Branch("mass", &mass, "mass/F");
  fT->Branch("pt", &pt, "pt/F");
  fT->Branch("abseta", &abseta, "abseta/F");
  fT->Branch("vertices", &vertices, "vertices/F");
  // for triggering eles
  fT->Branch("wp95trg", &WP95trg, "wp95trg/I");
  fT->Branch("wp90trg", &WP90trg, "wp90trg/I");
  fT->Branch("wp85trg", &WP85trg, "wp85trg/I");
  fT->Branch("wp80trg", &WP80trg, "wp80trg/I");
  fT->Branch("wp70trg", &WP70trg, "wp70trg/I");
  fT->Branch("wp80x70trg", &WP80x70trg, "wp80x70trg/I"); // 2012 WP?
  fT->Branch("hwwWP", &hwwWP, "hwwWP/I");  // 2011 WP
  // for non triggering eles
  fT->Branch("wp95notrg", &WP95notrg, "chwp95notrg/I");
  fT->Branch("wp90notrg", &WP90notrg, "chwp90notrg/I");
  fT->Branch("wp85notrg", &WP85notrg, "chwp85notrg/I");
  fT->Branch("wp80notrg", &WP80notrg, "chwp80notrg/I");
  fT->Branch("wp70notrg", &WP70notrg, "chwp70notrg/I");

  HZZEleIDSelector aSel;

  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     abseta = fabs(eta);

     hwwWP=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kBDTHWW2011_withIP) && 
	passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kIsoHWW2011)) hwwWP = 1;

     WP95trg=WP90trg=WP85trg=WP80trg=WP70trg=WP80x70trg=0;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVABiased)) WP95trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP90,HZZEleIDSelector::kMVABiased)) WP90trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP85,HZZEleIDSelector::kMVABiased)) WP85trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVABiased)) WP80trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVABiased)) WP70trg=1;
     // mixed 70 x 80
     if(pt<20 && aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVABiased)) WP80x70trg=1;
     if(pt>=20 && aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVABiased)) WP80x70trg=1;

     WP95notrg=WP90notrg=WP85notrg=WP80notrg=WP70notrg=0;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVAUnbiased)) WP95notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP90,HZZEleIDSelector::kMVAUnbiased)) WP90notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP85,HZZEleIDSelector::kMVAUnbiased)) WP85notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVAUnbiased)) WP80notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVAUnbiased)) WP70notrg=1;

     fT->Fill();
  }
  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF << endl;
}


int main(int argc, char* argv[]) {

  char files1[500], files2[500], fileb1[500], fileb2[500], fileb3[500];
  sprintf(files1,"macro/results_data/electrons.root");
  sprintf(files2,"macro/results_data/electrons_zeemc.root");
  sprintf(fileb1,"macro/results_data/fakes.root");
  sprintf(fileb2,"macro/results_data/fakes-unbiased-wlnu.root");
  sprintf(fileb3,"macro/results_data/fakes-zeeOneFake.root");

  cout << "\t===> DOING ISOLATION FRIEND TREES <===" << endl;
  // isolation
  makeFriendHZZIsolation(files1);
  makeFriendHZZIsolation(files2);
  makeFriendHZZIsolation(fileb1);
  makeFriendHZZIsolation(fileb2);
  makeFriendHZZIsolation(fileb3);

  cout << "\t===> DOING ID FRIEND TREES <===" << endl;
  // id bits
  makeFriendHZZIdBits(files1);
  makeFriendHZZIdBits(files2);
  makeFriendHZZIdBits(fileb1);
  makeFriendHZZIdBits(fileb2);
  makeFriendHZZIdBits(fileb3);
  
}
