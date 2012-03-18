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

#include "include/HZZEleIDSelector.hh"

using namespace std;


enum idType {
  kIsoHWW2011 = 0, // HWW cuts 2011
  kBDTHWW2011_withIP, // HWW cuts 2011
  kIso,
  kIsoEACorr,
  kBDTHZZ_withIP
};

bool passHWWID(Float_t eta, Float_t pt, Float_t bdthww, Float_t bdthzz, Float_t rho, Float_t combIso, Float_t combPFIsoHWW, idType type) {
  if(type == kIsoHWW2011) {
    if(fabs(eta)<1.479) return (combPFIsoHWW/pt < 0.13);
    else return (combPFIsoHWW/pt < 0.09);
  }

  if(type == kIsoEACorr) {
    // WP with same fake rate as HWW with IP
    if(fabs(eta) <  1.0) combIso -= 0.18 * rho;
    if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) combIso -= 0.19 * rho;
    if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) combIso -= 0.21 * rho;
    if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) combIso -= 0.38 * rho;
    if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) combIso -= 0.61 * rho;
    if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) combIso -= 0.73 * rho;
    if(fabs(eta) >=  2.4) combIso -= 0.78 * rho;

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

  Float_t chaPFIso, neuPFIso, phoPFIso, rho, eta;
  pT->SetBranchAddress("chaPFIso", &chaPFIso);
  pT->SetBranchAddress("neuPFIso", &neuPFIso);
  pT->SetBranchAddress("phoPFIso", &phoPFIso);
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

     combIsoNoEA = chaPFIso+neuPFIso+phoPFIso;
     combIso = combIsoNoEA;
     if(fabs(eta) <  1.0) combIso -= 0.18 * rho;
     if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) combIso -= 0.19 * rho;
     if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) combIso -= 0.21 * rho;
     if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) combIso -= 0.38 * rho;
     if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) combIso -= 0.61 * rho;
     if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) combIso -= 0.73 * rho;
     if(fabs(eta) >=  2.4) combIso -= 0.78 * rho;
     fT->Fill();
  }

  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF << endl;

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

  Float_t bdthww, combPFIsoHWW;
  Float_t iso, bdt, eta, abseta, pt, rho, vertices;
  Float_t mass; // not dummy only for TP trees
  pT->SetBranchAddress("bdthww", &bdthww);
  pT->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
  pT->SetBranchAddress("bdthzz",&bdt);
  pT->SetBranchAddress("combPFIsoHZZ",&iso);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("vertices", &vertices);
  if(!TString(file).Contains("fake")) pT->SetBranchAddress("mass", &mass);
  else mass=-1.0;

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  // the new WPs with full isolation
  Int_t WP95, WP90, WP85, WP80, WP70, WP80x70;
  // the new WPs with charged only isolation
  Int_t chWP95, chWP90, chWP85, chWP80, chWP70;
  // the hww2011 WP and the one with the same fake rate
  Int_t hwwWP, hzzwphww;
  // isolation only, three bits
  Int_t pfisohww, pfisohzz, pfisohzzEA; 
  // first 4 variables needed for TP
  fT->Branch("mass", &mass, "mass/F");
  fT->Branch("pt", &pt, "pt/F");
  fT->Branch("abseta", &abseta, "abseta/F");
  fT->Branch("vertices", &vertices, "vertices/F");
  fT->Branch("wp95", &WP95, "wp95/I");
  fT->Branch("wp90", &WP90, "wp90/I");
  fT->Branch("wp85", &WP85, "wp85/I");
  fT->Branch("wp80", &WP80, "wp80/I");
  fT->Branch("wp70", &WP70, "wp70/I");
  fT->Branch("wp70x80", &WP80x70, "wp80x70/I");
  fT->Branch("chwp95", &chWP95, "chwp95/I");
  fT->Branch("chwp90", &chWP90, "chwp90/I");
  fT->Branch("chwp85", &chWP85, "chwp85/I");
  fT->Branch("chwp80", &chWP80, "chwp80/I");
  fT->Branch("chwp70", &chWP70, "chwp70/I");
  fT->Branch("bdthww", &hwwWP, "bdthww/I");
  fT->Branch("hzzwphww", &hzzwphww, "hzzwphww/I");
  fT->Branch("pfisohww", &pfisohww, "pfisohww/I");
  fT->Branch("pfisohzz", &pfisohzz, "pfisohzz/I");
  fT->Branch("pfisohzzEA", &pfisohzzEA, "pfisohzzEA/I");

  HZZEleIDSelector aSel;

  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     abseta = fabs(eta);

     hwwWP=0;
     if(passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kBDTHWW2011_withIP) && 
	passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kIsoHWW2011)) hwwWP = 1;

     WP95=WP90=WP85=WP80=WP70=WP80x70=0;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP95)) WP95=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP90)) WP90=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP85)) WP85=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP80)) WP80=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP70)) WP70=1;
     // mixed 70 x 80
     if(pt<20 && aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP70)) WP80x70=1;
     if(pt>=20 && aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP80)) WP80x70=1;

     chWP95=chWP90=chWP85=0;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP95ChIso)) chWP95=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP90ChIso)) chWP90=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP85ChIso)) chWP85=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP80ChIso)) chWP80=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP70ChIso)) chWP70=1;
    
     hzzwphww=0;
     if(passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kBDTHZZ_withIP) && 
        passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kIsoEACorr)) hzzwphww = 1;

     pfisohww=pfisohzz=pfisohzzEA=0;
     if(passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kIsoHWW2011)) pfisohww = 1;
     if(passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kIso)) pfisohzz = 1;
     if(passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kIsoEACorr)) pfisohzzEA = 1;

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
