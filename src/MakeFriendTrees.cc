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
  kBDTHWW2011_withIP // HWW cuts 2011
};

bool passHWWID(Float_t eta, Float_t pt, Float_t bdthww, Float_t bdthzz, Float_t rho, Float_t combIso, Float_t combPFIsoHWW, idType type) {
  if(type == kIsoHWW2011) {
    if(fabs(eta)<1.479) return (combPFIsoHWW/pt < 0.13);
    else return (combPFIsoHWW/pt < 0.09);
  }

  if(type == kBDTHWW2011_withIP) {
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.884);
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
  Float_t combIso;
  fT->Branch("combPFIsoHZZ",&combIso,"combPFIsoHZZ/F");
  
  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     combIso = chaPFIso+neuPFIso+phoPFIso;
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
  Float_t iso, bdt, eta, pt, rho;
  Float_t mass; // not dummy only for TP trees
  pT->SetBranchAddress("bdthww", &bdthww);
  pT->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
  pT->SetBranchAddress("bdthzz",&bdt);
  pT->SetBranchAddress("combPFIsoHZZ",&iso);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("rho", &rho);
  if(!TString(file).Contains("fake")) pT->SetBranchAddress("mass", &mass);
  else mass=-1.0;

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  // the new WPs with full isolation
  Int_t WP95, WP90, WP85, WP80, WP70;
  // the new WPs with charged only isolation
  Int_t chWP95, chWP90, chWP85;
  // the hww2011 WP
  Int_t hwwWP;
  fT->Branch("mass", &mass, "mass/F");
  fT->Branch("wp95", &WP95, "wp95/I");
  fT->Branch("wp90", &WP90, "wp90/I");
  fT->Branch("wp85", &WP85, "wp85/I");
  fT->Branch("wp80", &WP80, "wp80/I");
  fT->Branch("wp70", &WP70, "wp70/I");
  fT->Branch("chwp95", &chWP95, "chwp95/I");
  fT->Branch("chwp90", &chWP90, "chwp90/I");
  fT->Branch("chwp85", &chWP85, "chwp85/I");
  fT->Branch("bdthww", &hwwWP, "bdthww/I");

  HZZEleIDSelector aSel;

  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     hwwWP=0;
     if(passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kBDTHWW2011_withIP) && 
	passHWWID(eta,pt,bdthww,bdt,rho,iso,combPFIsoHWW,kIsoHWW2011)) hwwWP = 1;

     WP95=WP90=WP85=WP80=WP70=0;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP95)) WP95=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP90)) WP90=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP85)) WP85=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP80)) WP80=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP70)) WP70=1;

     chWP95=chWP90=chWP85=0;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP95ChIso)) chWP95=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP90ChIso)) chWP90=1;
     if(aSel.output(pt,eta,bdt,iso,HZZEleIDSelector::kWP85ChIso)) chWP85=1;
     
     fT->Fill();
  }
  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF << endl;
}


int main(int argc, char* argv[]) {

  TFile *fileSig, *fileBkg;
  fileSig = fileBkg = 0;

  char files[500], fileb[500];
  sprintf(files,"macro/results_data/merged.root");
  sprintf(fileb,"macro/results_data_fakes/merged.root");

  cout << "\t===> DOING ISOLATION FRIEND TREES <===" << endl;
  // isolation
  makeFriendHZZIsolation(files);
  makeFriendHZZIsolation(fileb);

  cout << "\t===> DOING ID FRIEND TREES <===" << endl;
  // id bits
  makeFriendHZZIdBits(files);
  makeFriendHZZIdBits(fileb);
  
}
