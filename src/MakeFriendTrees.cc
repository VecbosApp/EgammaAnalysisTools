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
#include "include/eIDCiChzzSelector.hh"

using namespace std;

enum idType {
  kIsoHWW2011 = 0, // HWW cuts 2011
  kBDTHWW2011_withIP, // HWW cuts 2011
  kIso,
  kIsoEACorr,
  kBDTHZZ_withIP
};

float Aeff_neu_dr04[7], Aeff_pho_dr04[7];

float Aeff_neu_dr04_2011[7] = { 0.045, 0.065, 0.068, 0.057, 0.058, 0.061, 0.110 };
float Aeff_pho_dr04_2011[7] = { 0.140, 0.130, 0.079, 0.130, 0.150, 0.160, 0.180 };

float Aeff_neu_dr04_2012[7] = { 0.041, 0.068, 0.075, 0.068, 0.071, 0.078, 0.140 };
float Aeff_pho_dr04_2012[7] = { 0.144, 0.138, 0.084, 0.155, 0.201, 0.223, 0.265 };

// H->ZZ detector based effective areas
float Aeff_ecal_dr03[2] = { 0.078, 0.046 };
float Aeff_hcal_dr03[2] = { 0.026, 0.072 };

bool cicid(int *cic, int level) { return (cic[level]>>0)%2; }
bool ciciso(int *cic, int level) { return (cic[level]>>1)%2; }
bool cicconv(int *cic, int level) { return (cic[level]>>2)%2; }
bool cicip(int *cic, int level) { return (cic[level]>>3)%2; }

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
  Float_t trkIso,ecalIso,hcalIso;
  pT->SetBranchAddress("chaPFIso", chaPFIso);
  pT->SetBranchAddress("neuPFIso", neuPFIso);
  pT->SetBranchAddress("phoPFIso", phoPFIso);
  pT->SetBranchAddress("trkIso", &trkIso);
  pT->SetBranchAddress("ecalIso", &ecalIso);
  pT->SetBranchAddress("hcalIso", &hcalIso);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("eta", &eta);

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  Float_t combDetIso, combIso, combIsoNoEA;
  fT->Branch("combDetIsoHZZ",&combDetIso,"combDetIsoHZZ/F");
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

     int ieta = (fabs(eta)<1.479) ? 0 : 1;
     combDetIso = trkIso + ecalIso - rho * Aeff_ecal_dr03[ieta] + hcalIso - rho * Aeff_hcal_dr03[ieta];

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
  Int_t DenomFakeSmurf, misshits, ecalseed;
  Float_t eop,eseedopin,HoE,deta,dphi,see,fbrem,dist,dcot,d0,trkIso,ecalIso,hcalIso;
  pT->SetBranchAddress("bdthww", bdthww);
  pT->SetBranchAddress("newbdthww", newbdthww);
  pT->SetBranchAddress("bdthzz",bdthzz);
  pT->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
  pT->SetBranchAddress("combPFIsoHZZ",&combPFIsoHZZ);
  pT->SetBranchAddress("chaPFIso", chaPFIso);
  pT->SetBranchAddress("neuPFIso", neuPFIso);
  pT->SetBranchAddress("phoPFIso", phoPFIso);
  pT->SetBranchAddress("DenomFakeSmurf", &DenomFakeSmurf);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("vertices", &vertices);
  pT->SetBranchAddress("EoP", &eop);
  pT->SetBranchAddress("EseedoPin", &eseedopin);
  pT->SetBranchAddress("HoE", &HoE);
  pT->SetBranchAddress("deta", &deta);
  pT->SetBranchAddress("dphi", &dphi);
  pT->SetBranchAddress("see", &see);
  pT->SetBranchAddress("fbrem", &fbrem);
  pT->SetBranchAddress("missHits", &misshits);
  pT->SetBranchAddress("dist", &dist);
  pT->SetBranchAddress("dcot", &dcot);
  pT->SetBranchAddress("d0",&d0);
  pT->SetBranchAddress("ecaldriven", &ecalseed);
  pT->SetBranchAddress("trkIso", &trkIso);
  pT->SetBranchAddress("ecalIso", &ecalIso);
  pT->SetBranchAddress("hcalIso", &hcalIso);
  if(!TString(file).Contains("fake")) pT->SetBranchAddress("mass", &mass);
  else mass=-1.0;

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  // the new WPs with full isolation for triggering electrons
  Int_t WP95trg, WP90trg, WP85trg, WP80trg, WP70trg;
  // the new WPs with charged only isolation for non triggering electrons
  Int_t WP95notrg, WP90notrg, WP85notrg, WP80notrg, WP70notrg;
  // the hww2011 WP and the one with the same efficiency
  // hzz WP is using the MVA for the unbiased electrons
  Int_t hwwWP, newhwwWP, newhzzWP;
  Int_t hwwWPisoonly, newhwwWPisoonly, newhzzWPisoonly;
  Int_t hwwWPidonly, newhwwWPidonly, newhzzWPidonly;
  Int_t cicmedium, cicmediumid, cicmediumiso;
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
  fT->Branch("newhwwWP", &newhwwWP, "newhwwWP/I"); // 2012 WP?
  fT->Branch("newhwwWPisoonly", &newhwwWPisoonly, "newhwwWPisoonly/I"); // 2012 WP?
  fT->Branch("newhwwWPidonly", &newhwwWPidonly, "newhwwWPidonly/I"); // 2012 WP?
  fT->Branch("hwwWP", &hwwWP, "hwwWP/I");  // 2011 WP
  fT->Branch("hwwWPisoonly", &hwwWPisoonly, "hwwWPisoonly/I");  // 2011 WP
  fT->Branch("hwwWPidonly", &hwwWPidonly, "hwwWPidonly/I");  // 2011 WP
  // for non triggering eles
  fT->Branch("wp95notrg", &WP95notrg, "chwp95notrg/I");
  fT->Branch("wp90notrg", &WP90notrg, "chwp90notrg/I");
  fT->Branch("wp85notrg", &WP85notrg, "chwp85notrg/I");
  fT->Branch("wp80notrg", &WP80notrg, "chwp80notrg/I");
  fT->Branch("wp70notrg", &WP70notrg, "chwp70notrg/I");
  // same as HWW DenomFakeSmurf: change name for the friend tree
  fT->Branch("denom", &DenomFakeSmurf, "denom/I");
  // the cic used for HZZ
  fT->Branch("cicmedium", &cicmedium, "cicmedium/I");
  fT->Branch("cicmediumid", &cicmediumid, "cicmediumid/I");
  fT->Branch("cicmediumiso", &cicmediumiso, "cicmediumiso/I");
  fT->Branch("newhzzWP", &newhzzWP, "newhzzWP/I"); // 2012 WP?
  fT->Branch("newhzzWPisoonly", &newhzzWPisoonly, "newhzzWPisoonly/I"); // 2012 WP?
  fT->Branch("newhzzWPidonly", &newhzzWPidonly, "newhzzWPidonly/I"); // 2012 WP?

  HZZEleIDSelector aSel;

  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     abseta = fabs(eta);

     hwwWP=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kBDTHWW2011_withIP) && 
	passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kIsoHWW2011)) hwwWP = 1;
     hwwWPisoonly=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kIsoHWW2011)) hwwWPisoonly = 1;
     hwwWPidonly=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kBDTHWW2011_withIP)) hwwWPidonly = 1;

     WP95trg=WP90trg=WP85trg=WP80trg=WP70trg=0;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVABiased)) WP95trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP90,HZZEleIDSelector::kMVABiased)) WP90trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP85,HZZEleIDSelector::kMVABiased)) WP85trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVABiased)) WP80trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVABiased)) WP70trg=1;
     // same efficiency as 2011 WP
     newhwwWP=0;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWPHWW,HZZEleIDSelector::kMVABiased)) newhwwWP=1;
     newhwwWPisoonly=0;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWPHWW,HZZEleIDSelector::kMVABiased,HZZEleIDSelector::isoonly)) newhwwWPisoonly=1;
     newhwwWPidonly=0;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWPHWW,HZZEleIDSelector::kMVABiased,HZZEleIDSelector::idonly)) newhwwWPidonly=1;

     WP95notrg=WP90notrg=WP85notrg=WP80notrg=WP70notrg=0;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVAUnbiased)) WP95notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP90,HZZEleIDSelector::kMVAUnbiased)) WP90notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP85,HZZEleIDSelector::kMVAUnbiased)) WP85notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVAUnbiased)) WP80notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVAUnbiased)) WP70notrg=1;

     // CiC...
     eIDCiChzzSelector cicsel;
     int cic[5];
     int ieta = (fabs(eta)<1.479) ? 0 : 1;
     for(int i=0; i<5; i++) cic[i] = cicsel.ElectronId_V03(pt,eta,see,eop,eseedopin,fbrem,
							   trkIso,
							   ecalIso - rho * Aeff_ecal_dr03[ieta],
							   hcalIso - rho * Aeff_hcal_dr03[ieta],
							   d0,misshits,deta,dphi,HoE,dcot,
							   dist,!ecalseed,i,false);

     cicmedium=cicmediumid=cicmediumiso=0;
     if(cicid(cic,3) && ciciso(cic,3) && misshits<=1) cicmedium=1;
     if(cicid(cic,3) && misshits<=1) cicmediumid=1;
     if(ciciso(cic,3)) cicmediumiso=1;
     // same fake rate in Z+1 fake as 2011 CIC WP
     newhzzWP=0;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWPHZZ,HZZEleIDSelector::kMVAUnbiased)) newhzzWP=1;
     newhzzWPisoonly=0;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWPHZZ,HZZEleIDSelector::kMVAUnbiased,HZZEleIDSelector::isoonly)) newhzzWPisoonly=1;
     newhzzWPidonly=0;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWPHZZ,HZZEleIDSelector::kMVAUnbiased,HZZEleIDSelector::idonly)) newhzzWPidonly=1;
     
     fT->Fill();
  }
  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF << endl;
}


int main(int argc, char* argv[]) {

  int year = 2012;
  for(int i=0; i<7; i++) {
    if(year==2011) {
      Aeff_neu_dr04[i] = Aeff_neu_dr04_2011[i];
      Aeff_pho_dr04[i] = Aeff_pho_dr04_2011[i];
    } else if(year==2012) {
      Aeff_neu_dr04[i] = Aeff_neu_dr04_2012[i];
      Aeff_pho_dr04[i] = Aeff_pho_dr04_2012[i];
    } else {
      cout << "wrong year set. Returning 0." << endl;
      return 0;
    }
  }


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
