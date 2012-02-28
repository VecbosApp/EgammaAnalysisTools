#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>

using namespace std;

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
     if(fabs(eta) >=  2.4) combIso -= 0.90 * rho;
     fT->Fill();
  }

  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF << endl;

}
