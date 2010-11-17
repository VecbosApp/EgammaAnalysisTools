#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

using namespace std;

void addWeights(const char* filename, float weight) {

  cout << "Adding weight branch to file " << filename << " with weight " << weight << endl;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  if ( treeOrig ) {
  int nentriesOrig = treeOrig->GetEntries();

  TFile *fileNew = TFile::Open(filename,"recreate");
  TTree *treeNew = new TTree("T1","tree with only selected events");

  // Declaration of leaf types
  float EoPout;
  float EoP;
  float HoE;
  float deta;
  float dphi;
  float see;
  float spp;
  float fbrem;
  float dcot;
  float dist;
  int missHits;
  float eta;
  float pt;
  int nbrem;

  treeOrig->SetBranchAddress("EoPout", &EoPout);
  treeOrig->SetBranchAddress("EoP", &EoP);
  treeOrig->SetBranchAddress("HoE", &HoE);
  treeOrig->SetBranchAddress("deta", &deta);
  treeOrig->SetBranchAddress("dphi", &dphi);
  treeOrig->SetBranchAddress("see", &see);
  treeOrig->SetBranchAddress("spp", &spp);
  treeOrig->SetBranchAddress("fbrem", &fbrem);
  treeOrig->SetBranchAddress("dcot", &dcot);
  treeOrig->SetBranchAddress("dist", &dist);
  treeOrig->SetBranchAddress("missHits", &missHits);
  treeOrig->SetBranchAddress("eta", &eta);
  treeOrig->SetBranchAddress("pt", &pt);
  treeOrig->SetBranchAddress("nbrem", &nbrem);

  // copy branches
  treeNew->Branch("EoPout", &EoPout, "EoPout/F");
  treeNew->Branch("EoP", &EoP, "EoP/F");
  treeNew->Branch("HoE", &HoE, "HoE/F");
  treeNew->Branch("deta",     &deta, "deta/F");
  treeNew->Branch("dphi",     &dphi, "dphi/F");
  treeNew->Branch("see", &see, "see/F");
  treeNew->Branch("spp", &spp, "spp/F");
  treeNew->Branch("fbrem",    &fbrem, "fbrem/F");
  treeNew->Branch("dcot", &dcot, "dcot/F");
  treeNew->Branch("dist", &dist, "dist/F");
  treeNew->Branch("missHits", &missHits, "missHits/I");
  treeNew->Branch("eta", &eta, "eta/F");
  treeNew->Branch("pt",  &pt, "pt/F");
  treeNew->Branch("nbrem", &nbrem, "nbrem/I");
  treeNew->Branch("weight", &weight,  "weight/F");

  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%1000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentriesOrig << " entries" << std::endl;
    treeOrig->GetEntry(i);
    treeNew->Fill();
  }
  
  fileNew->cd();
  treeNew->Write();
  fileNew->Close();
  
  fileOrig->cd();
  fileOrig->Close();
  
  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }
}
