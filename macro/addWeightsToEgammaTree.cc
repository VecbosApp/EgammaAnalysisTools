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
  float deltaEta;
  float deltaEtaCorr;
  float deltaPhi;
  float deltaPhiCorr;
  float sigmaIEtaIEta;
  float fBrem;
  float convDcot;
  float convDist;
  int missingHits;
  float qcdEta;
  float qcdPt;
  float qcdDeltaphi;
  float qcdMet;
  int qcdNBrem;
  float absTrackerIsolGammaCand;
  float absEcalIsolGammaCand;
  float absHcalIsolGammaCand;
  int isGamma;

  treeOrig->SetBranchAddress("EoPout", &EoPout);
  treeOrig->SetBranchAddress("EoP", &EoP);
  treeOrig->SetBranchAddress("HoE", &HoE);
  treeOrig->SetBranchAddress("deltaEta", &deltaEta);
  treeOrig->SetBranchAddress("deltaEtaCorr", &deltaEtaCorr);
  treeOrig->SetBranchAddress("deltaPhi", &deltaPhi);
  treeOrig->SetBranchAddress("deltaPhiCorr", &deltaPhiCorr);
  treeOrig->SetBranchAddress("sigmaIEtaIEta", &sigmaIEtaIEta);
  treeOrig->SetBranchAddress("fBrem", &fBrem);
  treeOrig->SetBranchAddress("convDcot", &convDcot);
  treeOrig->SetBranchAddress("convDist", &convDist);
  treeOrig->SetBranchAddress("missingHits", &missingHits);
  treeOrig->SetBranchAddress("qcdEta", &qcdEta);
  treeOrig->SetBranchAddress("qcdPt", &qcdPt);
  treeOrig->SetBranchAddress("qcdDeltaphi", &qcdDeltaphi);
  treeOrig->SetBranchAddress("qcdMet", &qcdMet);
  treeOrig->SetBranchAddress("qcdNBrem", &qcdNBrem);
  treeOrig->SetBranchAddress("absTrackerIsolGammaCand", &absTrackerIsolGammaCand);
  treeOrig->SetBranchAddress("absEcalIsolGammaCand", &absEcalIsolGammaCand);
  treeOrig->SetBranchAddress("absHcalIsolGammaCand", &absHcalIsolGammaCand);
  treeOrig->SetBranchAddress("isGamma", &isGamma);

  // copy branches
  treeNew->Branch("EoPout", &EoPout, "EoPout/F");
  treeNew->Branch("EoP", &EoP, "EoP/F");
  treeNew->Branch("HoE", &HoE, "HoE/F");
  treeNew->Branch("deltaEta",     &deltaEta, "deltaEta/F");
  treeNew->Branch("deltaEtaCorr", &deltaEtaCorr, "deltaEtaCorr/F");
  treeNew->Branch("deltaPhi",     &deltaPhi, "deltaPhi/F");
  treeNew->Branch("deltaPhiCorr", &deltaPhiCorr, "deltaPhiCorr/F");
  treeNew->Branch("sigmaIEtaIEta", &sigmaIEtaIEta, "sigmaIEtaIEta/F");
  treeNew->Branch("fBrem",    &fBrem, "fBrem/F");
  treeNew->Branch("convDcot", &convDcot, "convDcot/F");
  treeNew->Branch("convDist", &convDist, "convDist/F");
  treeNew->Branch("missingHits", &missingHits, "missingHits/I");
  treeNew->Branch("qcdEta", &qcdEta, "qcdEta/F");
  treeNew->Branch("qcdPt",  &qcdPt, "qcdPt/F");
  treeNew->Branch("qcdDeltaphi", &qcdDeltaphi, "qcdDeltaphi/F");
  treeNew->Branch("qcdMet",   &qcdMet, "qcdMet/F");
  treeNew->Branch("qcdNBrem", &qcdNBrem, "qcdNBrem/I");
  treeNew->Branch("absTrackerIsolGammaCand", &absTrackerIsolGammaCand, "absTrackerIsolGammaCand/F");
  treeNew->Branch("absEcalIsolGammaCand", &absEcalIsolGammaCand, "absEcalIsolGammaCand/F");
  treeNew->Branch("absHcalIsolGammaCand", &absHcalIsolGammaCand, "absHcalIsolGammaCand/F");
  treeNew->Branch("isGamma", &isGamma);
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
