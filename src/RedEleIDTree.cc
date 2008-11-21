#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("EoPout",      &myEoPout,      "EoPout/F");
  myTree->Branch("HoE",         &myHoE,         "HoE/F");
  myTree->Branch("deltaEta",    &myDeltaEta,    "deltaEta/F");
  myTree->Branch("deltaPhi",    &myDeltaPhi,    "deltaPhi/F");
  myTree->Branch("s9s25",       &mys9s25,       "s9s25/F");
  myTree->Branch("sigmaEtaEta", &mySigmaEtaEta, "sigmaEtaEta/F");

}

RedEleIDTree::~RedEleIDTree() {
  delete myFile;
}

void RedEleIDTree::addAttributes() {
  myTree->Branch("charge",      &myCharge,      "charge/I");
  myTree->Branch("eta",         &myEta,         "eta/F");
  myTree->Branch("pt",          &myPt,          "pt/F");
  myTree->Branch("zmass",       &myZmass,       "zmass/F");
}

void RedEleIDTree::addCategories() {
  myTree->Branch("iecal",      &myiecal,      "iecal/I");
  myTree->Branch("iptbin",     &myiptbin,     "iptbin/I");
  myTree->Branch("iclass",     &myiclass,     "iclass/I");
}

void RedEleIDTree::store() {
  myTree->Fill();
}

void RedEleIDTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedEleIDTree::fillVariables(float EoPout, float HoE, float DeltaEta, float DeltaPhi, float s9s25, float SigmaEtaEta) {

  myEoPout=EoPout;
  myHoE=HoE;
  myDeltaEta=DeltaEta;
  myDeltaPhi=DeltaPhi;
  mys9s25=s9s25;
  mySigmaEtaEta=SigmaEtaEta;

}

void RedEleIDTree::fillAttributes(int charge, float eta, float pt, float zmass) {

  myCharge=charge;
  myEta=eta;
  myPt=pt;
  myZmass=zmass;

}

void RedEleIDTree::fillCategories(int iecal, int iptbin, int iclass) {

  myiecal=iecal;
  myiptbin=iptbin;
  myiclass=iclass;

}

