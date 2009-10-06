#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("EoPout",        &myEoPout,        "EoPout/F");
  myTree->Branch("EoP",           &myEoP,           "EoP/F");
  myTree->Branch("HoE",           &myHoE,           "HoE/F");
  myTree->Branch("deltaEta",      &myDeltaEta,      "deltaEta/F");
  myTree->Branch("deltaPhi",      &myDeltaPhi,      "deltaPhi/F");
  myTree->Branch("s9s25",         &mys9s25,         "s9s25/F");
  myTree->Branch("s1s9" ,         &mys1s9,          "s1s9/F");
  myTree->Branch("sigmaIEtaIEta", &mySigmaIEtaIEta, "sigmaIEtaIEta/F");

}

RedEleIDTree::~RedEleIDTree() {
  delete myFile;
}

void RedEleIDTree::addAttributesSignal() {
  myTree->Branch("charge",      &myCharge,      "charge/I");
  myTree->Branch("eta",         &myEta,         "eta/F");
  myTree->Branch("pt",          &myPt,          "pt/F");
  myTree->Branch("zmass",       &myZmass,       "zmass/F");
}

void RedEleIDTree::addAttributesBackground() {
  myTree->Branch("qcdCharge",      &myQCDCharge,      "qcdCharge/I");
  myTree->Branch("qcdEta",         &myQCDEta,         "qcdEta/F");
  myTree->Branch("qcdPt",          &myQCDPt,          "qcdPt/F");
  myTree->Branch("qcdDeltaphi",    &myQCDDeltaphi,    "qcdDeltaphi/F");
  myTree->Branch("qcdInvmass",     &myQCDInvmass,     "qcdInvmass/F");
  myTree->Branch("qcdMet",         &myQCDMet,         "qcdMet/F");
}

void RedEleIDTree::addCategories() {
  myTree->Branch("iecal",      &myiecal,      "iecal/I");
  myTree->Branch("iptbin",     &myiptbin,     "iptbin/I");
  myTree->Branch("iclass",     &myiclass,     "iclass/I");
}

void RedEleIDTree::addMore() {
  myTree->Branch("relIsolTag",   &myRelIsolTag,   "relIsolTag/F");
  myTree->Branch("relIsolProbe", &myRelIsolProbe, "relIsolProbe/F");
}

void RedEleIDTree::store() {
  myTree->Fill();
}

void RedEleIDTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedEleIDTree::fillVariables(float EoPout, float EoP, float HoE, float DeltaEta, float DeltaPhi, float s9s25, float s1s9, float SigmaIEtaIEta) {

  myEoPout=EoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeltaEta=DeltaEta;
  myDeltaPhi=DeltaPhi;
  mys9s25=s9s25;
  mys1s9=s1s9;
  mySigmaIEtaIEta=SigmaIEtaIEta;

}

void RedEleIDTree::fillAttributesSignal(int charge, float eta, float pt, float zmass) {

  myCharge=charge;
  myEta=eta;
  myPt=pt;
  myZmass=zmass;

}

void RedEleIDTree::fillAttributesBackground(int charge, float eta, float pt, float deltaphi, float invmass, float met) {

  myQCDCharge=charge;
  myQCDEta=eta;
  myQCDPt=pt;
  myQCDDeltaphi=deltaphi;
  myQCDInvmass=invmass;
  myQCDMet=met;

}

void RedEleIDTree::fillCategories(int iecal, int iptbin, int iclass) {

  myiecal=iecal;
  myiptbin=iptbin;
  myiclass=iclass;

}

void RedEleIDTree::fillMore(float rit, float rip) {
  myRelIsolTag=rit;
  myRelIsolProbe=rip;
}
