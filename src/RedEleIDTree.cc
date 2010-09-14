#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("EoPout",          &myEoPout,          "EoPout/F");
  myTree->Branch("EoP",             &myEoP,             "EoP/F");
  myTree->Branch("HoE",             &myHoE,             "HoE/F");
  myTree->Branch("deta",            &myDeta,            "deta/F");
  myTree->Branch("dphi",            &myDphi,            "dphi/F");
  myTree->Branch("detaUncorr",      &myDetaUncorr,      "detaUncorr/F");
  myTree->Branch("dphiUncorr",      &myDphiUncorr,      "dphiUncorr/F");
  myTree->Branch("s9s25",           &mys9s25,           "s9s25/F");
  myTree->Branch("s1s9" ,           &mys1s9,            "s1s9/F");
  myTree->Branch("see",             &mySee,             "see/F");
  myTree->Branch("spp",             &mySpp,             "spp/F");
  myTree->Branch("fbrem",           &myFbrem,           "fbrem/F");
  myTree->Branch("missHits",        &myMissHits,        "missHits/I");
  myTree->Branch("dist",            &myDist,            "dist/F");
  myTree->Branch("dcot",            &myDcot,            "dcot/F");
  myTree->Branch("pt",              &myPt,              "pt/F");
  myTree->Branch("eta",             &myEta,             "eta/F");
  myTree->Branch("charge",          &myCharge,          "charge/I");
}

RedEleIDTree::~RedEleIDTree() {
  delete myFile;
}

void RedEleIDTree::addAttributesSignal() {

  myTree->Branch("zmass",       &myZmass,       "zmass/F");
}

void RedEleIDTree::addAttributesBackground() {

  myTree->Branch("qcdDeltaphi",    &myQCDDeltaphi,    "qcdDeltaphi/F");
  myTree->Branch("qcdInvmass",     &myQCDInvmass,     "qcdInvmass/F");
  myTree->Branch("qcdMet",         &myQCDMet,         "qcdMet/F");
  myTree->Branch("qcdPtHat",       &myQCDPtHat,       "qcdPtHat/F");
}

void RedEleIDTree::addCategories() {

  myTree->Branch("iecal",      &myiecal,      "iecal/I");
  myTree->Branch("iptbin",     &myiptbin,     "iptbin/I");
  myTree->Branch("iclass",     &myiclass,     "iclass/I");
  myTree->Branch("nbrem",      &mynbrem,      "nbrem/I");
}

void RedEleIDTree::addMore() {
  
  myTree->Branch("relIsolTag",   &myRelIsolTag,   "relIsolTag/F");
  myTree->Branch("relIsolProbe", &myRelIsolProbe, "relIsolProbe/F");
}

void RedEleIDTree::addGamma() {

  myTree->Branch("absTrackerIsolGammaCand",&myAbsTrackerIsolGammaCand,"absTrackerIsolGammaCand/F");
  myTree->Branch("absEcalIsolGammaCand",   &myAbsEcalIsolGammaCand,   "absEcalIsolGammaCand/F");
  myTree->Branch("absHcalIsolGammaCand",   &myAbsHcalIsolGammaCand,   "absHcalIsolGammaCand/F");
  myTree->Branch("isGamma",                &myIsGamma,                "isGamma/I");
}

void RedEleIDTree::store() {

  myTree->Fill();
}

void RedEleIDTree::save() {

  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedEleIDTree::fillVariables(float EoPout, float EoP, float HoE, float Deta, float Dphi, float s9s25, float s1s9, float See, float Spp, float pt, float eta, int charge) {

  myEoPout=EoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeta=Deta;
  myDphi=Dphi;
  mys9s25=s9s25;
  mys1s9=s1s9;
  mySee=See;
  mySpp=Spp;
  myPt=pt;
  myEta=eta;
  myCharge=charge;
}

void RedEleIDTree::fillVariables(float EoPout, float EoP, float HoE, float Deta, float DetaUncorr, float Dphi, float DphiUncorr, float s9s25, float s1s9, float See, float Spp, float fbrem, int nHits, float dcot, float dist, float pt, float eta, int charge) {

  myEoPout=EoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeta=Deta;
  myDetaUncorr=DetaUncorr;
  myDphi=Dphi;
  myDphiUncorr=DphiUncorr;
  mys9s25=s9s25;
  mys1s9=s1s9;
  mySee=See;
  mySpp=Spp;
  myFbrem=fbrem;
  myMissHits=nHits;
  myDist=dist;
  myDcot=dcot;
  myPt=pt;
  myEta=eta;
  myCharge=charge;
}

void RedEleIDTree::fillAttributesSignal(float zmass) {

  myZmass=zmass;
}

void RedEleIDTree::fillAttributesBackground(float deltaphi, float invmass, float met, float pth) {

  myQCDDeltaphi=deltaphi;
  myQCDInvmass=invmass;
  myQCDMet=met;
  myQCDPtHat=pth;
}

void RedEleIDTree::fillCategories(int iecal, int iptbin, int iclass, int nbr) {

  myiecal=iecal;
  myiptbin=iptbin;
  myiclass=iclass;
  mynbrem=nbr;
}

void RedEleIDTree::fillMore(float rit, float rip) {

  myRelIsolTag=rit;
  myRelIsolProbe=rip;
}

void RedEleIDTree::fillGamma(float atg, float aeg, float ahg, int ig) {

  myAbsTrackerIsolGammaCand = atg;
  myAbsEcalIsolGammaCand    = aeg;
  myAbsHcalIsolGammaCand    = ahg;
  myIsGamma=ig;
}

