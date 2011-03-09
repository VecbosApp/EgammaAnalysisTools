#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("EoPout",          &myEoPout,          "EoPout/F");
  myTree->Branch("EoP",             &myEoP,             "EoP/F");
  myTree->Branch("HoE",             &myHoE,             "HoE/F");
  myTree->Branch("deta",            &myDeta,            "deta/F");
  myTree->Branch("dphi",            &myDphi,            "dphi/F");
  myTree->Branch("s9s25",           &mys9s25,           "s9s25/F");
  myTree->Branch("see",             &mySee,             "see/F");
  myTree->Branch("spp",             &mySpp,             "spp/F");
  myTree->Branch("fbrem",           &myFbrem,           "fbrem/F");
  myTree->Branch("nbrem",           &myNbrems,          "nbrems/I");
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

void RedEleIDTree::addElectronIdBits() {

  myTree->Branch("CutBasedId",         myCutBasedId,         "CutBasedId[4]/I");
  myTree->Branch("CutBasedIdOlyID",    myCutBasedIdOnlyID,   "CutBasedIdOnlyID[4]/I");
  myTree->Branch("CutBasedIdOnlyIso",  myCutBasedIdOnlyIso,  "CutBasedIdOnlyIso[4]/I");
  myTree->Branch("CutBasedIdOnlyConv", myCutBasedIdOnlyConv, "CutBasedIdOnlyConv[4]/I");
  myTree->Branch("LHBasedId",          myLHBasedId,          "LHBasedId[5]/I");
  myTree->Branch("LHBasedIdOlyID",     myLHBasedIdOnlyID,    "LHBasedIdOnlyID[5]/I");
  myTree->Branch("LHBasedIdOnlyIso",   myLHBasedIdOnlyIso,   "LHBasedIdOnlyIso[5]/I");
  myTree->Branch("LHBasedIdOnlyConv",  myLHBasedIdOnlyConv,  "LHBasedIdOnlyConv[5]/I");
  myTree->Branch("CiCBasedId",         myCiCBasedId,         "CiCBasedId[9]/I");
  myTree->Branch("CiCBasedIdOlyID",    myCiCBasedIdOnlyID,   "CiCBasedIdOnlyID[9]/I");
  myTree->Branch("CiCBasedIdOnlyIso",  myCiCBasedIdOnlyIso,  "CiCBasedIdOnlyIso[9]/I");
  myTree->Branch("CiCBasedIdOnlyConv", myCiCBasedIdOnlyConv, "CiCBasedIdOnlyConv[9]/I");

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

void RedEleIDTree::addIsolations() {
  myTree->Branch("trkIso",  &myTrkIso,    "trkIso/F");
  myTree->Branch("ecalIso", &myEcalIso,   "ecalIso/F");
  myTree->Branch("hcalIso", &myHcalIso,   "hcalIso/F");
}

void RedEleIDTree::addMore() {
  myTree->Branch("nVTx",   &myNVtx,   "nVtx/I");
  myTree->Branch("rho",    &myRho,    "rho/F");
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

void RedEleIDTree::fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float See, float Spp, float fbrem, int nbrems, float pt, float eta, int charge) {
  myEoPout=EoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeta=DEta;
  myDphi=DPhi;
  mys9s25=s9s25;
  mySee=See;
  mySpp=Spp;
  myFbrem=fbrem;
  myNbrems=nbrems;
  myPt=pt;
  myEta=eta;
  myCharge=charge;
}

void RedEleIDTree::fillVariables(float EoPout, float EoP, float HoE, float Deta, float Dphi, float s9s25, float s1s9, float See, float Spp, float fbrem, int nHits, float dcot, float dist, float pt, float eta, int charge) {

  myEoPout=EoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeta=Deta;
  myDphi=Dphi;
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

void RedEleIDTree::fillIsolations(float trkIso, float ecalIso, float hcalIso) {
  myTrkIso=trkIso;
  myEcalIso=ecalIso;
  myHcalIso=hcalIso;
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

void RedEleIDTree::fillMore(int nVtx, float rho) {
  myNVtx=nVtx;
  myRho=rho;
}

void RedEleIDTree::fillGamma(float atg, float aeg, float ahg, int ig) {

  myAbsTrackerIsolGammaCand = atg;
  myAbsEcalIsolGammaCand    = aeg;
  myAbsHcalIsolGammaCand    = ahg;
  myIsGamma=ig;
}

void RedEleIDTree::fillCutBasedIDBits(int CutBasedId[4], int CutBasedIdOnlyID[4], int CutBasedIdOnlyIso[4], int CutBasedIdOnlyConv[4]) {
  for(int i=0; i<4; i++) {
    myCutBasedId[i] = CutBasedId[i];
    myCutBasedIdOnlyID[i] = CutBasedIdOnlyID[i];
    myCutBasedIdOnlyIso[i] = CutBasedIdOnlyIso[i];
    myCutBasedIdOnlyConv[i] = CutBasedIdOnlyConv[i];
  }
}

void RedEleIDTree::fillLHBasedIDBits(int LHBasedId[4], int LHBasedIdOnlyID[4], int LHBasedIdOnlyIso[4], int LHBasedIdOnlyConv[4]) {
  for(int i=0; i<5; i++) {
    myLHBasedId[i] = LHBasedId[i];
    myLHBasedIdOnlyID[i] = LHBasedIdOnlyID[i];
    myLHBasedIdOnlyIso[i] = LHBasedIdOnlyIso[i];
    myLHBasedIdOnlyConv[i] = LHBasedIdOnlyConv[i];
  }
}

void RedEleIDTree::fillCiCBasedIDBits(int CiCBasedId[9], int CiCBasedIdOnlyID[9], int CiCBasedIdOnlyIso[9], int CiCBasedIdOnlyConv[9]) {
  for(int i=0; i<9; i++) {
    myCiCBasedId[i] = CiCBasedId[i];
    myCiCBasedIdOnlyID[i] = CiCBasedIdOnlyID[i];
    myCiCBasedIdOnlyIso[i] = CiCBasedIdOnlyIso[i];
    myCiCBasedIdOnlyConv[i] = CiCBasedIdOnlyConv[i];
  }
}
