#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myFile->mkdir("eleIDdir");
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

  myTree->Branch("mass",       &myZmass,       "mass/F");
}

void RedEleIDTree::addElectronIdBits() {

  myTree->Branch("WP95",         &myCutBasedId[0],         "WP95/I");
  myTree->Branch("WP90",         &myCutBasedId[1],         "WP90/I");
  myTree->Branch("WP85",         &myCutBasedId[2],         "WP85/I");
  myTree->Branch("WP80",         &myCutBasedId[3],         "WP80/I");
  myTree->Branch("WP70",         &myCutBasedId[4],         "WP70/I");
  myTree->Branch("WPSmurf",      &myCutBasedId[5],         "WPSmurf/I");
  myTree->Branch("CutBasedIdOlyID",    myCutBasedIdOnlyID,   "CutBasedIdOnlyID[6]/I");
  myTree->Branch("CutBasedIdOnlyIso",  myCutBasedIdOnlyIso,  "CutBasedIdOnlyIso[6]/I");
  myTree->Branch("CutBasedIdOnlyConv", myCutBasedIdOnlyConv, "CutBasedIdOnlyConv[6]/I");

  myTree->Branch("LHVeryLoose",  &myLHBasedId[0],            "LHVeryLoose/I");
  myTree->Branch("LHLoose",      &myLHBasedId[1],            "LHLoose/I");
  myTree->Branch("LHMedium",     &myLHBasedId[2],            "LHMedium/I");
  myTree->Branch("LHTight",      &myLHBasedId[3],            "LHTight/I");
  myTree->Branch("LHHyperTight", &myLHBasedId[4],            "LHHyperTight/I");
  myTree->Branch("LHBasedIdOlyID",     myLHBasedIdOnlyID,    "LHBasedIdOnlyID[5]/I");
  myTree->Branch("LHBasedIdOnlyIso",   myLHBasedIdOnlyIso,   "LHBasedIdOnlyIso[5]/I");
  myTree->Branch("LHBasedIdOnlyConv",  myLHBasedIdOnlyConv,  "LHBasedIdOnlyConv[5]/I");

  myTree->Branch("LHPFIsoVeryLoose",  &myLHBasedPFIsoId[0],            "LHPFIsoVeryLoose/I");
  myTree->Branch("LHPFIsoLoose",      &myLHBasedPFIsoId[1],            "LHPFIsoLoose/I");
  myTree->Branch("LHPFIsoMedium",     &myLHBasedPFIsoId[2],            "LHPFIsoMedium/I");
  myTree->Branch("LHPFIsoTight",      &myLHBasedPFIsoId[3],            "LHPFIsoTight/I");
  myTree->Branch("LHPFIsoHyperTight", &myLHBasedPFIsoId[4],            "LHPFIsoHyperTight/I");
  myTree->Branch("LHPFIsoBasedIdOlyID",     myLHBasedPFIsoIdOnlyID,    "LHPFIsoBasedIdOnlyID[5]/I");
  myTree->Branch("LHPFIsoBasedIdOnlyIso",   myLHBasedPFIsoIdOnlyIso,   "LHPFIsoBasedIdOnlyIso[5]/I");
  myTree->Branch("LHPFIsoBasedIdOnlyConv",  myLHBasedPFIsoIdOnlyConv,  "LHPFIsoBasedIdOnlyConv[5]/I");

  myTree->Branch("DenomFake",          &myDenomFake,             "DenomFake/I");
  myTree->Branch("DenomFakeSmurf",     &myDenomFakeSmurf,        "DenomFakeSmurf/I");
  myTree->Branch("BDTIdOnlyId",        &myBDTIdOnlyId,           "BDTIdOnlyId/I");
}

void RedEleIDTree::addDenominatorFakeBits() {
  myTree->Branch("DenomFake",          &myDenomFake,             "DenomFake/I");
  myTree->Branch("DenomFakeSmurf",     &myDenomFakeSmurf,        "DenomFakeSmurf/I");
}

void RedEleIDTree::addRunInfos() {
  myTree->Branch("run",     &myRun,     "run/I");
  myTree->Branch("lumi",    &myLS,      "lumi/I");
  myTree->Branch("event",   &myEvent,   "event/I");
  myTree->Branch("npu",      myNpu,     "npu[3]/I");
  myTree->Branch("mcmatch", &myMCMatch, "mcmatch/I");
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
  myTree->Branch("combPFIsoHWW", &myPFCandCombinedIsoHWW, "combPFIsoHWW/F");
  myTree->Branch("chaPFIso",     &myPFCandChargedIso,     "chPFIso/F");
  myTree->Branch("neuPFIso",     &myPFCandNeutralIso,     "neuPFIso/F");
  myTree->Branch("phoPFIso",     &myPFCandPhotonIso,      "phoPFIso/F");
}

void RedEleIDTree::addMore() {
  myTree->Branch("bdthww",   &myBdtHww, "bdthww/F");
  myTree->Branch("bdthzz",   &myBdtHzz, "bdthzz/F");
  myTree->Branch("vertices", &myNVtx, "vertices/F");
  myTree->Branch("rho",      &myRho,  "rho/F");
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

  myFile->cd("eleIDdir");
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

void RedEleIDTree::fillIsolations(float trkIso, float ecalIso, float hcalIso,
                                  float combPFiso,
                                  float chaPFiso, float neuPFiso, float phoPFiso) {
  myTrkIso=trkIso;
  myEcalIso=ecalIso;
  myHcalIso=hcalIso;
  myPFCandCombinedIsoHWW=combPFiso;
  myPFCandChargedIso=chaPFiso;
  myPFCandNeutralIso=neuPFiso;
  myPFCandPhotonIso=phoPFiso;
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

void RedEleIDTree::fillMore(float nVtx, float rho, float bdthww, float bdthzz) {
  myNVtx=nVtx;
  myRho=rho;
  myBdtHww=bdthww;
  myBdtHzz=bdthzz;
}

void RedEleIDTree::fillGamma(float atg, float aeg, float ahg, int ig) {

  myAbsTrackerIsolGammaCand = atg;
  myAbsEcalIsolGammaCand    = aeg;
  myAbsHcalIsolGammaCand    = ahg;
  myIsGamma=ig;
}

void RedEleIDTree::fillCutBasedIDBits(int CutBasedId[6], int CutBasedIdOnlyID[6], int CutBasedIdOnlyIso[6], int CutBasedIdOnlyConv[6]) {
  for(int i=0; i<6; i++) {
    myCutBasedId[i] = CutBasedId[i];
    myCutBasedIdOnlyID[i] = CutBasedIdOnlyID[i];
    myCutBasedIdOnlyIso[i] = CutBasedIdOnlyIso[i];
    myCutBasedIdOnlyConv[i] = CutBasedIdOnlyConv[i];
  }
}

void RedEleIDTree::fillLHBasedIDBits(int LHBasedId[5], int LHBasedIdOnlyID[5], int LHBasedIdOnlyIso[5], int LHBasedIdOnlyConv[5]) {
  for(int i=0; i<5; i++) {
    myLHBasedId[i] = LHBasedId[i];
    myLHBasedIdOnlyID[i] = LHBasedIdOnlyID[i];
    myLHBasedIdOnlyIso[i] = LHBasedIdOnlyIso[i];
    myLHBasedIdOnlyConv[i] = LHBasedIdOnlyConv[i];
  }
}

void RedEleIDTree::fillLHBasedPFIsoIDBits(int LHBasedPFIsoId[5], int LHBasedPFIsoIdOnlyID[5], int LHBasedPFIsoIdOnlyIso[5], int LHBasedPFIsoIdOnlyConv[5]) {
  for(int i=0; i<5; i++) {
    myLHBasedPFIsoId[i] = LHBasedPFIsoId[i];
    myLHBasedPFIsoIdOnlyID[i] = LHBasedPFIsoIdOnlyID[i];
    myLHBasedPFIsoIdOnlyIso[i] = LHBasedPFIsoIdOnlyIso[i];
    myLHBasedPFIsoIdOnlyConv[i] = LHBasedPFIsoIdOnlyConv[i];
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

void RedEleIDTree::fillFakeRateDenomBits(int isDenom, int isDenomSmurf) {
  myDenomFake = isDenom;
  myDenomFakeSmurf = isDenomSmurf;
}

void RedEleIDTree::fillBDTBasedIDBits(int isBDTOnlyId) {
  myBDTIdOnlyId = isBDTOnlyId;
}

void RedEleIDTree::fillRunInfos(int run, int lumi, int event, int npu[3], int mcmatch) {
  myRun = run;
  myLS = lumi;
  myEvent = event;
  for(int i=0;i<3;i++) myNpu[i]=npu[i];
  myMCMatch = mcmatch;
}
