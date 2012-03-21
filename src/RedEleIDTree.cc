#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myFile->mkdir("eleIDdir");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("eleEoPout",       &myEleEoPout,       "eleEoPout/F");
  myTree->Branch("EoPout",          &myEoPout,          "EoPout/F");
  myTree->Branch("EoP",             &myEoP,             "EoP/F");
  myTree->Branch("IoEmIoP",         &myIoEoIoP,         "IoEmIoP/F");
  myTree->Branch("HoE",             &myHoE,             "HoE/F");
  myTree->Branch("eledeta",         &myEleDeta,         "eledeta/F");
  myTree->Branch("deta",            &myDeta,            "deta/F");
  myTree->Branch("dphi",            &myDphi,            "dphi/F");
  myTree->Branch("detacalo",        &myDetaCalo,        "detacalo/F");
  myTree->Branch("dphicalo",        &myDphiCalo,        "dphicalo/F");
  myTree->Branch("s9s25",           &mys9s25,           "s9s25/F");
  myTree->Branch("phiwidth",        &myPhiWidth,        "phiwidth/F");
  myTree->Branch("etawidth",        &myEtaWidth,        "etawidth/F");
  myTree->Branch("see",             &mySee,             "see/F");
  myTree->Branch("sep",             &mySep,             "sep/F");
  myTree->Branch("spp",             &mySpp,             "spp/F");
  myTree->Branch("fbrem",           &myFbrem,           "fbrem/F");
  myTree->Branch("nbrem",           &myNbrems,          "nbrems/I");
  myTree->Branch("missHits",        &myMissHits,        "missHits/I");
  myTree->Branch("dist",            &myDist,            "dist/F");
  myTree->Branch("dcot",            &myDcot,            "dcot/F");
  myTree->Branch("d0",              &myD0,              "d0/F");
  myTree->Branch("dz",              &myDZ,              "dz/F");
  myTree->Branch("ip3d",            &myIP3d,            "ip3d/F");
  myTree->Branch("ip3ds",           &myIP3dSig,         "ip3ds/F");
  myTree->Branch("kfhits",          &myKFHits,          "kfhits/I");
  myTree->Branch("kfchi2",          &myKFChi2,          "kfchi2/F");
  myTree->Branch("gsfchi2",         &myGSFChi2,         "gsfchi2/F");
  myTree->Branch("e1x5e5x5",        &myE1x5E5x5,        "e1x5e5x5/F");
  myTree->Branch("SeedEMaxOverE",   &mySeedEMaxOverE,   "SeedEMaxOverE/F"); 
  myTree->Branch("SeedETopOverE",   &mySeedETopOverE,   "SeedETopOverE/F"); 
  myTree->Branch("SeedEBottomOverE",&mySeedEBottomOverE,"SeedEBottomOverE/F"); 
  myTree->Branch("SeedELeftOverE",  &mySeedELeftOverE,  "SeedELeftOverE/F"); 
  myTree->Branch("SeedERightOverE", &mySeedERightOverE, "SeedERightOverE/F"); 
  myTree->Branch("SeedE2ndOverE",   &mySeedE2ndOverE,   "SeedE2ndOverE/F"); 
  myTree->Branch("SeedE2x5RightOverE",  &mySeedE2x5RightOverE, "SeedE2x5RightOverE/F"); 
  myTree->Branch("SeedE2x5LeftOverE",   &mySeedE2x5LeftOverE,  "SeedE2x5LeftOverE/F"); 
  myTree->Branch("SeedE2x5TopOverE",&mySeedE2x5TopOverE,"SeedE2x5TopOverE/F"); 
  myTree->Branch("SeedE2x5BottomOverE", &mySeedE2x5BottomOverE, "SeedE2x5BottomOverE/F"); 
  myTree->Branch("SeedE2x5MaxOverE",&mySeedE2x5MaxOverE,"SeedE2x5MaxOverE/F"); 
  myTree->Branch("SeedE1x5OverE",   &mySeedE1x5OverE,   "SeedE1x5OverE/F"); 
  myTree->Branch("SeedE2x2OverE",   &mySeedE2x2OverE,   "SeedE2x2OverE/F"); 
  myTree->Branch("SeedE3x3OverE",   &mySeedE3x3OverE,   "SeedE3x3OverE/F"); 
  myTree->Branch("SeedE5x5OverE",   &mySeedE5x5OverE,   "SeedE5x5OverE/F"); 
  myTree->Branch("R9",              &myR9,              "R9/F");
  myTree->Branch("matchConv",       &myMatchConv,       "matchConv/I");
  myTree->Branch("ecaldriven",      &myEcalDriven,      "ecaldriven/I");
  myTree->Branch("scenergy",        &mySCEnergy,        "scenergy/F");
  myTree->Branch("scrawenergy",     &mySCRawEnergy,     "scrawenergy/F");
  myTree->Branch("scesenergy",      &mySCESEnergy,      "scesenergy/F");
  myTree->Branch("pt",              &myPt,              "pt/F");
  myTree->Branch("eta",             &myEta,             "eta/F");
  myTree->Branch("phi",             &myPhi,             "phi/F");
  myTree->Branch("charge",          &myCharge,          "charge/I");
  
}

RedEleIDTree::~RedEleIDTree() {
  delete myFile;
}

void RedEleIDTree::addAttributesSignal() {

  myTree->Branch("mass",       &myZmass,       "mass/F");
  myTree->Branch("zdec",       &myZDec,        "zdec/I");
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
  myTree->Branch("bdthww",     myBdtHww,    "bdthww[2]/F");
  myTree->Branch("newbdthww",  myNewBdtHww, "newbdthww[4]/F");
  myTree->Branch("bdthzz",     myBdtHzz,    "bdthzz[4]/F");
  myTree->Branch("lh",       &myLike,   "lh/F");
  myTree->Branch("pfmva",    &myPFMVA,  "pfmva/F");
  myTree->Branch("vertices", &myNVtx, "vertices/F");
  myTree->Branch("rho",      &myRho,  "rho/F");
}

void RedEleIDTree::addTrackMomenta() {
  myTree->Branch("pcomb",    &myPComb,    "pcomb/F");
  myTree->Branch("pmodegsf", &myPModeGsf, "pmodegsf/F");
  myTree->Branch("pmeangsf", &myPMeanGsf, "pmeangsf/F");
  myTree->Branch("pmeankf",  &myPKf,      "pmeankf/F");
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

void RedEleIDTree::fillVariables(float eleEoPout, float EoPout, float EoP, float HoE, float Deta, float Dphi, float s9s25, float s1s9, float See, float Spp, float fbrem, 
                                 int nbrems, int nHits, float dcot, float dist, float pt, float eta, int charge, float phiwidth, float etawidth,
                                 float IoEmIoP, float eledeta, float d0, float ip3d, float ip3ds, int kfhits, float kfchi2, float e1x5e5x5, int ecaldriven, int matchConv) {
  myEleEoPout=eleEoPout;
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
  myNbrems=nbrems;
  myMissHits=nHits;
  myDist=dist;
  myDcot=dcot;
  myPt=pt;
  myEta=eta;
  myCharge=charge;
  myPhiWidth=phiwidth;
  myEtaWidth=etawidth;
  myIoEoIoP=IoEmIoP;
  myEleDeta=eledeta;
  myD0=d0;
  myIP3d=ip3d;
  myIP3dSig=ip3ds;
  myKFHits=kfhits;
  myKFChi2=kfchi2;
  myE1x5E5x5=e1x5e5x5;
  myEcalDriven=ecaldriven;
  myMatchConv=matchConv;
}

void RedEleIDTree::fillVariables2(float detacalo, float dphicalo, float sep, float dz, float gsfchi2, float emaxovere, float etopovere, float ebottomovere, float eleftovere, float erightovere,
                                  float e2ndovere, float e2x5rightovere, float e2x5leftovere, float e2x5topovere, float e2x5bottomovere, 
                                  float e2x5maxovere, float e1x5overe, float e2x2overe, float e3x3overe, float e5x5overe, float r9,
                                  float phi, float scenergy, float scrawenergy, float scesenergy) {
  myDetaCalo=detacalo;
  myDphiCalo=dphicalo;
  mySep=sep;
  myDZ=dz;
  myGSFChi2=gsfchi2;
  mySeedEMaxOverE=emaxovere;
  mySeedETopOverE=etopovere;
  mySeedEBottomOverE=ebottomovere;
  mySeedELeftOverE=eleftovere;
  mySeedERightOverE=erightovere;
  mySeedE2ndOverE=e2ndovere;
  mySeedE2x5RightOverE=e2x5rightovere;
  mySeedE2x5LeftOverE=e2x5leftovere;
  mySeedE2x5TopOverE=e2x5topovere;
  mySeedE2x5BottomOverE=e2x5bottomovere;
  mySeedE2x5MaxOverE=e2x5maxovere;
  mySeedE1x5OverE=e1x5overe;
  mySeedE2x2OverE=e2x2overe;
  mySeedE3x3OverE=e3x3overe;
  mySeedE5x5OverE=e5x5overe;
  myR9=r9;
  myPhi=phi;
  mySCEnergy=scenergy;
  mySCRawEnergy=scrawenergy;
  mySCESEnergy=scesenergy;
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

void RedEleIDTree::fillAttributesSignal(float zmass, int zdec) {

  myZmass=zmass;
  myZDec=zdec;
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

void RedEleIDTree::fillMore(float nVtx, float rho, float bdthww[2], float newbdthww[4], float bdthzz[4], float pfmva, float like) {
  myNVtx=nVtx;
  myRho=rho;
  for(int i=0; i<2; i++) myBdtHww[i]=bdthww[i];
  for(int i=0; i<4; i++) myNewBdtHww[i]=newbdthww[i];
  for(int i=0; i<4; i++) myBdtHzz[i]=bdthzz[i];
  myPFMVA=pfmva;
  myLike=like;
}

void RedEleIDTree::fillTrackMomenta(float pcomb, float pmodegsf, float pmeangsf, float pkf) {
  myPComb=pcomb;
  myPModeGsf=pmodegsf;
  myPMeanGsf=pmeangsf;
  myPKf=pkf;
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
