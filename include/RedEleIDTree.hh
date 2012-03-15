#ifndef RedEleIDTree_H
#define RedEleIDTree_H

#include <TFile.h>
#include <TTree.h>

class RedEleIDTree {
public:
  
  RedEleIDTree(const char *filename);
  ~RedEleIDTree();

  //! add the electron attributes (see below)
  void addAttributesSignal();
  void addAttributesBackground();
  //! add the splitting categories (see below)
  void addCategories();
  void addMore();
  void addElectronIdBits();
  void addDenominatorFakeBits();
  void addIsolations();
  void addGamma();
  //! add run,lumi, event number (for data)
  void addRunInfos();

  //! fill the tree with electron id variables
  //! minimal set
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float See, float Spp, float fbrem, int nbrems, float pt, float eta, int charge);
  //! needed for HZZ BDT
  void fillVariables(float eleEoPout, float EoPout, float EoP, float HoE, float Deta, float Dphi, float s9s25, float s1s9, float See, float Spp, float fbrem, 
                     int nbrems, int nHits, float dcot, float dist, float pt, float eta, int charge, float phiwidth, float etawidth,
                     float IoEmIoP, float eledeta, float d0, float ip3d, float ip3ds, int kfhits, float kfchi2, float e1x5e5x5, int ecaldriven, int matchConv);
  //! additional needed for HWW BDT
  void fillVariables2(float detacalo, float dphicalo, float sep, float dz, float gsfchi2, float emaxovere, float etopovere, float ebottomovere, float eleftovere, float erightovere,
                      float e2ndovere, float e2x5rightovere, float e2x5leftovere, float e2x5topevere, float e2x5bottomovere, 
                      float e2x5maxovere, float e1x5overe, float e2x2overe, float e3x3overe, float e5x5overe, float r9,
                      float phi, float scenergy, float scrawenergy, float scesenergy);

  //! fill the tree with isolation variables
  void fillIsolations(float tkIso, float ecalIso, float hcalIso,
                      float combPFiso,
                      float chaPFiso, float neuPFiso, float phoPFiso);

  //! fill the electron ID bits
  void fillCutBasedIDBits(int CutBasedId[6], int CutBasedIdOnlyID[6], int CutBasedIdOnlyIso[6], int CutBasedIdOnlyConv[6]);
  void fillLHBasedIDBits(int LHBasedId[5], int LHBasedIdOnlyID[5], int LHBasedIdOnlyIso[5], int LHBasedIdOnlyConv[5]);
  void fillLHBasedPFIsoIDBits(int LHBasedPFIsoId[5], int LHBasedPFIsoIdOnlyID[5], int LHBasedPFIsoIdOnlyIso[5], int LHBasedPFIsoIdOnlyConv[5]);
  void fillCiCBasedIDBits(int CiCBasedId[9], int CiCBasedIdOnlyID[9], int CiCBasedIdOnlyIso[9], int CiCBasedIdOnlyConv[9]);
  void fillFakeRateDenomBits(int isDenom, int isDenomSmurf);
  void fillBDTBasedIDBits(int isBDTOnlyId);

  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributesSignal(float zmass, int zeeDec);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(float dphi, float invmass, float met, float pth);
  //! fill the splitting categories of the PDFs
  void fillCategories(int iecal, int iptbin, int iclass, int nbr);
  void fillMore(float nVtx, float rho, float bdthww, float bdthzz);
  void fillMore2(float bdthwwnoip, float bdthzznoip, float bdthzzmc, float pfmva, float like);
  void fillGamma(float atg, float aeg, float ahg, int ig);
  //! fill the run,lumi, event number, mc match
  void fillRunInfos(int run, int lumi, int event, int npu[3], int mcmatch);   

  void store();
  void save();

private:

  float myEleEoPout, myEoPout, myEoP,myHoE,myDeta,myDphi,mys9s25,mys1s9,mySee,mySpp,myFbrem, myPhiWidth, myEtaWidth;
  float myIoEoIoP, myEleDeta, myD0, myIP3d, myIP3dSig, myKFChi2, myE1x5E5x5;
  int myNbrems, myKFHits, myEcalDriven, myMissHits, myMatchConv;
  float myDetaCalo, myDphiCalo, mySep, myDZ, myGSFChi2;
  float mySeedEMaxOverE,mySeedETopOverE,mySeedEBottomOverE,mySeedELeftOverE,mySeedERightOverE,mySeedE2ndOverE,mySeedE2x5RightOverE,mySeedE2x5LeftOverE,mySeedE2x5TopOverE,mySeedE2x5BottomOverE;
  float mySeedE2x5MaxOverE,mySeedE1x5OverE,mySeedE2x2OverE,mySeedE3x3OverE,mySeedE5x5OverE,myR9;
  float myDist, myDcot;
  int myCharge;
  float myEta, myPhi, myPt, mySCEnergy, mySCRawEnergy, mySCESEnergy;
  int myNpu[3];
  int myRun, myLS, myEvent, myMCMatch;

  float myZmass;
  int myZDec;

  int myCutBasedId[6], myCutBasedIdOnlyID[6], myCutBasedIdOnlyIso[6], myCutBasedIdOnlyConv[6];
  int myLHBasedId[5], myLHBasedIdOnlyID[5], myLHBasedIdOnlyIso[5], myLHBasedIdOnlyConv[5];
  int myLHBasedPFIsoId[5], myLHBasedPFIsoIdOnlyID[5], myLHBasedPFIsoIdOnlyIso[5], myLHBasedPFIsoIdOnlyConv[5];
  int myCiCBasedId[9], myCiCBasedIdOnlyID[9], myCiCBasedIdOnlyIso[9], myCiCBasedIdOnlyConv[9];
  int myDenomFake, myDenomFakeSmurf;
  int myBDTIdOnlyId;

  float myQCDDeltaphi;
  float myQCDInvmass;
  float myQCDMet;
  float myQCDPtHat;

  int myiecal;
  int myiptbin;
  int myiclass;
  int mynbrem;

  float myTrkIso;
  float myEcalIso;
  float myHcalIso;
  float myPFCandCombinedIsoHWW;
  float myPFCandChargedIso, myPFCandNeutralIso, myPFCandPhotonIso;

  float myNVtx;
  float myRho;
  float myBdtHww, myBdtHwwNoIP, myBdtHzz, myBdtHzzNoIP, myBdtHzzMC, myPFMVA, myLike;

  float myAbsTrackerIsolGammaCand;
  float myAbsEcalIsolGammaCand;
  float myAbsHcalIsolGammaCand;
  int myIsGamma;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
