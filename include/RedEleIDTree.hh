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
  void addIsolations();
  void addGamma();
  //! add run,lumi, event number (for data)
  void addRunInfos();

  //! fill the tree with electron id variables
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float See, float Spp, float fbrem, int nbrems, float pt, float eta, int charge);
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float s1s9, float See, float Spp, float fBrem, int nHits, float dcot, float dist, float pt, float eta, int charge);

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
  void fillAttributesSignal(float zmass);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(float dphi, float invmass, float met, float pth);
  //! fill the splitting categories of the PDFs
  void fillCategories(int iecal, int iptbin, int iclass, int nbr);
  void fillMore(float nVtx, float rho, float bdthww, float bdthzz);
  void fillGamma(float atg, float aeg, float ahg, int ig);
  //! fill the run,lumi, event number, mc match
  void fillRunInfos(int run, int lumi, int event, int npu[3], int mcmatch);   

  void store();
  void save();

private:

  float myEoPout;
  float myEoP   ;
  float myHoE;
  float myDeta;
  float myDphi;
  float mys9s25;
  float mys1s9 ;
  float mySee;
  float mySpp;
  float myFbrem;
  int myNbrems;
  int myMissHits;
  float myDist, myDcot;
  int myCharge;
  float myEta;
  float myPt;
  int myNpu[3];
  int myRun, myLS, myEvent, myMCMatch;

  float myZmass;

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
  float myBdtHww, myBdtHzz;

  float myAbsTrackerIsolGammaCand;
  float myAbsEcalIsolGammaCand;
  float myAbsHcalIsolGammaCand;
  int myIsGamma;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
