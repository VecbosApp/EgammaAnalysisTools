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

  //! fill the tree with electron id variables
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float See, float Spp, float fbrem, int nbrems, float pt, float eta, int charge);
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float s1s9, float See, float Spp, float fBrem, int nHits, float dcot, float dist, float pt, float eta, int charge);

  //! fill the tree with isolation variables
  void fillIsolations(float tkIso, float ecalIso, float hcalIso);

  //! fill the electron ID bits
  void fillCutBasedIDBits(int CutBasedId[4], int CutBasedIdOnlyID[4], int CutBasedIdOnlyIso[4], int CutBasedIdOnlyConv[4]);
  void fillLHBasedIDBits(int LHBasedId[5], int LHBasedIdOnlyID[5], int LHBasedIdOnlyIso[5], int LHBasedIdOnlyConv[5]);
  void fillCiCBasedIDBits(int CiCBasedId[9], int CiCBasedIdOnlyID[9], int CiCBasedIdOnlyIso[9], int CiCBasedIdOnlyConv[9]);

  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributesSignal(float zmass);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(float dphi, float invmass, float met, float pth);
  //! fill the splitting categories of the PDFs
  void fillCategories(int iecal, int iptbin, int iclass, int nbr);
  void fillMore(int nVtx, float rho);
  void fillGamma(float atg, float aeg, float ahg, int ig);

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

  float myZmass;

  int myCutBasedId[4], myCutBasedIdOnlyID[4], myCutBasedIdOnlyIso[4], myCutBasedIdOnlyConv[4];
  int myLHBasedId[5], myLHBasedIdOnlyID[5], myLHBasedIdOnlyIso[5], myLHBasedIdOnlyConv[5];
  int myCiCBasedId[9], myCiCBasedIdOnlyID[9], myCiCBasedIdOnlyIso[9], myCiCBasedIdOnlyConv[9];

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

  int myNVtx;
  float myRho;

  float myAbsTrackerIsolGammaCand;
  float myAbsEcalIsolGammaCand;
  float myAbsHcalIsolGammaCand;
  int myIsGamma;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
