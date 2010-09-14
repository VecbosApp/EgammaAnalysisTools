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
  void addGamma();

  //! fill the tree with electron id variables
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float s1s9, float See, float Spp, float pt, float eta, int charge);
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DEtaUncorr, float DPhi, float DPhiUncorr, float s9s25, float s1s9, float See, float Spp, float fBrem, int nHits, float dcot, float dist, float pt, float eta, int charge);

  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributesSignal(float zmass);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(float dphi, float invmass, float met, float pth);
  //! fill the splitting categories of the PDFs
  void fillCategories(int iecal, int iptbin, int iclass, int nbr);
  void fillMore(float rit, float rip);
  void fillGamma(float atg, float aeg, float ahg, int ig);

  void store();
  void save();

private:

  float myEoPout;
  float myEoP   ;
  float myHoE;
  float myDeta;
  float myDphi;
  float myDetaUncorr;
  float myDphiUncorr;
  float mys9s25;
  float mys1s9 ;
  float mySee;
  float mySpp;
  float myFbrem;
  int myMissHits;
  float myDist, myDcot;
  int myCharge;
  float myEta;
  float myPt;

  float myZmass;

  float myQCDDeltaphi;
  float myQCDInvmass;
  float myQCDMet;
  float myQCDPtHat;

  int myiecal;
  int myiptbin;
  int myiclass;
  int mynbrem;

  float myRelIsolTag;
  float myRelIsolProbe;

  float myAbsTrackerIsolGammaCand;
  float myAbsEcalIsolGammaCand;
  float myAbsHcalIsolGammaCand;
  int myIsGamma;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
