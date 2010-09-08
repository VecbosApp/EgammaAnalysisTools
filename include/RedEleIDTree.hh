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
  void fillVariables(float EoPout, float EoP, float HoE, float DeltaEta, float DeltaPhi, float s9s25, float s1s9, float SigmaIEtaIEta);
  void fillVariables(float EoPout, float EoP, float HoE, float DeltaEta, float DeltaEtaCorr, float DeltaPhi, float DeltaPhiCorr, float s9s25, float s1s9, float SigmaIEtaIEta, float fBrem, int nHits, float dcot, float dist);

  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributesSignal(int charge, float eta, float pt, float zmass);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(int charge, float eta, float pt, float deltaphi, float invmass, float met, float pth);
  void fillAttributesBackground(int charge, float eta, float pt, float deltaphi, float invmass, float met, float pth, int nBrem);
  //! fill the splitting categories of the PDFs
  //! iclass: 0=non-showering, 1=showering
  //! iecal: 0=EB, 1=EE
  //! iptbin: 0=<15 GeV, 1=>15GeV
  void fillCategories(int iecal, int iptbin, int iclass); 
  void fillMore(float rit, float rip);
  void fillGamma(float atg, float aeg, float ahg, int ig);

  void store();
  void save();

private:

  float myEoPout;
  float myEoP   ;
  float myHoE;
  float myDeltaEta;
  float myDeltaPhi;
  float myDeltaEtaCorr;
  float myDeltaPhiCorr;
  float mys9s25;
  float mys1s9 ;
  float mySigmaIEtaIEta;
  float myFBrem;
  int myMissingHits;
  float myConvDcot, myConvDist;

  int myCharge;
  float myEta;
  float myPt;
  float myZmass;

  int myQCDCharge;
  float myQCDEta;
  float myQCDPt;
  float myQCDDeltaphi;
  float myQCDInvmass;
  float myQCDMet;
  float myQCDPtHat;
  int   myQCDNBrem;

  int myiecal;
  int myiptbin;
  int myiclass;

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
