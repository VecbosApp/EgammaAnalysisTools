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

  //! fill the tree with electron id variables
  void fillVariables(float EoPout, float EoP, float HoE, float DeltaEta, float DeltaPhi, float s9s25, float s1s9, float SigmaIEtaIEta);
  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributesSignal(int charge, float eta, float pt, float zmass);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(int charge, float eta, float pt, float deltaphi, float invmass, float met, float pth);
  //! fill the splitting categories of the PDFs
  //! iclass: 0=non-showering, 1=showering
  //! iecal: 0=EB, 1=EE
  //! iptbin: 0=<15 GeV, 1=>15GeV
  void fillCategories(int iecal, int iptbin, int iclass); 
  void fillMore(float rit, float rip);
  
  void store();
  void save();

private:

  float myEoPout;
  float myEoP   ;
  float myHoE;
  float myDeltaEta;
  float myDeltaPhi;
  float mys9s25;
  float mys1s9 ;
  float mySigmaIEtaIEta;

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

  int myiecal;
  int myiptbin;
  int myiclass;

  float myRelIsolTag;
  float myRelIsolProbe;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
