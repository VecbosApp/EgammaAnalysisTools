#ifndef RedEleIDTree_H
#define RedEleIDTree_H

#include <TFile.h>
#include <TTree.h>

class RedEleIDTree {
public:
  
  RedEleIDTree(const char *filename);
  ~RedEleIDTree();

  //! add the electron attributes (see below)
  void addAttributes();
  //! add the splitting categories (see below)
  void addCategories();

  //! fill the tree with electron id variables
  void fillVariables(float EoPout, float HoE, float DeltaEta, float DeltaPhi, float s9s25, float SigmaEtaEta);
  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributes(int charge, float eta, float pt, float zmass);
  //! fill the splitting categories of the PDFs
  //! iclass: 0=non-showering, 1=showering
  //! iecal: 0=EB, 1=EE
  //! iptbin: 0=<15 GeV, 1=>15GeV
  void fillCategories(int iecal, int iptbin, int iclass); 

  void store();
  void save();

private:

  float myEoPout;
  float myHoE;
  float myDeltaEta;
  float myDeltaPhi;
  float mys9s25;
  float mySigmaEtaEta;

  int myCharge;
  float myEta;
  float myPt;
  float myZmass;

  int myiecal;
  int myiptbin;
  int myiclass;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
