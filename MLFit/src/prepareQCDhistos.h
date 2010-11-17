//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 15 11:58:18 2010 by ROOT version 5.22/00d
// from TTree T1/eleID tree
// found on file: results_data/dataset_jetmettau_mergedTree.root
//////////////////////////////////////////////////////////

#ifndef prepareQCDhistos_h
#define prepareQCDhistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include "EgammaAnalysisTools/include/ElectronLikelihood.h"

class prepareQCDhistos {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Float_t         EoPout;
  Float_t         EoP;
  Float_t         HoE;
  Float_t         deta;
  Float_t         dphi;
  Float_t         s9s25;
  Float_t         s1s9;
  Float_t         see;
  Float_t         spp;
  Float_t         fbrem;
  Int_t           missHits;
  Float_t         dist;
  Float_t         dcot;
  Float_t         pt;
  Float_t         eta;
  Int_t           charge;
  Float_t         qcdDeltaphi;
  Float_t         qcdInvmass;
  Float_t         qcdMet;
  Float_t         qcdPtHat;
  Int_t           iecal;
  Int_t           iptbin;
  Int_t           iclass;
  Int_t           nbrem;
  Float_t         absTrackerIsolGammaCand;
  Float_t         absEcalIsolGammaCand;
  Float_t         absHcalIsolGammaCand;
  Int_t           isGamma;
  
  // List of branches
  TBranch        *b_EoPout;   //!
  TBranch        *b_EoP;   //!
  TBranch        *b_HoE;   //!
  TBranch        *b_deta;   //!
  TBranch        *b_dphi;   //!
  TBranch        *b_s9s25;   //!
  TBranch        *b_s1s9;   //!
  TBranch        *b_see;   //!
  TBranch        *b_spp;   //!
  TBranch        *b_fbrem;   //!
  TBranch        *b_missHits;   //!
  TBranch        *b_dist;   //!
  TBranch        *b_dcot;   //!
  TBranch        *b_pt;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_charge;   //!
  TBranch        *b_qcdDeltaphi;   //!
  TBranch        *b_qcdInvmass;   //!
  TBranch        *b_qcdMet;   //!
  TBranch        *b_qcdPtHat;   //!
  TBranch        *b_iecal;   //!
  TBranch        *b_iptbin;   //!
  TBranch        *b_iclass;   //!
  TBranch        *b_nbrem;   //!
  TBranch        *b_absTrackerIsolGammaCand;   //!
  TBranch        *b_absEcalIsolGammaCand;   //!
  TBranch        *b_absHcalIsolGammaCand;   //!
  TBranch        *b_isGamma;   //!
  
  prepareQCDhistos(TTree *tree=0);
  virtual ~prepareQCDhistos();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  void bookFullHistos();
  virtual float likelihoodRatio(ElectronLikelihood &lh);
  
  // histos [ecalsubdet][ptbin]
  TH1F *dPhiUnsplitEle[2][2];
  TH1F *dEtaUnsplitEle[2][2];
  TH1F *EoPUnsplitEle[2][2];
  TH1F *HoEUnsplitEle[2][2];  
  TH1F *sigmaIEtaIEtaUnsplitEle[2][2];
  TH1F *sigmaIPhiIPhiUnsplitEle[2][2];
  TH1F *fBremUnsplitEle[2][2];
  TH1F *lhUnsplitEle[2][2];
  
  // histos [ecalsubdet][ptbin][class]
  TH1F *dPhiClassEle[2][2][2];
  TH1F *dEtaClassEle[2][2][2];
  TH1F *EoPClassEle[2][2][2];
  TH1F *HoEClassEle[2][2][2];
  TH1F *sigmaIEtaIEtaClassEle[2][2][2];
  TH1F *sigmaIPhiIPhiClassEle[2][2][2];
  TH1F *fBremClassEle[2][2][2];
  TH1F *lhClassEle[2][2][2];

  // the likelihood algorithm
  ElectronLikelihood *LH;
};

#endif

#ifdef prepareQCDhistos_cxx
prepareQCDhistos::prepareQCDhistos(TTree *tree) {
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results/trees/QCD_Pt-20_TuneD6T_tree.root");
    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results_data/mergedTree.root");
    if (!f) {
      f = new TFile("results/trees/QCD_Pt-20_TuneD6T_tree.root");
      // f = new TFile("results_data/mergedTree.root");
    }
    tree = (TTree*)gDirectory->Get("T1");  
  }
  Init(tree);
}

prepareQCDhistos::~prepareQCDhistos() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t prepareQCDhistos::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t prepareQCDhistos::LoadTree(Long64_t entry) {

  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void prepareQCDhistos::Init(TTree *tree) {
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("EoPout", &EoPout, &b_EoPout);
  fChain->SetBranchAddress("EoP", &EoP, &b_EoP);
  fChain->SetBranchAddress("HoE", &HoE, &b_HoE);
  fChain->SetBranchAddress("deta", &deta, &b_deta);
  fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
  fChain->SetBranchAddress("s9s25", &s9s25, &b_s9s25);
  fChain->SetBranchAddress("s1s9", &s1s9, &b_s1s9);
  fChain->SetBranchAddress("see", &see, &b_see);
  fChain->SetBranchAddress("spp", &spp, &b_spp);
  fChain->SetBranchAddress("fbrem", &fbrem, &b_fbrem);
  fChain->SetBranchAddress("missHits", &missHits, &b_missHits);
  fChain->SetBranchAddress("dist", &dist, &b_dist);
  fChain->SetBranchAddress("dcot", &dcot, &b_dcot);
  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("charge", &charge, &b_charge);
  fChain->SetBranchAddress("qcdDeltaphi", &qcdDeltaphi, &b_qcdDeltaphi);
  fChain->SetBranchAddress("qcdInvmass", &qcdInvmass, &b_qcdInvmass);
  fChain->SetBranchAddress("qcdMet", &qcdMet, &b_qcdMet);
  fChain->SetBranchAddress("qcdPtHat", &qcdPtHat, &b_qcdPtHat);
  fChain->SetBranchAddress("iecal", &iecal, &b_iecal);
  fChain->SetBranchAddress("iptbin", &iptbin, &b_iptbin);
  fChain->SetBranchAddress("iclass", &iclass, &b_iclass);
  fChain->SetBranchAddress("nbrem", &nbrem, &b_nbrem);
  fChain->SetBranchAddress("absTrackerIsolGammaCand", &absTrackerIsolGammaCand, &b_absTrackerIsolGammaCand);
  fChain->SetBranchAddress("absEcalIsolGammaCand", &absEcalIsolGammaCand, &b_absEcalIsolGammaCand);
  fChain->SetBranchAddress("absHcalIsolGammaCand", &absHcalIsolGammaCand, &b_absHcalIsolGammaCand);
  fChain->SetBranchAddress("isGamma", &isGamma, &b_isGamma);
  Notify();
}

Bool_t prepareQCDhistos::Notify() {
  return kTRUE;
}

void prepareQCDhistos::Show(Long64_t entry) {
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t prepareQCDhistos::Cut(Long64_t entry){
  return 1;
}
#endif // #ifdef prepareQCDhistos_cxx
