//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  8 16:40:13 2010 by ROOT version 5.22/00d
// from TTree T1/tree with only selected events
// found on file: results/merged/PhotonJet_merged.root
//////////////////////////////////////////////////////////

#ifndef checkGammaJet_h
#define checkGammaJet_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class checkGammaJet {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         EoPout;
   Float_t         EoP;
   Float_t         HoE;
   Float_t         deltaEta;
   Float_t         deltaEtaCorr;
   Float_t         deltaPhi;
   Float_t         deltaPhiCorr;
   Float_t         sigmaIEtaIEta;
   Float_t         fBrem;
   Float_t         convDcot;
   Float_t         convDist;
   Int_t           missingHits;
   Float_t         qcdEta;
   Float_t         qcdPt;
   Float_t         qcdDeltaphi;
   Float_t         qcdMet;
   Int_t           qcdNBrem;
   Float_t         absTrackerIsolGammaCand;
   Float_t         absEcalIsolGammaCand;
   Float_t         absHcalIsolGammaCand;
   Int_t           isGamma;
   Float_t         weight;

   // List of branches
   TBranch        *b_EoPout;   //!
   TBranch        *b_EoP;   //!
   TBranch        *b_HoE;   //!
   TBranch        *b_deltaEta;   //!
   TBranch        *b_deltaEtaCorr;   //!
   TBranch        *b_deltaPhi;   //!
   TBranch        *b_deltaPhiCorr;   //!
   TBranch        *b_sigmaIEtaIEta;   //!
   TBranch        *b_fBrem;   //!
   TBranch        *b_convDcot;   //!
   TBranch        *b_convDist;   //!
   TBranch        *b_missingHits;   //!
   TBranch        *b_qcdEta;   //!
   TBranch        *b_qcdPt;   //!
   TBranch        *b_qcdDeltaphi;   //!
   TBranch        *b_qcdMet;   //!
   TBranch        *b_qcdNBrem;   //!
   TBranch        *b_absTrackerIsolGammaCand;   //!
   TBranch        *b_absEcalIsolGammaCand;   //!
   TBranch        *b_absHcalIsolGammaCand;   //!
   TBranch        *b_isGamma;   //!
   TBranch        *b_weight;   //!

   checkGammaJet(TTree *tree=0);
   virtual ~checkGammaJet();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef checkGammaJet_cxx
checkGammaJet::checkGammaJet(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results/merged/PhotonJet_merged.root");
      if (!f) {
         f = new TFile("results/merged/PhotonJet_merged.root");
      }
      tree = (TTree*)gDirectory->Get("T1");

   }
   Init(tree);
}

checkGammaJet::~checkGammaJet()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t checkGammaJet::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t checkGammaJet::LoadTree(Long64_t entry)
{
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

void checkGammaJet::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EoPout", &EoPout, &b_EoPout);
   fChain->SetBranchAddress("EoP", &EoP, &b_EoP);
   fChain->SetBranchAddress("HoE", &HoE, &b_HoE);
   fChain->SetBranchAddress("deltaEta", &deltaEta, &b_deltaEta);
   fChain->SetBranchAddress("deltaEtaCorr", &deltaEtaCorr, &b_deltaEtaCorr);
   fChain->SetBranchAddress("deltaPhi", &deltaPhi, &b_deltaPhi);
   fChain->SetBranchAddress("deltaPhiCorr", &deltaPhiCorr, &b_deltaPhiCorr);
   fChain->SetBranchAddress("sigmaIEtaIEta", &sigmaIEtaIEta, &b_sigmaIEtaIEta);
   fChain->SetBranchAddress("fBrem", &fBrem, &b_fBrem);
   fChain->SetBranchAddress("convDcot", &convDcot, &b_convDcot);
   fChain->SetBranchAddress("convDist", &convDist, &b_convDist);
   fChain->SetBranchAddress("missingHits", &missingHits, &b_missingHits);
   fChain->SetBranchAddress("qcdEta", &qcdEta, &b_qcdEta);
   fChain->SetBranchAddress("qcdPt", &qcdPt, &b_qcdPt);
   fChain->SetBranchAddress("qcdDeltaphi", &qcdDeltaphi, &b_qcdDeltaphi);
   fChain->SetBranchAddress("qcdMet", &qcdMet, &b_qcdMet);
   fChain->SetBranchAddress("qcdNBrem", &qcdNBrem, &b_qcdNBrem);
   fChain->SetBranchAddress("absTrackerIsolGammaCand", &absTrackerIsolGammaCand, &b_absTrackerIsolGammaCand);
   fChain->SetBranchAddress("absEcalIsolGammaCand", &absEcalIsolGammaCand, &b_absEcalIsolGammaCand);
   fChain->SetBranchAddress("absHcalIsolGammaCand", &absHcalIsolGammaCand, &b_absHcalIsolGammaCand);
   fChain->SetBranchAddress("isGamma", &isGamma, &b_isGamma);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t checkGammaJet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void checkGammaJet::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t checkGammaJet::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef checkGammaJet_cxx
