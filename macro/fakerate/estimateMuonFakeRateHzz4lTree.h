//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May  6 15:10:06 2012 by ROOT version 5.33/03
// from TTree probe_tree/probe_tree
// found on file: hzzTree.root
//////////////////////////////////////////////////////////

#ifndef estimateMuonFakeRateHzz4lTree_h
#define estimateMuonFakeRateHzz4lTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "../../../Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class estimateMuonFakeRateHzz4lTree {
public :

  enum wphzz {
    khzzPfIso = 0,
    khzzMvaPfIso
  };


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         bdtIdDz;
   Float_t         bdtIsoDz;
   Float_t         bdtIsoPnp;
   Float_t         ecalIso;
   Float_t         eta;
   Float_t         hcalIso;
   Float_t         pfIsoChHad04;
   Float_t         pfIsoNHad04_NoEA;
   Float_t         pfIsoPhoton04_NoEA;
   Float_t         phi;
   Float_t         pt;
   Float_t         rho;
   Float_t         sip3d;
   Float_t         tkIso;
   Float_t         vertices;
   Int_t           id;
   Int_t           id2012;
   Int_t           l1trig;
   Int_t           l2trig;
   Int_t           l3trig;
   Int_t           pfid;
   Int_t           globalmu;
   Int_t           trackermu;

   // List of branches
   TBranch        *b_bdtIdDz;   //!
   TBranch        *b_bdtIsoDz;   //!
   TBranch        *b_bdtIsoPnp;   //!
   TBranch        *b_ecalIso;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_pfIsoChHad04;   //!
   TBranch        *b_pfIsoNHad04_NoEA;   //!
   TBranch        *b_pfIsoPhoton04_NoEA;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_sip3d;   //!
   TBranch        *b_tkIso;   //!
   TBranch        *b_vertices;   //!
   TBranch        *b_id;   //!
   TBranch        *b_id2012;   //!
   TBranch        *b_l1trig;   //!
   TBranch        *b_l2trig;   //!
   TBranch        *b_l3trig;   //!
   TBranch        *b_pfid;   //!
   TBranch        *b_globalmu;   //!
   TBranch        *b_trackermu;   //!

   estimateMuonFakeRateHzz4lTree(TTree *tree=0);
   virtual ~estimateMuonFakeRateHzz4lTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char *filename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t   passRefMuSel();
};

#endif

#ifdef estimateMuonFakeRateHzz4lTree_cxx
estimateMuonFakeRateHzz4lTree::estimateMuonFakeRateHzz4lTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("hzzTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("hzzTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("hzzTree.root:/zllmtree");
      dir->GetObject("probe_tree",tree);

   }
   Init(tree);
}

estimateMuonFakeRateHzz4lTree::~estimateMuonFakeRateHzz4lTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t estimateMuonFakeRateHzz4lTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t estimateMuonFakeRateHzz4lTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void estimateMuonFakeRateHzz4lTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("bdtIdDz", &bdtIdDz, &b_bdtIdDz);
   fChain->SetBranchAddress("bdtIsoDz", &bdtIsoDz, &b_bdtIsoDz);
   fChain->SetBranchAddress("bdtIsoPnp", &bdtIsoPnp, &b_bdtIsoPnp);
   fChain->SetBranchAddress("ecalIso", &ecalIso, &b_ecalIso);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("hcalIso", &hcalIso, &b_hcalIso);
   fChain->SetBranchAddress("pfIsoChHad04", &pfIsoChHad04, &b_pfIsoChHad04);
   fChain->SetBranchAddress("pfIsoNHad04_NoEA", &pfIsoNHad04_NoEA, &b_pfIsoNHad04_NoEA);
   fChain->SetBranchAddress("pfIsoPhoton04_NoEA", &pfIsoPhoton04_NoEA, &b_pfIsoPhoton04_NoEA);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("sip3d", &sip3d, &b_sip3d);
   fChain->SetBranchAddress("tkIso", &tkIso, &b_tkIso);
   fChain->SetBranchAddress("vertices", &vertices, &b_vertices);
   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("id2012", &id2012, &b_id2012);
   fChain->SetBranchAddress("l1trig", &l1trig, &b_l1trig);
   fChain->SetBranchAddress("l2trig", &l2trig, &b_l2trig);
   fChain->SetBranchAddress("l3trig", &l3trig, &b_l3trig);
   fChain->SetBranchAddress("pfid", &pfid, &b_pfid);
   fChain->SetBranchAddress("globalmu", &globalmu, &b_globalmu);
   fChain->SetBranchAddress("trackermu", &trackermu, &b_trackermu);
   Notify();
}

Bool_t estimateMuonFakeRateHzz4lTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void estimateMuonFakeRateHzz4lTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t estimateMuonFakeRateHzz4lTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Bool_t estimateMuonFakeRateHzz4lTree::passRefMuSel() {
  
  float abseta=fabs(eta);

  MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = MuonEffectiveArea::kMuEAFall11MC;
  MuonEffectiveArea::MuonEffectiveAreaType effAreaGamma_   = MuonEffectiveArea::kMuGammaIso04;
  MuonEffectiveArea::MuonEffectiveAreaType effAreaNeutralHad_ = MuonEffectiveArea::kMuNeutralHadronIso04;

  float eff_area_ga  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);

  float iso = pfIsoChHad04;
  iso += max<float>(0.,pfIsoNHad04_NoEA-eff_area_nh*rho + pfIsoPhoton04_NoEA-eff_area_ga*rho);

  return (pfid && iso<0.4 && fabs(sip3d)<4.0);
} 
#endif // #ifdef estimateMuonFakeRateHzz4lTree_cxx
