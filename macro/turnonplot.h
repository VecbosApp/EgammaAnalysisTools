//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  7 18:20:23 2013 by ROOT version 5.30/02
// from TTree T1/tree for turn on curves
// found on file: all.root
//////////////////////////////////////////////////////////

#ifndef turnonplot_h
#define turnonplot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class turnonplot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         ptGamma;
   Float_t         etaGamma;
   Float_t         phiGamma;
   Int_t           idGamma;
   Int_t           hltmatchGamma;
   Float_t         mass;
   Float_t         ptEle1;
   Float_t         ptEle2;
   Float_t         etaEle1;
   Float_t         etaEle2;
   Int_t           tagEle1;
   Int_t           tagEle2;
   Int_t           isAZ;
   Int_t           isAZ_tight;
   Int_t           passHLT20;
   Int_t           passHLT30;
   Int_t           passHLT50;
   Int_t           passHLT75;
   Int_t           passHLT90;
   Float_t         vertices;
   Float_t         rho;
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;

   // List of branches
   TBranch        *b_ptGamma;   //!
   TBranch        *b_etaGamma;   //!
   TBranch        *b_phiGamma;   //!
   TBranch        *b_idGamma;   //!
   TBranch        *b_hltmatchGamma;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_ptEle1;   //!
   TBranch        *b_ptEle2;   //!
   TBranch        *b_etaEle1;   //!
   TBranch        *b_etaEle2;   //!
   TBranch        *b_tagEle1;   //!
   TBranch        *b_tagEle2;   //!
   TBranch        *b_isAZ;   //!
   TBranch        *b_isAZ_tight;   //!
   TBranch        *b_passHLT20;   //!
   TBranch        *b_passHLT30;   //!
   TBranch        *b_passHLT50;   //!
   TBranch        *b_passHLT75;   //!
   TBranch        *b_passHLT90;   //!
   TBranch        *b_vertices;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!

   turnonplot(TTree *tree=0);
   virtual ~turnonplot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef turnonplot_cxx
turnonplot::turnonplot(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("allDoubleElectron_base75.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("allDoubleElectron_base75.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("allDoubleElectron_base75.root:/turnoncurvedir");
      dir->GetObject("T1",tree);

   }
   Init(tree);
}

turnonplot::~turnonplot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t turnonplot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t turnonplot::LoadTree(Long64_t entry)
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

void turnonplot::Init(TTree *tree)
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

   fChain->SetBranchAddress("ptGamma", &ptGamma, &b_ptGamma);
   fChain->SetBranchAddress("etaGamma", &etaGamma, &b_etaGamma);
   fChain->SetBranchAddress("phiGamma", &phiGamma, &b_phiGamma);
   fChain->SetBranchAddress("idGamma", &idGamma, &b_idGamma);
   fChain->SetBranchAddress("hltmatchGamma", &hltmatchGamma, &b_hltmatchGamma);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("ptEle1", &ptEle1, &b_ptEle1);
   fChain->SetBranchAddress("ptEle2", &ptEle2, &b_ptEle2);
   fChain->SetBranchAddress("etaEle1", &etaEle1, &b_etaEle1);
   fChain->SetBranchAddress("etaEle2", &etaEle2, &b_etaEle2);
   fChain->SetBranchAddress("tagEle1", &tagEle1, &b_tagEle1);
   fChain->SetBranchAddress("tagEle2", &tagEle2, &b_tagEle2);
   fChain->SetBranchAddress("isAZ", &isAZ, &b_isAZ);
   fChain->SetBranchAddress("isAZ_tight", &isAZ_tight, &b_isAZ_tight);
   fChain->SetBranchAddress("passHLT20", &passHLT20, &b_passHLT20);
   fChain->SetBranchAddress("passHLT30", &passHLT30, &b_passHLT30);
   fChain->SetBranchAddress("passHLT50", &passHLT50, &b_passHLT50);
   fChain->SetBranchAddress("passHLT75", &passHLT75, &b_passHLT75);
   fChain->SetBranchAddress("passHLT90", &passHLT90, &b_passHLT90);
   fChain->SetBranchAddress("vertices", &vertices, &b_vertices);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   Notify();
}

Bool_t turnonplot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void turnonplot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t turnonplot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef turnonplot_cxx
