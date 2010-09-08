//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 15 01:31:13 2010 by ROOT version 5.22/00a
// from TTree dataset/dataset with sWeights
// found on file: ../results/sPlotsTree/sPlots_tree.root
//////////////////////////////////////////////////////////

#ifndef sPlotsPdfsComparison_h
#define sPlotsPdfsComparison_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <iostream>
#include <vector>

class sPlotsPdfsComparison {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types for data
   Double_t        N_sig_sw;
   Double_t        L_N_sig;
   Double_t        N_qcd_sw;
   Double_t        L_N_qcd;
   Double_t        trackerIso;
   Double_t        ecalJIso;
   Double_t        ecalGTIso;
   Double_t        hcalIso;
   Double_t        combinedIso;
   Double_t        classification;
   Double_t        deta;
   Double_t        dphi;
   Double_t        hoe;
   Double_t        see;
   Double_t        eop;
   Double_t        fbrem;
   Double_t        met;
   Double_t        tcmet;
   Double_t        pfmet;
   Double_t        mt;
   Double_t        tcmt;
   Double_t        pfmt;
   Double_t        pt;
   Double_t        eta;
   Double_t        phi;
   Double_t        charge;
   Double_t        weight;
   Double_t        nPFJets;
   Double_t        nJets;
   Double_t        event;

   Float_t       f_trackerIso;
   Float_t       f_ecalJIso;
   Float_t       f_ecalGTIso;
   Float_t       f_hcalIso;
   Float_t       f_combinedIso;
   Float_t       f_classification;
   Float_t       f_deta;
   Float_t       f_dphi;
   Float_t       f_hoe;
   Float_t       f_see;
   Float_t       f_eop;
   Float_t       f_fbrem;
   Float_t       f_met;
   Float_t       f_tcmet;
   Float_t       f_pfmet;
   Float_t       f_mt;
   Float_t       f_tcmt;
   Float_t       f_pfmt;
   Float_t       f_pt;
   Float_t       f_eta;
   Float_t       f_phi;
   Int_t         f_charge;
   Float_t       f_weight;
   Int_t         f_nPFJets;
   Int_t         f_nJets;
   Float_t       f_event;

   // List of branches
   TBranch        *b_N_sig_sw;   //!
   TBranch        *b_L_N_sig;   //!
   TBranch        *b_N_qcd_sw;   //!
   TBranch        *b_L_N_qcd;   //!
   TBranch        *b_trackerIso;   //!
   TBranch        *b_ecalJIso;   //!
   TBranch        *b_ecalGTIso;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_combinedIso;   //!
   TBranch        *b_classification;   //!
   TBranch        *b_deta;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_hoe;   //!
   TBranch        *b_see;   //!
   TBranch        *b_eop;   //!
   TBranch        *b_fbrem;   //!
   TBranch        *b_met;   //!
   TBranch        *b_tcmet;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_mt;   //!
   TBranch        *b_tcmt;   //!
   TBranch        *b_pfmt;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_nPFJets;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_event;   //!

   sPlotsPdfsComparison();
   virtual ~sPlotsPdfsComparison();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *treel, int isMC=1);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     bookHistosVariableBinning();
   virtual void     bookHistosFixedBinning();
   virtual void     InitCuts();
   virtual void     doSignalsPlots(bool what) { m_doSignal = what; }

 protected:

   bool m_isMC;
   bool m_doSignal;

   TH1F *etaClassEle;   
   TH1F *dPhiClassEle[2];
   TH1F *dEtaClassEle[2];
   TH1F *EoPClassEle[2];
   TH1F *HoEClassEle[2];
   TH1F *sigmaIEtaIEtaClassEle[2];
   TH1F *fbremClassEle[2];
   TH1F *phiClassEle[2];
   TH1F *chargeClassEle[2];

   // cuts
   std::vector<float> WP70_EB_sup, WP70_EB_inf, WP70_EE_sup, WP70_EE_inf;

};

#endif

#ifdef sPlotsPdfsComparison_cxx
sPlotsPdfsComparison::sPlotsPdfsComparison() {}

sPlotsPdfsComparison::~sPlotsPdfsComparison()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sPlotsPdfsComparison::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sPlotsPdfsComparison::LoadTree(Long64_t entry)
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

void sPlotsPdfsComparison::Init(TTree *tree, int isMC)
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

   if(isMC) {
     std::cout << "Setting branches for MC tree" << std::endl;
     m_isMC = 1;
     fChain->SetBranchAddress("trackerIso", &f_trackerIso, &b_trackerIso);
     fChain->SetBranchAddress("ecalJIso", &f_ecalJIso, &b_ecalJIso);
     fChain->SetBranchAddress("ecalGTIso", &f_ecalGTIso, &b_ecalGTIso);
     fChain->SetBranchAddress("hcalIso", &f_hcalIso, &b_hcalIso);
     fChain->SetBranchAddress("combinedIso", &f_combinedIso, &b_combinedIso);
     fChain->SetBranchAddress("classification", &f_classification, &b_classification);
     fChain->SetBranchAddress("deta", &f_deta, &b_deta);
     fChain->SetBranchAddress("dphi", &f_dphi, &b_dphi);
     fChain->SetBranchAddress("hoe", &f_hoe, &b_hoe);
     fChain->SetBranchAddress("see", &f_see, &b_see);
     fChain->SetBranchAddress("eop", &f_eop, &b_eop);
     fChain->SetBranchAddress("fbrem", &f_fbrem, &b_fbrem);
     fChain->SetBranchAddress("met", &f_met, &b_met);
     fChain->SetBranchAddress("tcmet", &f_tcmet, &b_tcmet);
     fChain->SetBranchAddress("pfmet", &f_pfmet, &b_pfmet);
     fChain->SetBranchAddress("mt", &f_mt, &b_mt);
     fChain->SetBranchAddress("tcmt", &f_tcmt, &b_tcmt);
     fChain->SetBranchAddress("pfmt", &f_pfmt, &b_pfmt);
     fChain->SetBranchAddress("pt", &f_pt, &b_pt);
     fChain->SetBranchAddress("eta", &f_eta, &b_eta);
     fChain->SetBranchAddress("phi", &f_phi, &b_phi);
     fChain->SetBranchAddress("charge", &f_charge, &b_charge);
     fChain->SetBranchAddress("weight", &f_weight, &b_weight);
     fChain->SetBranchAddress("event", &f_event, &b_event);
     fChain->SetBranchAddress("nPFJets", &f_nPFJets, &b_nPFJets);
     fChain->SetBranchAddress("nJets", &f_nJets, &b_nJets);
   } else {
     m_isMC = 0;
     fChain->SetBranchAddress("N_sig_sw", &N_sig_sw, &b_N_sig_sw);
     fChain->SetBranchAddress("L_N_sig", &L_N_sig, &b_L_N_sig);
     fChain->SetBranchAddress("N_qcd_sw", &N_qcd_sw, &b_N_qcd_sw);
     fChain->SetBranchAddress("L_N_qcd", &L_N_qcd, &b_L_N_qcd);
     fChain->SetBranchAddress("trackerIso", &trackerIso, &b_trackerIso);
     fChain->SetBranchAddress("ecalJIso", &ecalJIso, &b_ecalJIso);
     fChain->SetBranchAddress("ecalGTIso", &ecalGTIso, &b_ecalGTIso);
     fChain->SetBranchAddress("hcalIso", &hcalIso, &b_hcalIso);
     fChain->SetBranchAddress("combinedIso", &combinedIso, &b_combinedIso);
     fChain->SetBranchAddress("classification", &classification, &b_classification);
     fChain->SetBranchAddress("deta", &deta, &b_deta);
     fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
     fChain->SetBranchAddress("hoe", &hoe, &b_hoe);
     fChain->SetBranchAddress("see", &see, &b_see);
     fChain->SetBranchAddress("eop", &eop, &b_eop);
     fChain->SetBranchAddress("fbrem", &fbrem, &b_fbrem);
     fChain->SetBranchAddress("met", &met, &b_met);
     fChain->SetBranchAddress("tcmet", &tcmet, &b_tcmet);
     fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
     fChain->SetBranchAddress("mt", &mt, &b_mt);
     fChain->SetBranchAddress("tcmt", &tcmt, &b_tcmt);
     fChain->SetBranchAddress("pfmt", &pfmt, &b_pfmt);
     fChain->SetBranchAddress("pt", &pt, &b_pt);
     fChain->SetBranchAddress("eta", &eta, &b_eta);
     fChain->SetBranchAddress("phi", &phi, &b_phi);
     fChain->SetBranchAddress("charge", &charge, &b_charge);
     fChain->SetBranchAddress("weight", &weight, &b_weight);
     fChain->SetBranchAddress("event", &event, &b_event);
     fChain->SetBranchAddress("nPFJets", &nPFJets, &b_nPFJets);
     fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   }
   Notify();
}

Bool_t sPlotsPdfsComparison::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sPlotsPdfsComparison::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sPlotsPdfsComparison::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sPlotsPdfsComparison_cxx
