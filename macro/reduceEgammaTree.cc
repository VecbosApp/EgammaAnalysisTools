#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <PUWeight.C>

using namespace std;

void reduce(const char* filename) {

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("eleIDdir/T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  if ( treeOrig ) {
  int nentriesOrig = treeOrig->GetEntries();

  TFile *fileNew = TFile::Open(filename,"recreate");
  fileNew->mkdir("eleIDdir");
  TTree *treeNew = new TTree("T1","tree with only selected events");

  // Declaration of leaf types
   Float_t         puW;
   Float_t         mass;
   Float_t         pt;
   Float_t         eta;
   Float_t         vertices;
   Int_t           LHPFIsoBasedIdOlyID[5];
   Int_t           LHPFIsoBasedIdOnlyIso[5];
   Int_t           LHPFIsoBasedIdOnlyConv[5];
   Int_t           CutBasedIdOlyID[6];
   Int_t           WP95;
   Int_t           WP90;
   Int_t           WP85;
   Int_t           WP80;
   Int_t           WP70;
   Int_t           WPSmurf;
   Int_t           BDTIdOnlyId;
   Int_t           BDTWP;
   Int_t           DenomFake;
   Int_t           DenomFakeSmurf;
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Int_t           npu[3];
   Int_t           mcmatch;
   Float_t         bdthww;
   Float_t         bdthzz;
   Float_t         rho;
   Float_t         combPFIsoHWW;
   Float_t         chaPFIso;
   Float_t         neuPFIso;
   Float_t         phoPFIso;
   Float_t         EoPout;
   Float_t         EoP;
   Float_t         HoE;
   Float_t         deta;
   Float_t         dphi;
   Float_t         s9s25;
   Float_t         see;
   Float_t         spp;
   Float_t         fbrem;
   Int_t           nbrem;
   Int_t           missHits;
   Float_t         dist;
   Float_t         dcot;


   treeOrig->SetBranchAddress("mass", &mass);
   treeOrig->SetBranchAddress("pt", &pt);
   treeOrig->SetBranchAddress("eta", &eta);
   treeOrig->SetBranchAddress("vertices", &vertices);
   treeOrig->SetBranchAddress("LHPFIsoBasedIdOlyID", LHPFIsoBasedIdOlyID);
   treeOrig->SetBranchAddress("LHPFIsoBasedIdOnlyIso", LHPFIsoBasedIdOnlyIso);
   treeOrig->SetBranchAddress("LHPFIsoBasedIdOnlyConv", LHPFIsoBasedIdOnlyConv);
   treeOrig->SetBranchAddress("CutBasedIdOlyID", CutBasedIdOlyID);
   treeOrig->SetBranchAddress("WP95", &WP95);
   treeOrig->SetBranchAddress("WP90", &WP90);
   treeOrig->SetBranchAddress("WP85", &WP85);
   treeOrig->SetBranchAddress("WP80", &WP80);
   treeOrig->SetBranchAddress("WP70", &WP70);
   treeOrig->SetBranchAddress("WPSmurf", &WPSmurf);
   treeOrig->SetBranchAddress("BDTIdOnlyId", &BDTIdOnlyId);
   treeOrig->SetBranchAddress("DenomFake", &DenomFake);
   treeOrig->SetBranchAddress("DenomFakeSmurf", &DenomFakeSmurf);
   treeOrig->SetBranchAddress("run", &run);
   treeOrig->SetBranchAddress("lumi", &lumi);
   treeOrig->SetBranchAddress("event", &event);
   treeOrig->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
   treeOrig->SetBranchAddress("chaPFIso", &chaPFIso);
   treeOrig->SetBranchAddress("neuPFIso", &neuPFIso);
   treeOrig->SetBranchAddress("phoPFIso", &phoPFIso);
   treeOrig->SetBranchAddress("npu", npu);
   treeOrig->SetBranchAddress("mcmatch", &mcmatch);
   treeOrig->SetBranchAddress("bdthww", &bdthww);
   treeOrig->SetBranchAddress("bdthzz", &bdthzz);
   treeOrig->SetBranchAddress("rho", &rho);
   treeOrig->SetBranchAddress("EoPout", &EoPout);
   treeOrig->SetBranchAddress("EoP", &EoP);
   treeOrig->SetBranchAddress("HoE", &HoE);
   treeOrig->SetBranchAddress("deta", &deta);
   treeOrig->SetBranchAddress("dphi", &dphi);
   treeOrig->SetBranchAddress("s9s25", &s9s25);
   treeOrig->SetBranchAddress("see", &see);
   treeOrig->SetBranchAddress("spp", &spp);
   treeOrig->SetBranchAddress("fbrem", &fbrem);
   treeOrig->SetBranchAddress("nbrem", &nbrem);
   treeOrig->SetBranchAddress("missHits", &missHits);
   treeOrig->SetBranchAddress("dist", &dist);
   treeOrig->SetBranchAddress("dcot", &dcot);

   Float_t abseta;

   // copy branches
   treeNew->Branch("mass", &mass, "mass/F");
   treeNew->Branch("pt", &pt, "pt/F");
   treeNew->Branch("eta", &eta, "eta/F");
   treeNew->Branch("abseta", &abseta, "abseta/F");
   treeNew->Branch("vertices", &vertices, "vertices/F");
   treeNew->Branch("PFIso", &LHPFIsoBasedIdOnlyIso[4], "PFIso/I");
   treeNew->Branch("ConvRej", &LHPFIsoBasedIdOnlyConv[4], "ConvRej/I");
   treeNew->Branch("CutBasedIdOlyID", CutBasedIdOlyID, "CutBasedIdOlyID[6]/I");
   treeNew->Branch("WP95", &WP95, "WP95/I");
   treeNew->Branch("WP90", &WP90, "WP90/I");
   treeNew->Branch("WP85", &WP85, "WP85/I");
   treeNew->Branch("WP80", &WP80, "WP80/I");
   treeNew->Branch("WP70", &WP70, "WP70/I");
   treeNew->Branch("WPSmurf", &WPSmurf, "WPSmurf/I");
   treeNew->Branch("BDTIdOnlyId", &BDTIdOnlyId, "BDTIdOnlyId/I");
   treeNew->Branch("BDTWP", &BDTWP, "BDTWP/I");
   treeNew->Branch("DenomFake", &DenomFake, "DenomFake/I");
   treeNew->Branch("DenomFakeSmurf", &DenomFakeSmurf, "DenomFakeSmurf/I");
   treeNew->Branch("run", &run, "run/I");
   treeNew->Branch("lumi", &lumi, "lumi/I");
   treeNew->Branch("event", &event, "event/I");
   treeNew->Branch("combPFIsoHWW", &combPFIsoHWW, "combPFIsoHWW/F");
   treeNew->Branch("chaPFIso", &chaPFIso, "chaPFIso/F");
   treeNew->Branch("neuPFIso", &neuPFIso, "neuPFIso/F");
   treeNew->Branch("phoPFIso", &phoPFIso, "phoPFIso/F");
   treeNew->Branch("npu", npu, "npu[3]/I");
   treeNew->Branch("mcmatch", &mcmatch, "mcmatch/I");
   treeNew->Branch("bdthww", &bdthww, "bdthww/F");
   treeNew->Branch("bdthzz", &bdthzz, "bdthzz/F");
   treeNew->Branch("rho", &rho, "rho/F");
   treeNew->Branch("puW", &puW, "puW/F");
   treeNew->Branch("EoPout", &EoPout, "EoPout/F");
   treeNew->Branch("EoP", &EoP, "EoP/F");
   treeNew->Branch("HoE", &HoE, "HoE/F");
   treeNew->Branch("deta", &deta, "deta/F");
   treeNew->Branch("dphi", &dphi, "dphi/F");
   treeNew->Branch("s9s25", &s9s25, "s9s25/F");
   treeNew->Branch("see", &see, "see/F");
   treeNew->Branch("spp", &spp, "spp/F");
   treeNew->Branch("fbrem", &fbrem, "fbrem/F");
   treeNew->Branch("nbrem", &nbrem, "nbrem/I");
   treeNew->Branch("missHits", &missHits, "missHits/I");
   treeNew->Branch("dist", &dist, "dist/F");
   treeNew->Branch("dcot", &dcot, "dcot/F");

   // used for PU reweighting
   PUWeight* fPUWeightFull2011 = new PUWeight("summer11","DY",-1,"Full2011",-1); 

   for(int i=0; i<nentriesOrig; i++) {
    
     if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << nentriesOrig << " entries" << std::endl;
     treeOrig->GetEntry(i);

     abseta = fabs(eta);

     BDTWP = 0;
     if(BDTIdOnlyId && LHPFIsoBasedIdOnlyIso[4] && LHPFIsoBasedIdOnlyConv[4]) BDTWP = 1;

     puW = fPUWeightFull2011->GetWeight(npu[1]);

     // customize: apply here the cut
     // if(DenomFakeSmurf==1) treeNew->Fill();
     //if(pt<20.0) treeNew->Fill();
     // if(run<170000) treeNew->Fill(); 
     treeNew->Fill();

   }
  
   fileNew->cd();
   fileNew->cd("eleIDdir");
   treeNew->Write();
   fileNew->Close();
   
   fileOrig->cd();
   fileOrig->Close();
   
  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }
}
