#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include "../PUWeight.C"

using namespace std;

enum idType {
  kIsoHWW2011 = 0, // HWW cuts 2011
  kIsoEACorr, // EA corrected isolation HZZ with WP equal to kIsoHWW2011
  kBDTHWW2011_withIP, // HWW cuts 2011
  kBDTHWW2011_noIP, // HWW BDT w/o IP (not optimized cuts)
  kBDTHZZ_withIP, // HZZ BDT with IP
  kBDTHZZ_noIP 
};

bool passID(Float_t eta, Float_t pt, Float_t bdthww, Float_t bdthzz, Float_t rho, Float_t combIso, Float_t combPFIsoHWW, idType type) {
  if(type == kIsoHWW2011) {
    if(fabs(eta)<1.479) return (combPFIsoHWW/pt < 0.13);
    else return (combPFIsoHWW/pt < 0.09);
  }

  if(type == kIsoEACorr) {

    if(fabs(eta) <  1.0) combIso -= 0.18 * rho;
    if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) combIso -= 0.19 * rho;
    if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) combIso -= 0.21 * rho;
    if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) combIso -= 0.38 * rho;
    if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) combIso -= 0.61 * rho;
    if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) combIso -= 0.73 * rho;
    if(fabs(eta) >=  2.4) combIso -= 0.90 * rho;

    if(pt>20) {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.23); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.20);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.12);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.11);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.049);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.070);
      if(fabs(eta) >=  2.4) return (combIso/pt < 0.010);
    } else {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.20); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.21);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.13);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.10);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.04);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.03);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.03);
    }
  }

  if(type == kBDTHWW2011_withIP) {
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.884);
  }

  if(type == kBDTHZZ_withIP) {
    // WP with same fake rate as HWW with IP
    if(pt < 20 && fabs(eta) < 1.0) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.091);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthzz > 0.064);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.071);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.067);
  }

  return false;
}

void makeFitTree(const char* filename, float weight) {
  cout << "Adding weight and ID bits to file " << filename << " with weight " << weight << endl;

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

  TFile *fileNew = TFile::Open("tptree.root","recreate");
  fileNew->mkdir("eleIDdir");
  TTree *treeNew = new TTree("T1","tree with id bits");

   // Declaration of leaf types
  Float_t         mass;
  Int_t           ecaldriven;
   Float_t         pt;
   Float_t         eta;
   Float_t         phi;
   Int_t           DenomFakeSmurf;
   Float_t         combPFIsoHWW;
   Float_t         chaPFIso;
   Float_t         neuPFIso;
   Float_t         phoPFIso;
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Int_t           mcmatch;
   Float_t         bdthww;
   Float_t         bdthwwnoip;
   Float_t         bdthzz;
   Float_t         bdthzznoip;
   Float_t         bdthzzmc;
   Float_t         vertices;
   Float_t         rho;
   Int_t           npu[3];

   treeOrig->SetBranchAddress("mass", &mass);
   treeOrig->SetBranchAddress("ecaldriven", &ecaldriven);
   treeOrig->SetBranchAddress("pt", &pt);
   treeOrig->SetBranchAddress("eta", &eta);
   treeOrig->SetBranchAddress("phi", &phi);
   treeOrig->SetBranchAddress("DenomFakeSmurf", &DenomFakeSmurf);
   treeOrig->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
   treeOrig->SetBranchAddress("chaPFIso", &chaPFIso);
   treeOrig->SetBranchAddress("neuPFIso", &neuPFIso);
   treeOrig->SetBranchAddress("phoPFIso", &phoPFIso);
   treeOrig->SetBranchAddress("run", &run);
   treeOrig->SetBranchAddress("lumi", &lumi);
   treeOrig->SetBranchAddress("event", &event);
   treeOrig->SetBranchAddress("mcmatch", &mcmatch);
   treeOrig->SetBranchAddress("bdthww", &bdthww);
   treeOrig->SetBranchAddress("bdthwwnoip", &bdthwwnoip);
   treeOrig->SetBranchAddress("bdthzz", &bdthzz);
   treeOrig->SetBranchAddress("bdthzznoip", &bdthzznoip);
   treeOrig->SetBranchAddress("bdthzzmc", &bdthzzmc);
   treeOrig->SetBranchAddress("vertices", &vertices);
   treeOrig->SetBranchAddress("rho", &rho);
   treeOrig->SetBranchAddress("npu", npu);

   // the bits
   Float_t abseta, puW;
   Int_t pass_bdthww, pass_bdthzz;
   treeNew->Branch("mass", &mass, "mass/F");
   treeNew->Branch("ecaldriven", &ecaldriven, "ecaldriven/I");
   treeNew->Branch("pt", &pt, "pt/F");
   treeNew->Branch("abseta", &abseta, "abseta/F");
   treeNew->Branch("phi", &phi, "phi/F");
   treeNew->Branch("vertices", &vertices, "vertices/I");
   treeNew->Branch("bdthww", &pass_bdthww, "bdthww/I");
   treeNew->Branch("bdthzz", &pass_bdthzz, "bdthzz/I");
   treeNew->Branch("puW", &puW, "puW/F");

   // used for PU reweighting
   PUWeight* fPUWeightFull2011 = new PUWeight("summer11","DY",-1,"Full2011",-1); 
   
   for(int i=0; i<nentriesOrig; i++) {
    
     if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << nentriesOrig << " entries" << std::endl;
     treeOrig->GetEntry(i);

     abseta = fabs(eta);
     puW = fPUWeightFull2011->GetWeight(npu[1]);

     Float_t combIso = chaPFIso+neuPFIso+phoPFIso;
     pass_bdthww=0;
     if(passID(eta,pt,bdthww,bdthzz,rho,combIso,combPFIsoHWW,kBDTHWW2011_withIP) && 
	passID(eta,pt,bdthww,bdthzz,rho,combIso,combPFIsoHWW,kIsoHWW2011)) pass_bdthww = 1;
     pass_bdthzz=0;
     if(passID(eta,pt,bdthww,bdthzz,rho,combIso,combPFIsoHWW,kBDTHZZ_withIP) && 
	passID(eta,pt,bdthww,bdthzz,rho,combIso,combPFIsoHWW,kIsoEACorr)) pass_bdthzz = 1;

     // customize: apply here the cut
     // if(DenomFakeSmurf==1) treeNew->Fill();
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
