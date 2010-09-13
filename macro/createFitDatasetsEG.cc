#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>

#include <iostream>

using namespace std;

void createSingleDataset(const char *treefile, const char *roofitfile);

void createAll() {

  cout << "CREATING ROODATASETS..." << endl;
  createSingleDataset("results/trees_reduced/WJetsMADGRAPH_tree_50k.root","results/datasets/WJetsMADGRAPH_tree_50k.root");
  cout << "done with W" << endl;
  createSingleDataset("results/trees_reduced/ZJetsMADGRAPH_tree_fullStat.root","results/datasets/ZJetsMADGRAPH_tree_fullStat.root");
  cout << "done with Z" << endl;
  createSingleDataset("results/trees_reduced/GammaJets_tree_50k.root",    "results/datasets/GammaJets_tree_50k.root");
  cout << "done with gamma + jets" << endl;
  createSingleDataset("results/trees_reduced/QCDandGJ_tree_50k.root",     "results/datasets/QCDandGJ_tree_50k.root");
  cout << "done with QCD and gamma+jets" << endl;
  createSingleDataset("results/trees_reduced/QCD_tree_50k.root",          "results/datasets/QCD_tree_50k.root");
  cout << "done with QCD" << endl;
  createSingleDataset("results/trees_reduced/TOP_tree_50k.root",          "results/datasets/TOP_tree_50k.root");
  cout << "done with top" << endl;

  // createSingleDataset("results_data/data_Wenu.root",              "results/datasets/data_Wenu"); 
}

void createSingleDataset(const char *treefile, const char *roofitfile) {

  cout << "roofitting file " << treefile << " in " << roofitfile << endl;

  // kinematic and ID variables
  RooRealVar *EoPout = new RooRealVar("EoPout", "EoPout", 0, 0, 1000);
  RooRealVar *EoP    = new RooRealVar("EoP",    "EoP",    0, 0, 1000);
  RooRealVar *HoE    = new RooRealVar("HoE",    "HoE",    0, -100, 100.0);
  RooRealVar *deltaEta      = new RooRealVar("deltaEta",      "deltaEta",      0, -1000., 1000.);
  RooRealVar *deltaEtaCorr  = new RooRealVar("deltaEtaCorr",  "deltaEtaCorr",  0, -1000., 1000.);
  RooRealVar *deltaPhi      = new RooRealVar("deltaPhi",      "deltaPhi",      0, -1000., 1000.);
  RooRealVar *deltaPhiCorr  = new RooRealVar("deltaPhiCorr",  "deltaPhiCorr",  0, -1000., 1000.);
  RooRealVar *sigmaIEtaIEta = new RooRealVar("sigmaIEtaIEta", "sigmaIEtaIEta", 0.01, 0.0, 1000.);
  RooRealVar *charge = new RooRealVar("qcdCharge", "qcdCharge", 0);
  RooRealVar *eta    = new RooRealVar("qcdEta", "qcdEta", 0, -5., 5.);
  RooRealVar *pt     = new RooRealVar("qcdPt", "qcdPt", 0, 0., 1000,"GeV");
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi", "qcddeltaPhi", 0, 0., 3.1415);
  RooRealVar *mass     = new RooRealVar("qcdInvmass","qcdInvmass",1,0,150,"GeV");
  RooRealVar *met      = new RooRealVar("qcdMet", "qcdMet", 0, 0.,100,"GeV");
  RooRealVar *iecal  = new RooRealVar("iecal",  "iecal",  0, 1);
  RooRealVar *iptbin = new RooRealVar("iptbin", "iptbin", 0, 1);

  // other variables - chiara, commenta per dati
  RooRealVar *weight   = new RooRealVar("weight","weight",0,1000);

  RooArgSet setTagAndProbe(*EoPout,*EoP,*HoE,*deltaEta,*deltaPhi,*deltaEtaCorr,*deltaPhiCorr,*sigmaIEtaIEta);
  setTagAndProbe.add(*eta);
  setTagAndProbe.add(*pt);
  setTagAndProbe.add(*charge);
  setTagAndProbe.add(*deltaphi);
  setTagAndProbe.add(*mass);
  setTagAndProbe.add(*met);
  setTagAndProbe.add(*iecal);
  setTagAndProbe.add(*iptbin);
  setTagAndProbe.add(*weight); // chiara, commenta x dati

  TFile *file = TFile::Open(treefile);
  TTree *tree = (TTree*)file->Get("T1");
  
  RooDataSet *data = new RooDataSet("T1","dataset",tree,setTagAndProbe);
  
  TFile *roofitFile = TFile::Open(roofitfile,"recreate");
  data->Write();
  roofitFile->Close();
}

