#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>

void createPdfsDataset(const char *treefile, const char *roofitfile) {

  gSystem->Load("libRooFit");
  
  RooRealVar *EoPout = new RooRealVar("EoPout", "EoPout", 0, 0, 8);
  RooRealVar *HoE = new RooRealVar("HoE", "HoE", 0, 0, 0.1);
  RooRealVar *deltaEta = new RooRealVar("deltaEta", "deltaEta", 0, -0.02, 0.02);
  RooRealVar *deltaPhi = new RooRealVar("deltaPhi", "deltaPhi", 0, -0.1, 0.1);
  RooRealVar *s9s25 = new RooRealVar("s9s25", "s9s25", 0.7, 0.5, 1.0);
  RooRealVar *sigmaEtaEta = new RooRealVar("sigmaEtaEta", "sigmaEtaEta", 0.01, 0.0, 0.05);
  RooRealVar *charge = new RooRealVar("charge", "charge", 0);
  RooRealVar *eta = new RooRealVar("eta", "eta", 0, -2.5, 2.5);
  RooRealVar *pt = new RooRealVar("pt", "pt", 0, 0., 1000);
  RooRealVar *zmass = new RooRealVar("zmass","zmass",90,40,120);

//   RooCategory *iecal = new RooCategory("iecal", "iecal");
//   iecal->defineType("barrel",0);
//   iecal->defineType("endcap",1);
//   RooCategory *iptbin = new RooCategory("iptbin", "iptbin");
//   iptbin->defineType("lowpt",0);
//   iptbin->defineType("highpt",1);
//   RooCategory *iclass = new RooCategory("iclass", "iclass");
//   iclass->defineType("nonshowering",0);
//   iclass->defineType("showering",1); 

  RooArgSet setTagAndProbe(*EoPout,*HoE,*deltaEta,*deltaPhi,*s9s25,*sigmaEtaEta,*zmass);
  setTagAndProbe.add(*charge);
  setTagAndProbe.add(*eta);
  setTagAndProbe.add(*pt);
//   setTagAndProbe.add(*iecal);
//   setTagAndProbe.add(*iptbin);
//   setTagAndProbe.add(*iclass);

  TFile *file = TFile::Open(treefile);
  TTree *tree = (TTree*)file->Get("T1");

  RooDataSet *data = new RooDataSet("T1","dataset",tree,setTagAndProbe);

  TFile *roofitFile = TFile::Open(roofitfile,"recreate");
  data->Write();
  roofitFile->Close();

}
