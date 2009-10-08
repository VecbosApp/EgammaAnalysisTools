#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>

float lumi = 10.0; // pb-1

void createAll() {
  
  double lumiZee = lumi( 2675110, 1944, 1.0); 
  double lumiWenu = lumi( 2157227, 11840, 1.0);
  double lumiWgamma = lumi( 100480, 11960, 1.0);
  double lumiTTbar = lumi( 528940, 375, 1.0);
  double lumiQCD_BCtoE_Pt20to30 = lumi( 2468398, 400e+6, 0.00048);
  double lumiQCD_BCtoE_Pt30to80 = lumi( 2041296, 100e+6, 0.0024);
  double lumiQCD_BCtoE_Pt80to170 = lumi( 1042477, 1.9e+6, 0.012);

  // --- now create the datasets useful for signal PDFs (Z -> ee ML fit) ---
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/Zee_zTandP_tree.root","datasets_ZTaP/zee.root", lumi/lumiZee);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/TTbar_zTandP_tree.root","datasets_ZTaP/ttbar.root", lumi/lumiTTbar);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/QCD_BCtoE_Pt20to30_zTandP_tree.root","datasets_ZTaP/QCD_BCtoE_Pt20to30.root", lumi/lumiQCD_BCtoE_Pt20to30);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/QCD_BCtoE_Pt30to80_zTandP_tree.root","datasets_ZTaP/QCD_BCtoE_Pt30to80.root", lumi/lumiQCD_BCtoE_Pt30to80);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/QCD_BCtoE_Pt80to170_zTandP_tree.root","datasets_ZTaP/QCD_BCtoE_Pt80to170.root", lumi/lumiQCD_BCtoE_Pt20to30);

  // merge the backgrounds to Z->ee
  TFile *ttbar = TFile::Open("datasets_ZTaP/ttbar.root");
  RooDataSet *data_ttbar = (RooDataSet*) ttbar->Get("T1");
  ttbar->Close();

  TFile *QCD_BCtoE_Pt20to30 = TFile::Open("datasets_ZTaP/QCD_BCtoE_Pt20to30.root");
  RooDataSet *data_QCD_BCtoE_Pt20to30 = (RooDataSet*) QCD_BCtoE_Pt20to30->Get("T1");
  QCD_BCtoE_Pt20to30->Close();

  TFile *QCD_BCtoE_Pt30to80 = TFile::Open("datasets_ZTaP/QCD_BCtoE_Pt30to80.root");
  RooDataSet *data_QCD_BCtoE_Pt30to80 = (RooDataSet*) QCD_BCtoE_Pt30to80->Get("T1");
  QCD_BCtoE_Pt30to80->Close();
  
  TFile *QCD_BCtoE_Pt80to170 = TFile::Open("datasets_ZTaP/QCD_BCtoE_Pt80to170.root");
  RooDataSet *data_QCD_BCtoE_Pt80to170 = (RooDataSet*) QCD_BCtoE_Pt80to170->Get("T1");
  QCD_BCtoE_Pt80to170->Close();

  RooDataSet *bkg = new RooDataSet(*data_ttbar);
  bkg->append(*data_QCD_BCtoE_Pt20to30);
  bkg->append(*data_QCD_BCtoE_Pt30to80);
  bkg->append(*data_QCD_BCtoE_Pt80to170);

  TFile *bkgFile = TFile::Open("datasets_ZTaP/background.root","recreate");
  bkg->Write();
  bkgFile->Close();

  // merge the signal + bkg to data-like dataset
  TFile *zee = TFile::Open("datasets_ZTaP/zee.root");
  RooDataSet *data_zee = (RooDataSet*) zee->Get("T1");
  zee->Close();
  
  RooDataSet *data = new RooDataSet(*data_zee);
  data->append(*bkg);
  
  TFile *datafile = TFile::Open("datasets_ZTaP/data.root","recreate");
  data->Write();
  datafile->Close();

  //  createPdfsDataset_QCDTaP("tmp/QCD_Pt15_mergedTree.root","datasets/qcd.root");

}

void createPdfsDataset_ZTaP(const char *treefile, const char *roofitfile, double weightVal=1.0) {

  gSystem->Load("libRooFit");
  
  RooRealVar *EoPout = new RooRealVar("EoPout", "EoPout", 0, 0, 100);
  RooRealVar *EoP = new RooRealVar("EoP", "EoP", 0, 0, 20);
  RooRealVar *HoE = new RooRealVar("HoE", "HoE", 0, 0, 1.0);
  RooRealVar *deltaEta = new RooRealVar("deltaEta", "deltaEta", 0, -0.02, 0.02);
  RooRealVar *deltaPhi = new RooRealVar("deltaPhi", "deltaPhi", 0, -0.1, 0.1);
  RooRealVar *s9s25 = new RooRealVar("s9s25", "s9s25", 0.7, 0.0, 1.0);
  RooRealVar *s1s9 = new RooRealVar("s1s9", "s1s9", 0.7, 0.0, 1.0);
  RooRealVar *sigmaIEtaIEta = new RooRealVar("sigmaIEtaIEta", "sigmaIEtaIEta", 0.01, 0.0, 0.1);
  RooRealVar *charge = new RooRealVar("charge", "charge", 0);
  RooRealVar *eta = new RooRealVar("eta", "eta", 0, -2.5, 2.5);
  RooRealVar *pt = new RooRealVar("pt", "pt", 0, 0., 1000,"GeV");
  RooRealVar *zmass = new RooRealVar("zmass","zmass",90,60,110,"GeV");

  RooCategory *iecal = new RooCategory("iecal", "iecal");
  iecal->defineType("barrel",0);
  iecal->defineType("endcap",1);
  RooCategory *iptbin = new RooCategory("iptbin", "iptbin");
  iptbin->defineType("lowpt",0);
  iptbin->defineType("highpt",1);
  RooCategory *iclass = new RooCategory("iclass", "iclass");
  iclass->defineType("nonshowering",0);
  iclass->defineType("showering",1); 

  RooArgSet setTagAndProbe(*EoPout,*EoP,*HoE,*deltaEta,*deltaPhi,*s9s25,*s1s9,*sigmaIEtaIEta);
  setTagAndProbe.add(*zmass);
  setTagAndProbe.add(*charge);
  setTagAndProbe.add(*eta);
  setTagAndProbe.add(*pt);
  setTagAndProbe.add(*iecal);
  setTagAndProbe.add(*iptbin);
  setTagAndProbe.add(*iclass);

  TFile *file = TFile::Open(treefile);
  TTree *tree = (TTree*)file->Get("T1");

  RooDataSet *data = new RooDataSet("T1","dataset",tree,setTagAndProbe);

  int numEntries = data->numEntries();

  RooRealVar *weight = new RooRealVar("weight","weight",0.0,100000);
  RooDataSet *weightColumn = new RooDataSet("weightData","weightData",RooArgSet(*weight));

  for(int i=0; i<numEntries; i++) {
    *weight = weightVal;
    weightColumn->add(RooArgSet(*weight));
  }

  data->merge(weightColumn);

  TFile *roofitFile = TFile::Open(roofitfile,"recreate");
  data->Write();
  roofitFile->Close();

}

void createPdfsDataset_QCDTaP(const char *treefile, const char *roofitfile, double weightVal=1.0) {

  gSystem->Load("libRooFit");
  
  RooRealVar *EoPout = new RooRealVar("EoPout", "EoPout", 0, 0, 100);
  RooRealVar *EoP = new RooRealVar("EoP", "EoP", 0, 0, 20);
  RooRealVar *HoE = new RooRealVar("HoE", "HoE", 0, 0, 1.0);
  RooRealVar *deltaEta = new RooRealVar("deltaEta", "deltaEta", 0, -0.02, 0.02);
  RooRealVar *deltaPhi = new RooRealVar("deltaPhi", "deltaPhi", 0, -0.1, 0.1);
  RooRealVar *s9s25 = new RooRealVar("s9s25", "s9s25", 0.7, 0.0, 1.0);
  RooRealVar *s1s9 = new RooRealVar("s1s9", "s1s9", 0.7, 0.0, 1.0);
  RooRealVar *sigmaIEtaIEta = new RooRealVar("sigmaIEtaIEta", "sigmaIEtaIEta", 0.01, 0.0, 0.1);
  RooRealVar *charge = new RooRealVar("qcdCharge", "qcdCharge", 0);
  RooRealVar *eta = new RooRealVar("qcdEta", "qcdEta", 0, -2.5, 2.5);
  RooRealVar *pt = new RooRealVar("qcdPt", "qcdPt", 0, 0., 1000,"GeV");
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi", "qcddeltaPhi", 0, 0., TMath::Pi());
  RooRealVar *mass = new RooRealVar("qcdInvmass","qcdInvmass",1,0,1000,"GeV");
  RooRealVar *met = new RooRealVar("qcdMet", "qcdMet", 0, 0.,100,"GeV");

  RooCategory *iecal = new RooCategory("iecal", "iecal");
  iecal->defineType("barrel",0);
  iecal->defineType("endcap",1);
  RooCategory *iptbin = new RooCategory("iptbin", "iptbin");
  iptbin->defineType("lowpt",0);
  iptbin->defineType("highpt",1);
  RooCategory *iclass = new RooCategory("iclass", "iclass");
  iclass->defineType("nonshowering",0);
  iclass->defineType("showering",1); 

  RooArgSet setTagAndProbe(*EoPout,*EoP,*HoE,*deltaEta,*deltaPhi,*s9s25,*s1s9,*sigmaIEtaIEta);
  setTagAndProbe.add(*eta);
  setTagAndProbe.add(*pt);
  setTagAndProbe.add(*deltaphi);
  setTagAndProbe.add(*mass);
  setTagAndProbe.add(*met);
  setTagAndProbe.add(*iecal);
  setTagAndProbe.add(*iptbin);
  setTagAndProbe.add(*iclass);

  TFile *file = TFile::Open(treefile);
  TTree *tree = (TTree*)file->Get("T1");

  RooDataSet *data = new RooDataSet("T1","dataset",tree,setTagAndProbe);

  int numEntries = data->numEntries();

  RooRealVar *weight = new RooRealVar("weight","weight",0.0,100000);
  RooDataSet *weightColumn = new RooDataSet("weightData","weightData",RooArgSet(*weight));

  for(int i=0; i<numEntries; i++) {
    *weight = weightVal;
    weightColumn->add(RooArgSet(*weight));
  }

  data->merge(weightColumn);

  TFile *roofitFile = TFile::Open(roofitfile,"recreate");
  data->Write();
  roofitFile->Close();

}

double lumi(double numgen, double sigma, double filtereff) {

  return numgen / sigma / filtereff; 

} 
