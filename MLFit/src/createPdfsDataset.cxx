#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>

float lumi = 10.0; // pb-1

void createAll() {
  
  double lumiZee    = lumi( 2675110, 1944, 1.0); 
  double lumiWenu   = lumi( 2157227, 11840, 1.0);
  double lumiWgamma = lumi( 100480, 11960, 1.0);
  double lumiTTbar  = lumi( 528940, 375, 1.0);
  double lumiQCD_BCtoE_Pt20to30   = lumi( 2468398, 400e+6, 0.00048);
  double lumiQCD_BCtoE_Pt30to80   = lumi( 2041296, 100e+6, 0.0024);
  double lumiQCD_BCtoE_Pt80to170  = lumi( 1042477, 1.9e+6, 0.012);
  double lumiQCD80                = lumi( 341147,  1.9e+6,  1.0 );
  double lumiPhotonJet_Pt0to15    = lumi( 109595, 100.5e+6, 1.0 );
  double lumiPhotonJet_Pt15to20   = lumi( 107480, 168.0e+3, 1.0 );
  double lumiPhotonJet_Pt20to30   = lumi( 108084, 84870,    1.0 );
  double lumiPhotonJet_Pt30to50   = lumi( 124040, 26320,    1.0 );
  double lumiPhotonJet_Pt50to80   = lumi( 112160, 4589,     1.0 );
  double lumiPhotonJet_Pt80to120  = lumi( 110000, 786.4,    1.0 );
  double lumiPhotonJet_Pt120to170 = lumi( 110000, 164.8,    1.0 );
  double lumiPhotonJet_Pt170to300 = lumi( 106143, 45.96,    1.0 );
  double lumiPhotonJet_Pt300to500 = lumi( 108290, 3.708,    1.0 );
  double lumiPhotonJet_Pt500toInf = lumi( 104735, 0.3285,   1.0 );

  /*
  // --- now create the datasets useful for signal PDFs (Z -> ee ML fit) ---
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/Zee_zTandP_tree.root","datasets_ZTaP/zee.root", lumi/lumiZee);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/TTbar_zTandP_tree.root","datasets_ZTaP/ttbar.root", lumi/lumiTTbar);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/QCD_BCtoE_Pt20to30_zTandP_tree.root","datasets_ZTaP/QCD_BCtoE_Pt20to30.root", lumi/lumiQCD_BCtoE_Pt20to30);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/QCD_BCtoE_Pt30to80_zTandP_tree.root","datasets_ZTaP/QCD_BCtoE_Pt30to80.root", lumi/lumiQCD_BCtoE_Pt30to80);
  createPdfsDataset_ZTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsSIGNAL/QCD_BCtoE_Pt80to170_zTandP_tree.root","datasets_ZTaP/QCD_BCtoE_Pt80to170.root", lumi/lumiQCD_BCtoE_Pt20to30);
  */

  // --- now create the datasets useful for background PDFs (QCD ML fit) ---
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/QCD_Pt80_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/qcd_pt80.root", lumi/lumiQCD80);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/TTbar_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root", lumi/lumiTTbar);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/Wenu_qcdTandP_5jobs_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu_5jobs.root", lumi/lumiWenu);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/Zee_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zee.root", lumi/lumiZee);  
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt0to15_qcdTandP_tree.root",   "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet015.root",    lumi/lumiPhotonJet_Pt0to15);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt15to20_qcdTandP_tree.root",  "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet1520.root",   lumi/lumiPhotonJet_Pt15to20);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt20to30_qcdTandP_tree.root",  "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet2030.root",   lumi/lumiPhotonJet_Pt20to30);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt30to50_qcdTandP_tree.root",  "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet3050.root",   lumi/lumiPhotonJet_Pt30to50);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt50to80_qcdTandP_tree.root",  "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet5080.root",   lumi/lumiPhotonJet_Pt50to80);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt80to120_qcdTandP_tree.root", "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet80120.root",  lumi/lumiPhotonJet_Pt80to120);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt120to170_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet120170.root", lumi/lumiPhotonJet_Pt120to170);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt170to300_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet170300.root", lumi/lumiPhotonJet_Pt170to300);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt300to500_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet300500.root", lumi/lumiPhotonJet_Pt300to500);
  createPdfsDataset_QCDTaP("/cmsrm/pc21/crovelli/data/Like3.2.X/resultsQCD/PhotonJet_Pt500toInf_qcdTandP_tree.root","/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet500Inf.root", lumi/lumiPhotonJet_Pt500toInf);

  /*

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


  // merge the signal + bkg to data-like dataset for Z
  TFile *zee = TFile::Open("datasets_ZTaP/zee.root");
  RooDataSet *data_zee = (RooDataSet*) zee->Get("T1");
  zee->Close();
  
  RooDataSet *data = new RooDataSet(*data_zee);
  data->append(*bkg);
  
  TFile *datafile = TFile::Open("datasets_ZTaP/data.root","recreate");
  data->Write();
  datafile->Close();
  */

  // merge some backgrounds to QCD: gammaJet + Z
  TFile *ForQCD_GammaJet_015_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet015.root");
  RooDataSet *dataQCD_GammaJet015 = (RooDataSet*) ForQCD_GammaJet_015_file->Get("T1");
  ForQCD_GammaJet_015_file->Close();

  TFile *ForQCD_GammaJet_1520_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet1520.root");
  RooDataSet *dataQCD_GammaJet1520 = (RooDataSet*) ForQCD_GammaJet_1520_file->Get("T1");
  ForQCD_GammaJet_1520_file->Close();

  TFile *ForQCD_GammaJet_2030_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet2030.root");
  RooDataSet *dataQCD_GammaJet2030 = (RooDataSet*) ForQCD_GammaJet_2030_file->Get("T1");
  ForQCD_GammaJet_2030_file->Close();

  TFile *ForQCD_GammaJet_3050_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet3050.root");
  RooDataSet *dataQCD_GammaJet3050 = (RooDataSet*) ForQCD_GammaJet_3050_file->Get("T1");
  ForQCD_GammaJet_3050_file->Close();

  TFile *ForQCD_GammaJet_5080_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet5080.root");
  RooDataSet *dataQCD_GammaJet5080 = (RooDataSet*) ForQCD_GammaJet_5080_file->Get("T1");
  ForQCD_GammaJet_5080_file->Close();

  TFile *ForQCD_GammaJet_80120_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet80120.root");
  RooDataSet *dataQCD_GammaJet80120 = (RooDataSet*) ForQCD_GammaJet_80120_file->Get("T1");
  ForQCD_GammaJet_80120_file->Close();

  TFile *ForQCD_GammaJet_120170_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet120170.root");
  RooDataSet *dataQCD_GammaJet120170 = (RooDataSet*) ForQCD_GammaJet_120170_file->Get("T1");
  ForQCD_GammaJet_120170_file->Close();

  TFile *ForQCD_GammaJet_170300_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet170300.root");
  RooDataSet *dataQCD_GammaJet170300 = (RooDataSet*) ForQCD_GammaJet_170300_file->Get("T1");
  ForQCD_GammaJet_170300_file->Close();

  TFile *ForQCD_GammaJet_300500_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet300500.root");
  RooDataSet *dataQCD_GammaJet300500 = (RooDataSet*) ForQCD_GammaJet_300500_file->Get("T1");
  ForQCD_GammaJet_300500_file->Close();

  TFile *ForQCD_GammaJet_500Inf_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/photonJet500Inf.root");
  RooDataSet *dataQCD_GammaJet500Inf = (RooDataSet*) ForQCD_GammaJet_500Inf_file->Get("T1");
  ForQCD_GammaJet_500Inf_file->Close();

  TFile *ForQCD_Zee_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zee.root");
  RooDataSet *dataQCD_Zee = (RooDataSet*) ForQCD_Zee_file->Get("T1");
  ForQCD_Zee_file->Close();

  RooDataSet *bkgZandPhotonJetForQcd = new RooDataSet(*dataQCD_Zee);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet015);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet1520);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet2030);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet3050);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet5080);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet80120);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet120170);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet170300);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet300500);
  bkgZandPhotonJetForQcd->append(*dataQCD_GammaJet500Inf);

  TFile *bkgZandPhotonJetForQcd_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zAndPhotonJet.root","recreate");
  bkgZandPhotonJetForQcd->Write();
  bkgZandPhotonJetForQcd_file->Close();

  
  // merge the signal + bkg to data-like dataset for QCD
  TFile *qcdPt80ForQcd_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/qcd_pt80.root");
  RooDataSet *dataQcd_qcd80 = (RooDataSet*) qcdPt80ForQcd_file->Get("T1");
  qcdPt80ForQcd_file->Close();

  TFile *ttbarForQcd_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root");
  RooDataSet *dataQcd_ttbar = (RooDataSet*) ttbarForQcd_file->Get("T1");
  ttbarForQcd_file->Close();  

  TFile *wenuForQcd_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu_5jobs.root");
  RooDataSet *dataQcd_wenu = (RooDataSet*) wenuForQcd_file->Get("T1");
  wenuForQcd_file->Close();    

  RooDataSet *dataForQcd = new RooDataSet(*dataQcd_qcd80);
  dataForQcd->append(*bkgZandPhotonJetForQcd);
  dataForQcd->append(*dataQcd_ttbar);
  dataForQcd->append(*dataQcd_wenu);
  
  TFile *datafileForQcd_file = TFile::Open("/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/data.root","recreate");
  dataForQcd->Write();
  datafileForQcd_file->Close();
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
  RooRealVar *zmass = new RooRealVar("zmass","zmass", 90, 60, 110, "GeV");
  RooRealVar *iecal = new RooRealVar("iecal", "iecal", 0, 1);
  RooRealVar *iptbin = new RooRealVar("iptbin", "iptbin", 0, 1);
  RooRealVar *iclass = new RooRealVar("iclass", "iclass", 0, 1);

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
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi", "qcddeltaPhi", 0, 0.5, 3.1415);
  RooRealVar *mass = new RooRealVar("qcdInvmass","qcdInvmass",1,0,600,"GeV");
  RooRealVar *met = new RooRealVar("qcdMet", "qcdMet", 0, 0.,100,"GeV");
  RooRealVar *iecal = new RooRealVar("iecal", "iecal", 0, 1);
  RooRealVar *iptbin = new RooRealVar("iptbin", "iptbin", 0, 1);
  RooRealVar *iclass = new RooRealVar("iclass", "iclass", 0, 1);

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
