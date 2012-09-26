#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"

#include <sstream>
#include <iostream>

using namespace std;
using namespace RooFit;

class RhhCruijffPdf : public RooAbsPdf {
public:
  RhhCruijffPdf() { } ;
  RhhCruijffPdf(const char *name, const char *title, RooAbsReal& _m,
		RooAbsReal& _m0, 
		RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
		RooAbsReal& _alphaL, RooAbsReal& _alphaR) ;
  
  RhhCruijffPdf(const RhhCruijffPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new RhhCruijffPdf(*this,newname); }

  inline virtual ~RhhCruijffPdf() { }

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy alphaR;

  Double_t evaluate() const;

private:
  
  ClassDef(RhhCruijffPdf,0)
};

RhhCruijffPdf::RhhCruijffPdf(const char *name, const char *title,
	     RooAbsReal& _m, RooAbsReal& _m0, 
	     RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
	     RooAbsReal& _alphaL, RooAbsReal& _alphaR)
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigmaL("sigmaL", "SigmaL", this, _sigmaL),
  sigmaR("sigmaR", "SigmaR", this, _sigmaR),
  alphaL("alphaL", "AlphaL", this, _alphaL),
  alphaR("alphaR", "AlphaR", this, _alphaR)
{
}

RhhCruijffPdf::RhhCruijffPdf(const RhhCruijffPdf& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigmaL("sigmaL", this, other.sigmaL), sigmaR("sigmaR", this, other.sigmaR), 
  alphaL("alphaL", this, other.alphaL), alphaR("alphaR", this, other.alphaR)
{
}

Double_t RhhCruijffPdf::evaluate() const 
{
  double dx = (m-m0) ;
  double sigma = dx<0 ? sigmaL: sigmaR ;
  double alpha = dx<0 ? alphaL: alphaR ;
  double f = 2*sigma*sigma + alpha*dx*dx ;
  return exp(-dx*dx/f) ;
}


Double_t cryball( Double_t *x, Double_t * par) {

  // par[0] = normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n


  double cryball;
  double test = par[1] + par[2]*par[3];
  double xp = (x[0] - par[1]) / par[2];
  double A = TMath::Power((par[4]/par[3]),par[4]) * TMath::Exp(-0.5*par[3]*par[3]);
  double  B = ( par[4]/par[3] ) - par[3];


  if (x[0] < test)

    {
      cryball = par[0]*TMath::Exp(-0.5*xp*xp);

    } else {

      cryball =par[0]* A * TMath::Power(( B + xp ),-par[4]);

   }
  return cryball;
}

Double_t cruijff( Double_t *x, Double_t * par) {

  // par[0] = normalization
  // par[1] = mean
  // par[2] = sigmaL
  // par[3] = sigmaR
  // par[4] = alphaL
  // par[5] = alphaR

  double dx = (x[0]-par[1]) ;
  double sigma = dx<0 ? par[2]: par[3] ;
  double alpha = dx<0 ? par[4]: par[5] ;
  double f = 2*sigma*sigma + alpha*dx*dx ;
  return par[0] * exp(-dx*dx/f) ;

}

void makeResolutionFriendTree(const char *file) {
  TFile *pF = TFile::Open(file);
  TTree *pT = (TTree*)pF->Get("electronTree/probe_tree");
  
  float ecalE, scE, scrawE, p, genp;
  float pt, eta, rho, classification;

  pT->SetBranchAddress("ecalE", &ecalE);
  pT->SetBranchAddress("scE", &scE);
  pT->SetBranchAddress("scrawE", &scrawE);
  pT->SetBranchAddress("p", &p);
  pT->SetBranchAddress("genp", &genp);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("classification", &classification);
  
  TString nF(file);
  nF.ReplaceAll(".root","_resolutionFriend.root");
  TFile *fF = TFile::Open(nF,"recreate");

  fF->mkdir("electronTree");
  TTree *fT = new TTree("probe_tree","tree with energy resolution");

  float resecalE, resscE, resscrawE, resp;
  fT->Branch("resecalE", &resecalE);
  fT->Branch("resscE", &resscE);
  fT->Branch("resscrawE", &resscrawE);
  fT->Branch("resp", &resp);
  fT->Branch("genp", &genp);
  fT->Branch("pt", &pt);
  fT->Branch("eta", &eta);
  fT->Branch("rho", &rho);
  fT->Branch("classification", &classification);

 
 for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);
     resecalE=(ecalE-genp)/genp;
     resscE=(scE-genp)/genp;
     resscrawE=(scrawE-genp)/genp;
     resp=(p-genp)/genp;
     fT->Fill();
 }

 fF->cd("electronTree");
 fT->Write();
 fF->Close();

 cout << "DONE. Friend tree is in file: " << nF.Data() << endl;


}

void makeDependencyPlot() {

 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  TFile *file = TFile::Open("/Users/emanuele/Work/data/hzz4l/electronreg/HZZ4L_53X_S1_V10_S2_V01/MC/EleRegr1/ggH125.root");
  TTree *tree = (TTree*)file->Get("electronTree/probe_tree");

  TH1F *EoEt = new TH1F("EoEt","",101,-0.5,0.5);
  EoEt->GetXaxis()->SetTitle("(E-p_{true})/p_{true}");
  
  float ptbins[12] = {7,10,15,20,25,30,35,40,45,50,70,100};
  float etabins[7] = {0,0.5,0.8,1.2,1.479,2.0,2.5};
  float vtxbins[11] = {0,5,10,15,17,19,21,23,25,30,50};
  float classificationbins[5] = {0,1,2,3,4};
  
  // mean
  TH1F *ptM = new TH1F("ptM","",11,ptbins);
  TH1F *etaM = new TH1F("etaM","",6,etabins);
  TH1F *vtxM = new TH1F("vtxM","",10,vtxbins);
  TH1F *classM = new TH1F("classM","",4,classificationbins);
 
  ptM->GetXaxis()->SetTitle("p_{T} [GeV]");
  etaM->GetXaxis()->SetTitle("#eta");
  vtxM->GetXaxis()->SetTitle("#rho");
  classM->GetXaxis()->SetTitle("class");

  ptM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");
  etaM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");
  vtxM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");
  classM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");

  ptM->GetYaxis()->SetTitleOffset(1.8);
  etaM->GetYaxis()->SetTitleOffset(1.8);
  vtxM->GetYaxis()->SetTitleOffset(1.8);
  classM->GetYaxis()->SetTitleOffset(1.8);

  ptM->SetMarkerStyle(8);
  ptM->SetMarkerSize(1);
  etaM->SetMarkerStyle(8);
  etaM->SetMarkerSize(1);
  vtxM->SetMarkerStyle(8);
  vtxM->SetMarkerSize(1);
  classM->SetMarkerStyle(8);
  classM->SetMarkerSize(1);


  ptM->SetMinimum(-0.05);
  ptM->SetMaximum(0.05);


  // peak
  TH1F *ptP = (TH1F*)ptM->Clone("ptP");
  TH1F *etaP = (TH1F*)etaM->Clone("etaP");
  TH1F *vtxP = (TH1F*)vtxM->Clone("vtxP");
  TH1F *classP = (TH1F*)classM->Clone("classP");
  ptP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");
  etaP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");
  vtxP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");
  classP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");

  ptP->SetMinimum(-0.05);
  ptP->SetMaximum(0.05);
 
  // RMS
  TH1F *ptRMS = (TH1F*)ptM->Clone("ptRMS");
  TH1F *etaRMS = (TH1F*)etaM->Clone("etaRMS");
  TH1F *vtxRMS = (TH1F*)vtxM->Clone("vtxRMS");
  TH1F *classRMS = (TH1F*)classM->Clone("classRMS");
  ptRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  etaRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  vtxRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  classRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  ptRMS->SetMinimum(0);
  ptRMS->SetMaximum(0.2);


  // Sigma
  TH1F *ptSigma = (TH1F*)ptM->Clone("ptSigma");
  TH1F *etaSigma = (TH1F*)etaM->Clone("etaSigma");
  TH1F *vtxSigma = (TH1F*)vtxM->Clone("vtxSigma");
  TH1F *classSigma = (TH1F*)classM->Clone("classSigma");
  ptSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  etaSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  vtxSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  classSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  ptSigma->SetMinimum(0);
  ptSigma->SetMaximum(0.2);

  TF1 *func = new TF1("cruijff",cruijff,-1,1,6); 

  TCanvas *c1 = new TCanvas("c1","c1",600,600);


  
  // ============ PT ==============
  cout << "===> RUNNING VS PT " << endl;
  for(int i=0;i<11;++i) {
    stringstream cut;
    cut << "pt>" << ptbins[i] << "&& pt<" << ptbins[i+1] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_pt_" << ptbins[i] << "To" << ptbins[i+1] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(1,-0.05,0.05);
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.08,0.08);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    ptM->SetBinContent(i+1,mean);    
    ptM->SetBinError(i+1,meanerr);    
    ptP->SetBinContent(i+1,peak);    
    ptP->SetBinError(i+1,peakerr);    
    ptRMS->SetBinContent(i+1,rms);    
    ptRMS->SetBinError(i+1,rmserr);    
    ptSigma->SetBinContent(i+1,sigma);    
    ptSigma->SetBinError(i+1,sigmaerr);    

  }
  
  ptM->Draw();  c1->SaveAs("ptM.png");
  ptP->Draw();  c1->SaveAs("ptP.png");
  ptRMS->Draw();  c1->SaveAs("ptRMS.png");
  ptSigma->Draw();  c1->SaveAs("ptSigma.png");
  




  // ============ ETA ==============
  cout << "===> RUNNING VS ETA " << endl;
  etaM->SetMaximum(0.1);
  etaP->SetMaximum(0.1);
  etaRMS->SetMaximum(0.1);
  etaSigma->SetMaximum(0.1);

  for(int i=0;i<6;++i) {
    stringstream cut;
    cut << "abs(eta)>" << etabins[i] << "&& abs(eta)<" << etabins[i+1] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_eta_" << etabins[i] << "To" << etabins[i+1] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(1,-0.05,0.05);
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.08,0.08);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    etaM->SetBinContent(i+1,mean);    
    etaM->SetBinError(i+1,meanerr);    
    etaP->SetBinContent(i+1,peak);    
    etaP->SetBinError(i+1,peakerr);    
    etaRMS->SetBinContent(i+1,rms);    
    etaRMS->SetBinError(i+1,rmserr);    
    etaSigma->SetBinContent(i+1,sigma);    
    etaSigma->SetBinError(i+1,sigmaerr);    

  }
  
  etaM->Draw();  c1->SaveAs("etaM.png");
  etaP->Draw();  c1->SaveAs("etaP.png");
  etaRMS->Draw();  c1->SaveAs("etaRMS.png");
  etaSigma->Draw();  c1->SaveAs("etaSigma.png");



  // ============ VTX ==============
  cout << "===> RUNNING VS NVTX " << endl;
  vtxM->SetMaximum(0.05);
  vtxP->SetMaximum(0.05);
  vtxM->SetMinimum(-0.05);
  vtxP->SetMinimum(-0.05);
  vtxRMS->SetMaximum(0.1);
  vtxSigma->SetMaximum(0.05);

  for(int i=0;i<10;++i) {
    stringstream cut;
    cut << "rho>" << vtxbins[i] << "&& rho<" << vtxbins[i+1] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_vtx_" << vtxbins[i] << "To" << vtxbins[i+1] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(1,-0.05,0.05);
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.08,0.08);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    vtxM->SetBinContent(i+1,mean);    
    vtxM->SetBinError(i+1,meanerr);    
    vtxP->SetBinContent(i+1,peak);    
    vtxP->SetBinError(i+1,peakerr);    
    vtxRMS->SetBinContent(i+1,rms);    
    vtxRMS->SetBinError(i+1,rmserr);    
    vtxSigma->SetBinContent(i+1,sigma);    
    vtxSigma->SetBinError(i+1,sigmaerr);    

  }
  
  vtxM->Draw();  c1->SaveAs("vtxM.png");
  vtxP->Draw();  c1->SaveAs("vtxP.png");
  vtxRMS->Draw();  c1->SaveAs("vtxRMS.png");
  vtxSigma->Draw();  c1->SaveAs("vtxSigma.png");





  // ============ ELE CLASS ==============
  cout << "===> RUNNING VS ELE CLASS " << endl;
  classM->SetMaximum(0.5);
  classP->SetMaximum(0.05);
  classRMS->SetMaximum(0.2);
  classSigma->SetMaximum(0.05);

  for(int i=0;i<4;++i) {
    stringstream cut;
    cut << "classification==" << classificationbins[i] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_class_" << classificationbins[i] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(1,-0.05,0.05);
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.08,0.08);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    classM->SetBinContent(i+1,mean);    
    classM->SetBinError(i+1,meanerr);    
    classP->SetBinContent(i+1,peak);    
    classP->SetBinError(i+1,peakerr);    
    classRMS->SetBinContent(i+1,rms);    
    classRMS->SetBinError(i+1,rmserr);    
    classSigma->SetBinContent(i+1,sigma);    
    classSigma->SetBinError(i+1,sigmaerr);    

  }
  
  classM->Draw();  c1->SaveAs("classM.png");
  classP->Draw();  c1->SaveAs("classP.png");
  classRMS->Draw();  c1->SaveAs("classRMS.png");
  classSigma->Draw();  c1->SaveAs("classSigma.png");



  
  TFile *resultfile = TFile::Open("results_elereg_ggH125.root","recreate");
  ptM->Write(); ptP->Write(); ptRMS->Write(); ptSigma->Write();
  etaM->Write(); etaP->Write(); etaRMS->Write(); etaSigma->Write();
  vtxM->Write(); vtxP->Write(); vtxRMS->Write(); vtxSigma->Write();
  classM->Write(); classP->Write(); classRMS->Write(); classSigma->Write();
  resultfile->Close();


}


void makeDependencyPlotConv() {

 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  TFile *file = TFile::Open("/Users/emanuele/Work/data/hzz4l/electronreg/HZZ4L_53X_S1_V10_S2_V01/MC/EleRegr1/ggH125_resolutionFriend.root");
  TTree *tree = (TTree*)file->Get("electronTree/probe_tree");

  RooRealVar x("resecalE","E resolution",-0.5,0.5);
  RooRealVar genp("genp","genp",-100,10000,"GeV");
  RooRealVar pt("pt","pt",7,100,"GeV");
  RooRealVar eta("eta","eta",-2.5,2.5);
  RooRealVar rho("rho","rho",-0.5,50);
  RooRealVar classification("classification","classification",-0.5,4.5);

  RooArgSet varset(x,genp,pt,eta,rho,classification);

  float ptbins[12] = {7,10,15,20,25,30,35,40,45,50,70,100};
  float etabins[7] = {0,0.5,0.8,1.2,1.479,2.0,2.5};
  float vtxbins[11] = {0,5,10,15,17,19,21,23,25,30,50};
  float classificationbins[5] = {0,1,2,3,4};
  
  TH1F *EoEt = new TH1F("EoEt","",100,-0.5,0.5);

  // mean
  TH1F *ptM = new TH1F("ptM","",11,ptbins);
  TH1F *etaM = new TH1F("etaM","",6,etabins);
  TH1F *vtxM = new TH1F("vtxM","",10,vtxbins);
  TH1F *classM = new TH1F("classM","",4,classificationbins);
 
  ptM->GetXaxis()->SetTitle("p_{T} [GeV]");
  etaM->GetXaxis()->SetTitle("#eta");
  vtxM->GetXaxis()->SetTitle("#rho");
  classM->GetXaxis()->SetTitle("class");

  ptM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");
  etaM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");
  vtxM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");
  classM->GetYaxis()->SetTitle("(E-p_{true})/p_{true} mean");

  ptM->GetYaxis()->SetTitleOffset(1.8);
  etaM->GetYaxis()->SetTitleOffset(1.8);
  vtxM->GetYaxis()->SetTitleOffset(1.8);
  classM->GetYaxis()->SetTitleOffset(1.8);

  ptM->SetMarkerStyle(8);
  ptM->SetMarkerSize(1);
  etaM->SetMarkerStyle(8);
  etaM->SetMarkerSize(1);
  vtxM->SetMarkerStyle(8);
  vtxM->SetMarkerSize(1);
  classM->SetMarkerStyle(8);
  classM->SetMarkerSize(1);


  ptM->SetMinimum(-0.05);
  ptM->SetMaximum(0.05);


  // peak
  TH1F *ptP = (TH1F*)ptM->Clone("ptP");
  TH1F *etaP = (TH1F*)etaM->Clone("etaP");
  TH1F *vtxP = (TH1F*)vtxM->Clone("vtxP");
  TH1F *classP = (TH1F*)classM->Clone("classP");
  ptP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");
  etaP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");
  vtxP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");
  classP->GetYaxis()->SetTitle("(E-p_{true})/p_{true} peak");

  ptP->SetMinimum(-0.05);
  ptP->SetMaximum(0.05);
 
  // RMS
  TH1F *ptRMS = (TH1F*)ptM->Clone("ptRMS");
  TH1F *etaRMS = (TH1F*)etaM->Clone("etaRMS");
  TH1F *vtxRMS = (TH1F*)vtxM->Clone("vtxRMS");
  TH1F *classRMS = (TH1F*)classM->Clone("classRMS");
  ptRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  etaRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  vtxRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  classRMS->GetYaxis()->SetTitle("(E-p_{true})/p_{true} RMS");
  ptRMS->SetMinimum(0);
  ptRMS->SetMaximum(0.2);


  // Sigma
  TH1F *ptSigma = (TH1F*)ptM->Clone("ptSigma");
  TH1F *etaSigma = (TH1F*)etaM->Clone("etaSigma");
  TH1F *vtxSigma = (TH1F*)vtxM->Clone("vtxSigma");
  TH1F *classSigma = (TH1F*)classM->Clone("classSigma");
  ptSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  etaSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  vtxSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  classSigma->GetYaxis()->SetTitle("(E-p_{true})/p_{true} #sigma");
  ptSigma->SetMinimum(0);
  ptSigma->SetMaximum(0.2);

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  
  // ============ PT ==============
  cout << "===> RUNNING VS PT " << endl;
  for(int i=0;i<11;++i) {

    //--- simple CrystalBall
    RooRealVar mean1("mean1","mean of gaussian",0.0,-0.1,0.1) ;
    RooRealVar sigma1("sigma1","width of gaussian",0.01,0.005,0.10); 
    RooRealVar sigma2("sigma2","width of gaussian",0.01,0.005,0.10); 
    RooRealVar a("a","a",1.46,0.,10.);
    RooRealVar n("n","n",1.92,0.,20.);   
    RooRealVar alphaL("alphaL","alphaL",1.46,0.,10.);
    RooRealVar alphaR("alphaR","alphaR",1.46,0.,10.);
    RooCBShape CBall("CBall","Crystal ball",x, mean1,sigma1, a,n);
    RhhCruijffPdf cruijff("Cruijff","cruijff",x,mean1,sigma1,sigma2,alphaL,alphaR);
    
    //--- Breit-Wigner
    RooRealVar mean3("mean3","mean3",0.0) ;
    RooRealVar sigma3("sigma3","width3",0.005); 
    RooBreitWigner bw("bw","bw",x,mean3,sigma3);
    
    x.setBins(10000,"fft");
    RooFFTConvPdf model("model","model",x,bw,cruijff);

    stringstream cut;
    cut << "pt>" << ptbins[i] << "&& pt<" << ptbins[i+1] << "&& abs(resecalE)<0.5 && genp>0";
    TCut thecut(cut.str().c_str());

    RooDataSet *dataset = new RooDataSet("electrons","electrons",varset,Import(*tree),Cut(thecut));
    cout << "@@@@ dataset has " << dataset->numEntries() << "  entries @@@@" << endl;
    
    model.fitTo(*dataset,SumW2Error(1),Range(-0.1,0.1),Strategy(2),NumCPU(8));

    stringstream resfile;
    resfile << "res_pt_" << ptbins[i] << "To" << ptbins[i+1] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    RooPlot* xframe = x.frame(Title("Energy resolution")) ;
    dataset->plotOn(xframe,DataError(RooAbsData::SumW2) );
    model.plotOn(xframe);
    //    model.paramOn(xframe);
    xframe->Draw(); gPad->Update(); 
    c1->SaveAs(resfile.str().c_str());

    float peak = mean1.getVal();
    float peakerr = mean1.getError();
    float sigma = sigma1.getVal();
    float sigmaerr = sigma1.getError();

    ptM->SetBinContent(i+1,mean);    
    ptM->SetBinError(i+1,meanerr);    
    ptP->SetBinContent(i+1,peak);    
    ptP->SetBinError(i+1,peakerr);    
    ptRMS->SetBinContent(i+1,rms);    
    ptRMS->SetBinError(i+1,rmserr);    
    ptSigma->SetBinContent(i+1,sigma);    
    ptSigma->SetBinError(i+1,sigmaerr);    

    delete dataset;

  }
  
  ptM->Draw();  c1->SaveAs("ptM.png");
  ptP->Draw();  c1->SaveAs("ptP.png");
  ptRMS->Draw();  c1->SaveAs("ptRMS.png");
  ptSigma->Draw();  c1->SaveAs("ptSigma.png");
  



  /*
  // ============ ETA ==============
  cout << "===> RUNNING VS ETA " << endl;
  etaM->SetMaximum(0.1);
  etaP->SetMaximum(0.1);
  etaRMS->SetMaximum(0.1);
  etaSigma->SetMaximum(0.1);

  for(int i=0;i<6;++i) {
    stringstream cut;
    cut << "abs(eta)>" << etabins[i] << "&& abs(eta)<" << etabins[i+1] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_eta_" << etabins[i] << "To" << etabins[i+1] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.2,0.2);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    etaM->SetBinContent(i+1,mean);    
    etaM->SetBinError(i+1,meanerr);    
    etaP->SetBinContent(i+1,peak);    
    etaP->SetBinError(i+1,peakerr);    
    etaRMS->SetBinContent(i+1,rms);    
    etaRMS->SetBinError(i+1,rmserr);    
    etaSigma->SetBinContent(i+1,sigma);    
    etaSigma->SetBinError(i+1,sigmaerr);    

  }
  
  etaM->Draw();  c1->SaveAs("etaM.png");
  etaP->Draw();  c1->SaveAs("etaP.png");
  etaRMS->Draw();  c1->SaveAs("etaRMS.png");
  etaSigma->Draw();  c1->SaveAs("etaSigma.png");



  // ============ VTX ==============
  cout << "===> RUNNING VS NVTX " << endl;
  vtxM->SetMaximum(0.05);
  vtxP->SetMaximum(0.05);
  vtxM->SetMinimum(-0.05);
  vtxP->SetMinimum(-0.05);
  vtxRMS->SetMaximum(0.1);
  vtxSigma->SetMaximum(0.05);

  for(int i=0;i<10;++i) {
    stringstream cut;
    cut << "rho>" << vtxbins[i] << "&& rho<" << vtxbins[i+1] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_vtx_" << vtxbins[i] << "To" << vtxbins[i+1] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.2,0.2);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    vtxM->SetBinContent(i+1,mean);    
    vtxM->SetBinError(i+1,meanerr);    
    vtxP->SetBinContent(i+1,peak);    
    vtxP->SetBinError(i+1,peakerr);    
    vtxRMS->SetBinContent(i+1,rms);    
    vtxRMS->SetBinError(i+1,rmserr);    
    vtxSigma->SetBinContent(i+1,sigma);    
    vtxSigma->SetBinError(i+1,sigmaerr);    

  }
  
  vtxM->Draw();  c1->SaveAs("vtxM.png");
  vtxP->Draw();  c1->SaveAs("vtxP.png");
  vtxRMS->Draw();  c1->SaveAs("vtxRMS.png");
  vtxSigma->Draw();  c1->SaveAs("vtxSigma.png");





  // ============ ELE CLASS ==============
  cout << "===> RUNNING VS ELE CLASS " << endl;
  classM->SetMaximum(0.5);
  classP->SetMaximum(0.05);
  classRMS->SetMaximum(0.2);
  classSigma->SetMaximum(0.05);

  for(int i=0;i<4;++i) {
    stringstream cut;
    cut << "classification==" << classificationbins[i] << "&& abs(p-genp)/genp<0.5 && genp>0";

    stringstream resfile;
    resfile << "res_class_" << classificationbins[i] << ".png";

    tree->Project("EoEt","(p-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    func->SetParLimits(2,0.005,0.05);
    func->SetParLimits(3,0.005,0.05);
    func->SetParLimits(4,0,0.6);
    func->SetParLimits(5,0,0.6);
    func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

    EoEt->Fit("cruijff","","same",-0.2,0.2);

    c1->SaveAs(resfile.str().c_str());

    float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
    float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
    float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
    float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

    classM->SetBinContent(i+1,mean);    
    classM->SetBinError(i+1,meanerr);    
    classP->SetBinContent(i+1,peak);    
    classP->SetBinError(i+1,peakerr);    
    classRMS->SetBinContent(i+1,rms);    
    classRMS->SetBinError(i+1,rmserr);    
    classSigma->SetBinContent(i+1,sigma);    
    classSigma->SetBinError(i+1,sigmaerr);    

  }
  
  classM->Draw();  c1->SaveAs("classM.png");
  classP->Draw();  c1->SaveAs("classP.png");
  classRMS->Draw();  c1->SaveAs("classRMS.png");
  classSigma->Draw();  c1->SaveAs("classSigma.png");

  */

  
  TFile *resultfile = TFile::Open("results_elereg_ggH125.root","recreate");
  ptM->Write(); ptP->Write(); ptRMS->Write(); ptSigma->Write();
  etaM->Write(); etaP->Write(); etaRMS->Write(); etaSigma->Write();
  vtxM->Write(); vtxP->Write(); vtxRMS->Write(); vtxSigma->Write();
  classM->Write(); classP->Write(); classRMS->Write(); classSigma->Write();
  resultfile->Close();


}


void compareResults() {

 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetOptStat("");
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  // comparison on ECAL energy
  TFile *fileScRaw = TFile::Open("elereg/plots/NoRegr/ErawoGenP/results_elereg_ggH125.root");
  TFile *fileSc = TFile::Open("elereg/plots/NoRegr/EoGenP/results_elereg_ggH125.root");
  TFile *fileRegr = TFile::Open("elereg/plots/Regr1/EoGenP/results_elereg_ggH125.root");

  // comparison on P
  // TFile *fileScRaw = TFile::Open("elereg/plots/Regr1/EoGenP/results_elereg_ggH125.root");
  // TFile *fileSc = TFile::Open("elereg/plots/NoRegr/PoGenP/results_elereg_ggH125.root");
  // TFile *fileRegr = TFile::Open("elereg/plots/Regr1/PoGenP/results_elereg_ggH125.root");

  vector<string> histos;
  histos.push_back("ptM");	
  histos.push_back("ptP");	
  histos.push_back("ptRMS");	
  histos.push_back("ptSigma");	
  histos.push_back("etaM");	
  histos.push_back("etaP");	
  histos.push_back("etaRMS");	
  histos.push_back("etaSigma");	
  histos.push_back("vtxM");	
  histos.push_back("vtxP");	
  histos.push_back("vtxRMS");	
  histos.push_back("vtxSigma");	
  histos.push_back("classM");	
  histos.push_back("classP");	
  histos.push_back("classRMS");	
  histos.push_back("classSigma");	


  vector<TFile*> files;
  files.push_back(fileScRaw);
  files.push_back(fileSc);
  files.push_back(fileRegr);

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  for(int h=0;h<(int)histos.size();++h) {
    cout << "Plotting now " << histos[h] << endl;
    c1->Clear();
    TLegend* legend = new TLegend(0.24, 0.70, 0.47, 0.85);
    
    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (0.05);

    Double_t min=100;
    Double_t max=-100;
    for(int i=0;i<3;++i) {
      TH1F *histo = (TH1F*)files[i]->Get(histos[h].c_str());
      min=TMath::Min(min,histo->GetMinimum());
      max=TMath::Max(max,histo->GetMaximum());
    }

    bool reso=(histos[h].find("RMS")!=string::npos || histos[h].find("Sigma")!=string::npos);
    if(reso) c1->Divide(1,2);

    for(int i=0;i<3;++i) {

      TH1F *histo = (TH1F*)files[i]->Get(histos[h].c_str());
      histo->SetLineColor(i+1);
      histo->SetMarkerColor(i+1);
      histo->SetMaximum(max);
      histo->SetMinimum(min);
      if(reso) { 
	histo->SetMinimum(0); 
	c1->cd(1);
      }
      histo->Draw(i==0 ? "pe" : "samepe");

      // comparison in E
      if(i==0) legend->AddEntry(histo,"raw SC");
      if(i==1) legend->AddEntry(histo,"std SC");
      if(i==2) legend->AddEntry(histo,"reg SC");
      // comparison in P
      // if(i==0) legend->AddEntry(histo,"reg SC");
      // if(i==1) legend->AddEntry(histo,"E(std)--p comb.");
      // if(i==2) legend->AddEntry(histo,"E(reg)--p comb.");
      legend->Draw();
    }

    // in case of resolution, make ratio plot
    if(reso) {
      TH1F *histo1 = (TH1F*)files[1]->Get(histos[h].c_str());
      TH1F *histo2 = (TH1F*)files[2]->Get(histos[h].c_str());
      TH1F *histoD = (TH1F*)histo1->Clone();
      histoD->Add(histo2,-1.);
      histoD->SetMarkerStyle(24);
      histoD->SetMarkerSize(2);
      histoD->SetMarkerColor(kBlack);
      histoD->SetLineColor(kBlack);

      stringstream ytitle;
      if(histos[h].find("RMS")!=string::npos) ytitle << "RMS(std SC)-RMS(reg SC)";
      if(histos[h].find("Sigma")!=string::npos) ytitle << "#sigma(std SC)-#sigma(reg SC)";
      histoD->GetYaxis()->SetTitle(ytitle.str().c_str());

      c1->cd(2);
      histoD->Draw();
    }

    stringstream namefile;
    namefile << histos[h] << "_comp.pdf";
    c1->SaveAs(namefile.str().c_str());

  }

}
