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

  TFile *file = TFile::Open("/Users/emanuele/Work/data/hzz4l/electronreg/HZZ4L_53X_S1_V10_S2_V01/MC/EleRegr1/ggH125_newtraining_nosmear.root");
  // TFile *file = TFile::Open("/Users/emanuele/Work/data/hzz4l/electronreg/HZZ4L_42X_S1_V07_S2_V00/MC/EleRegr1/ggH125.root");
  TTree *tree = (TTree*)file->Get("electronTree/probe_tree");

  TH1F *EoEt = new TH1F("EoEt","",151,-0.5,0.5);
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

  TH1F* ptMEB     = (TH1F*)ptM->Clone("ptMEB");          TH1F* ptMEE     = (TH1F*)ptM->Clone("ptMEE"); 
  TH1F* ptPEB     = (TH1F*)ptP->Clone("ptPEB");          TH1F* ptPEE     = (TH1F*)ptP->Clone("ptPEE"); 
  TH1F* ptRMSEB   = (TH1F*)ptRMS->Clone("ptRMSEB");      TH1F* ptRMSEE   = (TH1F*)ptRMS->Clone("ptRMSEE"); 
  TH1F* ptSigmaEB = (TH1F*)ptSigma->Clone("ptSigmaEB");  TH1F* ptSigmaEE = (TH1F*)ptSigma->Clone("ptSigmaEE"); 
  TH1F* vtxMEB     = (TH1F*)vtxM->Clone("vtxMEB");          TH1F* vtxMEE     = (TH1F*)vtxM->Clone("vtxMEE"); 
  TH1F* vtxPEB     = (TH1F*)vtxP->Clone("vtxPEB");          TH1F* vtxPEE     = (TH1F*)vtxP->Clone("vtxPEE"); 
  TH1F* vtxRMSEB   = (TH1F*)vtxRMS->Clone("vtxRMSEB");      TH1F* vtxRMSEE   = (TH1F*)vtxRMS->Clone("vtxRMSEE"); 
  TH1F* vtxSigmaEB = (TH1F*)vtxSigma->Clone("vtxSigmaEB");  TH1F* vtxSigmaEE = (TH1F*)vtxSigma->Clone("vtxSigmaEE"); 
  TH1F* classMEB     = (TH1F*)classM->Clone("classMEB");          TH1F* classMEE     = (TH1F*)classM->Clone("classMEE"); 
  TH1F* classPEB     = (TH1F*)classP->Clone("classPEB");          TH1F* classPEE     = (TH1F*)classP->Clone("classPEE"); 
  TH1F* classRMSEB   = (TH1F*)classRMS->Clone("classRMSEB");      TH1F* classRMSEE   = (TH1F*)classRMS->Clone("classRMSEE"); 
  TH1F* classSigmaEB = (TH1F*)classSigma->Clone("classSigmaEB");  TH1F* classSigmaEE = (TH1F*)classSigma->Clone("classSigmaEE"); 

  for(int e=0;e<2;++e) {

    stringstream etacut;
    etacut << ((e==0) ? "abs(eta)<1.479" : "abs(eta)>1.479");

    stringstream etabin;
    etabin << ((e==0) ? "_EB_" : "_EE_");
    
    // ============ P ==============
    cout << "===> RUNNING VS P " << endl;
    for(int i=0;i<11;++i) {
      stringstream cut;
      cut << "pt>" << ptbins[i] << "&& pt<" << ptbins[i+1] << "&& abs((ecalE-genp)/genp)<0.5 && genp>0 && pt>7 && " << etacut.str();
      cout << cut.str() << endl;

      stringstream resfile;
      resfile << "res_pt" << etabin.str() << ptbins[i] << "To" << ptbins[i+1] << ".png";

      tree->Project("EoEt","(ecalE-genp)/genp",cut.str().c_str());
      EoEt->Draw();

      float mean = EoEt->GetMean();
      float meanerr = EoEt->GetMeanError();
      float rms = EoEt->GetRMS();
      float rmserr = EoEt->GetRMSError();
    
      // fit the Gaussian core
      func->SetParameter(1,EoEt->GetMean());
      func->SetParameter(2,EoEt->GetRMS());
      func->SetParameter(3,EoEt->GetRMS());
      if(e==0) func->SetParLimits(1,-0.05,0.05);
      if(e==1) func->SetParLimits(1,-0.08,0.05);
      func->SetParLimits(2,0.007,0.2);
      func->SetParLimits(3,0.007,0.2);
      func->SetParLimits(4,0.01,0.4);
      func->SetParLimits(5,0.01,0.4);
      func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

      EoEt->Fit("cruijff","","same",-0.08,0.08);

      c1->SaveAs(resfile.str().c_str());

      float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
      float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
      float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
      float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

      if(e==0) {
	ptMEB->SetBinContent(i+1,mean);    
	ptMEB->SetBinError(i+1,meanerr);    
	ptPEB->SetBinContent(i+1,peak);    
	ptPEB->SetBinError(i+1,peakerr);    
	ptRMSEB->SetBinContent(i+1,rms);    
	ptRMSEB->SetBinError(i+1,rmserr);    
	ptSigmaEB->SetBinContent(i+1,sigma);    
	ptSigmaEB->SetBinError(i+1,sigmaerr);    
      } else { 
	ptMEE->SetBinContent(i+1,mean);    
	ptMEE->SetBinError(i+1,meanerr);    
	ptPEE->SetBinContent(i+1,peak);    
	ptPEE->SetBinError(i+1,peakerr);    
	ptRMSEE->SetBinContent(i+1,rms);    
	ptRMSEE->SetBinError(i+1,rmserr);    
	ptSigmaEE->SetBinContent(i+1,sigma);    
	ptSigmaEE->SetBinError(i+1,sigmaerr);    
      }
    }
  
    
    
    
    
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
      cut << "rho>" << vtxbins[i] << "&& rho<" << vtxbins[i+1] << "&& abs((ecalE-genp)/genp)<0.5 && genp>0 && pt>7 && " << etacut.str();
      
      stringstream resfile;
      resfile << "res_vtx" << etabin.str() << vtxbins[i] << "To" << vtxbins[i+1] << ".png";
      
      tree->Project("EoEt","(ecalE-genp)/genp",cut.str().c_str());
      EoEt->Draw();
      
      float mean = EoEt->GetMean();
      float meanerr = EoEt->GetMeanError();
      float rms = EoEt->GetRMS();
      float rmserr = EoEt->GetRMSError();
      
      // fit the Gaussian core
      func->SetParameter(1,EoEt->GetMean());
      func->SetParameter(2,EoEt->GetRMS());
      func->SetParameter(3,EoEt->GetRMS());
      if(e==0) func->SetParLimits(1,-0.05,0.05);
      if(e==1) func->SetParLimits(1,-0.08,0.05);
      func->SetParLimits(2,0.007,0.2);
      func->SetParLimits(3,0.007,0.2);
      func->SetParLimits(4,0.01,0.4);
      func->SetParLimits(5,0.01,0.4);
      func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 
      
      EoEt->Fit("cruijff","","same",-0.08,0.08);
      
      c1->SaveAs(resfile.str().c_str());
      
      float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
      float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
      float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
      float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;
     
      if(e==0) {
	vtxMEB->SetBinContent(i+1,mean);    
	vtxMEB->SetBinError(i+1,meanerr);    
	vtxPEB->SetBinContent(i+1,peak);    
	vtxPEB->SetBinError(i+1,peakerr);    
	vtxRMSEB->SetBinContent(i+1,rms);    
	vtxRMSEB->SetBinError(i+1,rmserr);    
	vtxSigmaEB->SetBinContent(i+1,sigma);    
	vtxSigmaEB->SetBinError(i+1,sigmaerr);    
      } else {
	vtxMEE->SetBinContent(i+1,mean);    
	vtxMEE->SetBinError(i+1,meanerr);    
	vtxPEE->SetBinContent(i+1,peak);    
	vtxPEE->SetBinError(i+1,peakerr);    
	vtxRMSEE->SetBinContent(i+1,rms);    
	vtxRMSEE->SetBinError(i+1,rmserr);    
	vtxSigmaEE->SetBinContent(i+1,sigma);    
	vtxSigmaEE->SetBinError(i+1,sigmaerr);    
      }
    }



    // ============ ELE CLASS ==============
    cout << "===> RUNNING VS ELE CLASS " << endl;
    classM->SetMaximum(0.5);
    classP->SetMaximum(0.05);
    classRMS->SetMaximum(0.2);
    classSigma->SetMaximum(0.05);

    for(int i=0;i<4;++i) {
      stringstream cut;
      cut << "classification==" << classificationbins[i] << "&& abs((ecalE-genp)/genp)<0.5 && genp>0 && pt>7 && " << etacut.str();

      stringstream resfile;
      resfile << "res_class" << etabin.str() << classificationbins[i] << ".png";

      tree->Project("EoEt","(ecalE-genp)/genp",cut.str().c_str());
      EoEt->Draw();

      float mean = EoEt->GetMean();
      float meanerr = EoEt->GetMeanError();
      float rms = EoEt->GetRMS();
      float rmserr = EoEt->GetRMSError();
    
      // fit the Gaussian core
      func->SetParameter(1,EoEt->GetMean());
      func->SetParameter(2,EoEt->GetRMS());
      func->SetParameter(3,EoEt->GetRMS());
      if(e==0) func->SetParLimits(1,-0.05,0.05);
      if(e==1) func->SetParLimits(1,-0.08,0.05);
      func->SetParLimits(2,0.007,0.2);
      func->SetParLimits(3,0.007,0.2);
      func->SetParLimits(4,0.01,0.4);
      func->SetParLimits(5,0.01,0.4);
      func->SetParNames ("Constant","Mean","sigmaL","sigmaR","alphaL","alphaR"); 

      EoEt->Fit("cruijff","","same",-0.08,0.08);

      c1->SaveAs(resfile.str().c_str());

      float peak = EoEt->GetFunction("cruijff")->GetParameter(1);
      float peakerr = EoEt->GetFunction("cruijff")->GetParError(1);
      float sigma = (EoEt->GetFunction("cruijff")->GetParameter(2) +  EoEt->GetFunction("cruijff")->GetParameter(3))/2.;
      float sigmaerr = (EoEt->GetFunction("cruijff")->GetParError(2) + EoEt->GetFunction("cruijff")->GetParError(3))/2.;

      if(e==0) {
	classMEB->SetBinContent(i+1,mean);    
	classMEB->SetBinError(i+1,meanerr);    
	classPEB->SetBinContent(i+1,peak);    
	classPEB->SetBinError(i+1,peakerr);    
	classRMSEB->SetBinContent(i+1,rms);    
	classRMSEB->SetBinError(i+1,rmserr);    
	classSigmaEB->SetBinContent(i+1,sigma);    
	classSigmaEB->SetBinError(i+1,sigmaerr);    
      } else {
	classMEE->SetBinContent(i+1,mean);    
	classMEE->SetBinError(i+1,meanerr);    
	classPEE->SetBinContent(i+1,peak);    
	classPEE->SetBinError(i+1,peakerr);    
	classRMSEE->SetBinContent(i+1,rms);    
	classRMSEE->SetBinError(i+1,rmserr);    
	classSigmaEE->SetBinContent(i+1,sigma);    
	classSigmaEE->SetBinError(i+1,sigmaerr);    
      }
    }

  }


  ptMEB->Draw();  c1->SaveAs("ptM_EB.png");         ptMEE->Draw();  c1->SaveAs("ptM_EE.png");
  ptPEB->Draw();  c1->SaveAs("ptP_EB.png");         ptPEE->Draw();  c1->SaveAs("ptP_EE.png");
  ptRMSEB->Draw();  c1->SaveAs("ptRMS_EB.png");     ptRMSEE->Draw();  c1->SaveAs("ptRMS_EE.png");
  ptSigmaEB->Draw();  c1->SaveAs("ptSigma_EB.png"); ptSigmaEE->Draw();  c1->SaveAs("ptSigma_EE.png");

  vtxMEB->Draw();  c1->SaveAs("vtxM_EB.png");         vtxMEE->Draw();  c1->SaveAs("vtxM_EE.png");
  vtxPEB->Draw();  c1->SaveAs("vtxP_EB.png");         vtxPEE->Draw();  c1->SaveAs("vtxP_EE.png");
  vtxRMSEB->Draw();  c1->SaveAs("vtxRMS_EB.png");     vtxRMSEE->Draw();  c1->SaveAs("vtxRMS_EE.png");
  vtxSigmaEB->Draw();  c1->SaveAs("vtxSigma_EB.png"); vtxSigmaEE->Draw();  c1->SaveAs("vtxSigma_EE.png");

  classMEB->Draw();  c1->SaveAs("classM_EB.png");         classMEE->Draw();  c1->SaveAs("classM_EE.png");
  classPEB->Draw();  c1->SaveAs("classP_EB.png");         classPEE->Draw();  c1->SaveAs("classP_EE.png");
  classRMSEB->Draw();  c1->SaveAs("classRMS_EB.png");     classRMSEE->Draw();  c1->SaveAs("classRMS_EE.png");
  classSigmaEB->Draw();  c1->SaveAs("classSigma_EB.png"); classSigmaEE->Draw();  c1->SaveAs("classSigma_EE.png");


  // ============ ETA ==============
  cout << "===> RUNNING VS ETA " << endl;
  etaM->SetMaximum(0.1);
  etaP->SetMaximum(0.1);
  etaRMS->SetMaximum(0.1);
  etaSigma->SetMaximum(0.1);

  for(int i=0;i<6;++i) {
    stringstream cut;
    cut << "abs(eta)>" << etabins[i] << "&& abs(eta)<" << etabins[i+1] << "&& abs((ecalE-genp)/genp)<0.5 && genp>0 && pt>7";

    stringstream resfile;
    resfile << "res_eta_" << etabins[i] << "To" << etabins[i+1] << ".png";

    tree->Project("EoEt","(ecalE-genp)/genp",cut.str().c_str());
    EoEt->Draw();

    float mean = EoEt->GetMean();
    float meanerr = EoEt->GetMeanError();
    float rms = EoEt->GetRMS();
    float rmserr = EoEt->GetRMSError();
    
    // fit the Gaussian core
    func->SetParameter(1,EoEt->GetMean());
    func->SetParameter(2,EoEt->GetRMS());
    func->SetParameter(3,EoEt->GetRMS());
    if(i<4) func->SetParLimits(1,-0.05,0.05);
    else func->SetParLimits(1,-0.08,0.05);
    func->SetParLimits(2,0.007,0.2);
    func->SetParLimits(3,0.01,0.2);
    func->SetParLimits(4,0.01,0.4);
    func->SetParLimits(5,0.01,0.4);
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


  
  TFile *resultfile = TFile::Open("results_elereg_ggH125.root","recreate");
  ptMEB->Write(); ptPEB->Write(); ptRMSEB->Write(); ptSigmaEB->Write();
  ptMEE->Write(); ptPEE->Write(); ptRMSEE->Write(); ptSigmaEE->Write();
  vtxMEB->Write(); vtxPEB->Write(); vtxRMSEB->Write(); vtxSigmaEB->Write();
  vtxMEE->Write(); vtxPEE->Write(); vtxRMSEE->Write(); vtxSigmaEE->Write();
  classMEB->Write(); classPEB->Write(); classRMSEB->Write(); classSigmaEB->Write();
  classMEE->Write(); classPEE->Write(); classRMSEE->Write(); classSigmaEE->Write();
  etaM->Write(); etaP->Write(); etaRMS->Write(); etaSigma->Write();
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
    for(int e=0;e<2;++e) {
      stringstream ecal;
      if(!(histos[h].find("eta")!=string::npos)) ecal << ((e==0) ? "EB" : "EE");
      cout << "Plotting now " << histos[h] << ecal.str() << endl;
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
	TH1F *histo = (TH1F*)files[i]->Get((histos[h]+ecal.str()).c_str());
	min=TMath::Min(min,histo->GetMinimum()-1.5 * histo->GetBinError(1));
	max=TMath::Max(max,histo->GetMaximum()+1.5 * histo->GetBinError(1));
      }
      
      bool reso=(histos[h].find("RMS")!=string::npos || histos[h].find("Sigma")!=string::npos);
      if(reso) c1->Divide(1,2);
      
      for(int i=0;i<3;++i) {
	
	TH1F *histo = (TH1F*)files[i]->Get((histos[h]+ecal.str()).c_str());
	histo->SetLineColor(i+1);
	histo->SetMarkerColor(i+1);
	cout << histo->GetName() << " min = " << min << "   max = " << max << endl;
	if(reso) { 
	  histo->SetMinimum(0); 
	  c1->cd(1);
	}
	histo->GetYaxis()->SetRangeUser(min,max);
	histo->Draw(i==0 ? "pe" : "samepe");
	histo->GetYaxis()->SetRangeUser(min,max);
	
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
	TH1F *histo1 = (TH1F*)files[1]->Get((histos[h]+ecal.str()).c_str());
	TH1F *histo2 = (TH1F*)files[2]->Get((histos[h]+ecal.str()).c_str());
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
      if(!(histos[h].find("eta")!=string::npos)) namefile << histos[h] << "_" << ecal.str() << "_comp.png";
      else namefile << histos[h] << "_comp.png";
      c1->SaveAs(namefile.str().c_str());
      
    }
  }

}

