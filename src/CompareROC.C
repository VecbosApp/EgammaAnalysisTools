// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>

// Offline analysis includes
#include "CommonTools/include/FiguresOfMeritEvaluator.h"

#include <string>
#include <sstream>
#include <fstream>

using namespace std;

void makeCurve(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile);
void make3Curves(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile);

int main(int argc, char* argv[]) {

  TFile *fileSig, *fileBkg;
  TTree *treeSig, *treeBkg;
  fileSig = fileBkg = 0;
  treeSig = treeBkg = 0;

  fileSig = TFile::Open("macro/results_data/merged.root");
  fileBkg = TFile::Open("macro/results_data_fakes/merged.root");
  if( fileSig && fileBkg) {
    fileSig->cd();
    treeSig = (TTree*)fileSig->Get("eleIDdir/T1");
    fileBkg->cd();
    treeBkg = (TTree*)fileBkg->Get("eleIDdir/T1");
  } else {
    cout << "File " << fileSig << " or " << fileBkg
	 << " not existing !" << endl;
    return 0;
  }

  if(!treeSig || !treeBkg) {
    cout << "Tree eleIDdir/T1 not existing inside signal or background files!" << endl;
    return 0;
  }

  cout << "All files OK!" << endl;

  vector<TString> cutBase;
  cutBase.push_back(TString("abs(eta)<1.0                   && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt<20"));
  cutBase.push_back(TString("abs(eta)<1.0                   && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt>20"));

  vector<TString> cutSignal;
  for(int i=0;i<(int)cutBase.size();++i) 
    cutSignal.push_back(cutBase[i]+TString("&& DenomFakeSmurf && abs(mass-91.1876)<7.5"));

  vector<TString> cutBackground;
  for(int i=0;i<(int)cutBase.size();++i)
    cutBackground.push_back(cutBase[i]+TString("&& DenomFakeSmurf"));

  vector<TString> id;
  id.push_back(TString("ROC_Data_inEB_LowPt.pdf"));
  id.push_back(TString("ROC_Data_outEB_LowPt.pdf"));
  id.push_back(TString("ROC_Data_EE_LowPt.pdf"));
  id.push_back(TString("ROC_Data_inEB_HighPt.pdf"));
  id.push_back(TString("ROC_Data_outEB_HighPt.pdf"));
  id.push_back(TString("ROC_Data_EE_HighPt.pdf"));

  for(int i=0;i<(int)cutBase.size();++i) 
    makeCurve(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
  // make3Curves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);

  return 0;
}

void makeCurve(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",100,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-0.4,0.3);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-0.4,0.3);
  
  treeSig->Project("bdthww_sig","bdthww",cutSig);
  treeBkg->Project("bdthww_bkg","bdthww",cutBkg);
  treeSig->Project("bdthzz_sig","bdthzz",cutSig);
  treeBkg->Project("bdthzz_bkg","bdthzz",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.addSignal("H #rightarrow WW BDT", bdthww_sig);
  roc.addBackgrounds(bdthww_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow ZZ BDT", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");

  roc.drawResults(namefile.Data());
}

void make3Curves(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",100,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-0.4,0.3);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-0.4,0.3);
  TH1F *bdthzznoip_sig = new TH1F("bdthzznoip_sig","",100,-0.4,0.3);
  TH1F *bdthzznoip_bkg = new TH1F("bdthzznoip_bkg","",100,-0.4,0.3);

  cout << "cut sig = " << cutSig << " bkg = " << cutBkg << endl;
  
  treeSig->Project("bdthww_sig","bdthww",cutSig);
  treeBkg->Project("bdthww_bkg","bdthww",cutBkg);
  treeSig->Project("bdthzz_sig","bdthzz",cutSig);
  treeBkg->Project("bdthzz_bkg","bdthzz",cutBkg);
  treeSig->Project("bdthzznoip_sig","bdthzznoip",cutSig);
  treeBkg->Project("bdthzznoip_bkg","bdthzznoip",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.addSignal("H #rightarrow WW BDT", bdthww_sig);
  roc.addBackgrounds(bdthww_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow ZZ BDT", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow ZZ BDT (no IP)", bdthzznoip_sig);
  roc.addBackgrounds(bdthzznoip_bkg);
  roc.setCutDirection(">");

  roc.drawResults(namefile.Data());
}
