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


void makeIdCurve(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",100,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-1.0,1.0);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-1.0,1.0);
  
  treeSig->Project("bdthww_sig","bdthww[0]",cutSig);
  treeBkg->Project("bdthww_bkg","bdthww[0]",cutBkg);
  treeSig->Project("bdthzz_sig","bdthzz[3]",cutSig);
  treeBkg->Project("bdthzz_bkg","bdthzz[3]",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("H #rightarrow WW BDT", bdthww_sig);
  roc.addBackgrounds(bdthww_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow ZZ BDT (DanSiV0)", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");

  roc.drawResults(namefile.Data(),1);
}

void make3IdCurves(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",100,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-1.0,1.0);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-1.0,1.0);
  TH1F *bdthzznoip_sig = new TH1F("bdthzznoip_sig","",100,-1.0,1.0);
  TH1F *bdthzznoip_bkg = new TH1F("bdthzznoip_bkg","",100,-1.0,1.0);

  cout << "cut sig = " << cutSig << " bkg = " << cutBkg << endl;
  
  treeSig->Project("bdthww_sig","bdthww[0]",cutSig);
  treeBkg->Project("bdthww_bkg","bdthww[0]",cutBkg);
  treeSig->Project("bdthzz_sig","bdthzz[3]",cutSig);
  treeBkg->Project("bdthzz_bkg","bdthzz[3]",cutBkg);
  treeSig->Project("bdthzznoip_sig","newbdthww[3]",cutSig);
  treeBkg->Project("bdthzznoip_bkg","newbdthww[3]",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,0.5);
  roc.addSignal("H #rightarrow WW BDT (2011)", bdthww_sig);
  roc.addBackgrounds(bdthww_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow ZZ BDT (DanSiV0)", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow WW BDT (DanSiV0)", bdthzznoip_sig);
  roc.addBackgrounds(bdthzznoip_bkg);
  roc.setCutDirection(">");

  roc.drawResults(namefile.Data(),1);
}


void makeIsolationCurve(TTree *treeSig, TTree* treeBkg, 
			TString cutSig, TString cutBkg, TString namefile) {

  TH1F *isohww_sig = new TH1F("isohww_sig","",100,0.0,2.0);
  TH1F *isohww_bkg = new TH1F("isohww_bkg","",100,0.0,2.0);
  TH1F *isohzz_sig = new TH1F("isohzz_sig","",100,-1.0,2.0);
  TH1F *isohzz_bkg = new TH1F("isohzz_bkg","",100,-1.0,2.0);
  
  treeSig->Project("isohww_sig","combPFIsoHWW/pt",cutSig);
  treeBkg->Project("isohww_bkg","combPFIsoHWW/pt",cutBkg);
  treeSig->Project("isohzz_sig","combPFIsoHZZ/pt",cutSig);
  treeBkg->Project("isohzz_bkg","combPFIsoHZZ/pt",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("H #rightarrow WW iso", isohww_sig);
  roc.addBackgrounds(isohww_bkg);
  roc.setCutDirection("<");
  roc.addSignal("H #rightarrow ZZ iso", isohzz_sig);
  roc.addBackgrounds(isohzz_bkg);
  roc.setCutDirection("<");

  roc.drawResults(namefile.Data(),1);
}

void make3IsolationCurves(TTree *treeSig, TTree* treeBkg, 
			  TString cutSig, TString cutBkg, TString namefile) {

  TH1F *isohww_sig = new TH1F("isohww_sig","",100,0.0,2.0);
  TH1F *isohww_bkg = new TH1F("isohww_bkg","",100,0.0,2.0);
  TH1F *isohzz_sig = new TH1F("isohzz_sig","",100,-1.0,2.0);
  TH1F *isohzz_bkg = new TH1F("isohzz_bkg","",100,-1.0,2.0);
  TH1F *isohzznoEA_sig = new TH1F("isohzznoEA_sig","",100,-1.0,2.0);
  TH1F *isohzznoEA_bkg = new TH1F("isohzznoEA_bkg","",100,-1.0,2.0);
  
  treeSig->Project("isohww_sig","combPFIsoHWW/pt",cutSig);
  treeBkg->Project("isohww_bkg","combPFIsoHWW/pt",cutBkg);
  treeSig->Project("isohzz_sig","combPFIsoHZZ/pt",cutSig);
  treeBkg->Project("isohzz_bkg","combPFIsoHZZ/pt",cutBkg);
  treeSig->Project("isohzznoEA_sig","combPFIsoHZZNoEA/pt",cutSig);
  treeBkg->Project("isohzznoEA_bkg","combPFIsoHZZNoEA/pt",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("PF iso", isohww_sig);
  roc.addBackgrounds(isohww_bkg);
  roc.setCutDirection("<");
  roc.addSignal("PF iso, new vetoes", isohzznoEA_sig);
  roc.addBackgrounds(isohzznoEA_bkg);
  roc.setCutDirection("<");
  roc.addSignal("PF iso, new vetoes, EA", isohzz_sig);
  roc.addBackgrounds(isohzz_bkg);
  roc.setCutDirection("<");
  roc.drawResults(namefile.Data(),1);
}


int main(int argc, char* argv[]) {

  TFile *fileSig, *fileBkg;
  TTree *treeSig, *treeBkg;
  fileSig = fileBkg = 0;
  treeSig = treeBkg = 0;

  fileSig = TFile::Open("macro/results_data/electrons_zeemc.root");
  fileBkg = TFile::Open("macro/results_data/fakes.root");
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
  // cutBase.push_back(TString("abs(eta)<1.0                   && pt<10"));
  // cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt<10"));
  // cutBase.push_back(TString("abs(eta)>1.479                 && pt<10"));
  cutBase.push_back(TString("abs(eta)<1.0                   && pt>10 && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt>10 && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt>10 && pt<20"));
  cutBase.push_back(TString("abs(eta)<1.0                   && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt>20"));

  vector<TString> cutSignal;
  for(int i=0;i<(int)cutBase.size();++i) 
    cutSignal.push_back(cutBase[i]+TString("&& DenomFakeSmurf && mcmatch && pt<35"));

  vector<TString> cutBackground;
  for(int i=0;i<(int)cutBase.size();++i)
    cutBackground.push_back(cutBase[i]+TString("&& DenomFakeSmurf && pt<35 && !(run<=173692 && event%2==0)")); // HWW used the even events to train

  vector<TString> id;
  // id.push_back(TString("ROC_IdOnly_Data_inEB_VeryLowPt.pdf"));
  // id.push_back(TString("ROC_IdOnly_Data_outEB_VeryLowPt.pdf"));
  // id.push_back(TString("ROC_IdOnly_Data_EE_VeryLowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_inEB_LowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_outEB_LowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_EE_LowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_inEB_HighPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_outEB_HighPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_EE_HighPt.pdf"));

  for(int i=0;i<(int)cutBase.size();++i) {
    // makeIdCurve(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
    make3IdCurves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
  }

  // HZZ isolations are in friend trees
  treeSig->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/electrons_zeemc_hzzisoFriend.root");
  treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/fakes_hzzisoFriend.root");

  for(int i=0;i<(int)cutBase.size();++i) {
    id[i].ReplaceAll("IdOnly","IsoOnly");
    make3IsolationCurves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
  }

  // and now in 2 big bins in nvtx
  vector<TString> cutSignalHiPU, cutSignalLoPU, cutBackgroundHiPU, cutBackgroundLoPU;

  for(int i=0;i<(int)cutBase.size();++i) {
    cutSignalLoPU.push_back(cutSignal[i]+TString(" && vertices<8"));
    cutBackgroundLoPU.push_back(cutBackground[i]+TString(" && vertices<8"));
    id[i].ReplaceAll("IsoOnly","IsoOnlyVtx1To7");
    make3IsolationCurves(treeSig,treeBkg,cutSignalLoPU[i],cutBackgroundLoPU[i],id[i]);
    id[i].ReplaceAll("IsoOnlyVtx1To7","IdOnlyVtx1To7");
    make3IdCurves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
  }

  for(int i=0;i<(int)cutBase.size();++i) {
    cutSignalHiPU.push_back(cutSignal[i]+TString(" && vertices>=8"));
    cutBackgroundHiPU.push_back(cutBackground[i]+TString(" && vertices>=8"));
    id[i].ReplaceAll("IdOnlyVtx1To7","IsoOnlyVtx8ToInf");
    make3IsolationCurves(treeSig,treeBkg,cutSignalHiPU[i],cutBackgroundHiPU[i],id[i]);
    id[i].ReplaceAll("IsoOnlyVtx8ToInf","IdOnlyVtx8ToInf");
    make3IdCurves(treeSig,treeBkg,cutSignalHiPU[i],cutBackgroundHiPU[i],id[i]);
  }


  return 0;
}
