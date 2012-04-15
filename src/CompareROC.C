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


void makeIdCurveHWW(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",100,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-1.0,1.0);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-1.0,1.0);
  
  treeSig->Project("bdthww_sig","bdthww[0]",cutSig);
  treeBkg->Project("bdthww_bkg","bdthww[0]",cutBkg);
  treeSig->Project("bdthzz_sig","newbdthww[3]",cutSig);
  treeBkg->Project("bdthzz_bkg","newbdthww[3]",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("H #rightarrow WW BDT (2011)", bdthww_sig);
  roc.addBackgrounds(bdthww_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow WW BDT (2012)", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");

  roc.drawResults(namefile.Data(),1);
}

void makeIdCurveHZZ(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *cichzz_sig = new TH1F("cichzz_sig","",100,-1.0,1.0);
  TH1F *cichzz_bkg = new TH1F("cichzz_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-1.0,1.0);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-1.0,1.0);

  treeSig->Project("cichzz_sig","cicmediumid",cutSig);
  treeBkg->Project("cichzz_bkg","cicmediumid",cutBkg);
  treeSig->Project("bdthzz_sig","newbdthww[3]",cutSig);
  treeBkg->Project("bdthzz_bkg","newbdthww[3]",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("H #rightarrow ZZ CiC (2011)", cichzz_sig);
  roc.addBackgrounds(cichzz_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow ZZ BDT (2012)", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");

  roc.drawResults(namefile.Data(),1);
}

void make3IdCurves(TTree *treeSig, TTree* treeBkg, TString cutSig, TString cutBkg, TString namefile) {

  TH1F *bdthww_sig = new TH1F("bdthww_sig","",100,-1.0,1.0);
  TH1F *bdthww_bkg = new TH1F("bdthww_bkg","",100,-1.0,1.0);
  TH1F *bdthwwiso_sig = new TH1F("bdthwwiso_sig","",100,-1.0,1.0);
  TH1F *bdthwwiso_bkg = new TH1F("bdthwwiso_bkg","",100,-1.0,1.0);
  TH1F *bdthzz_sig = new TH1F("bdthzz_sig","",100,-1.0,1.0);
  TH1F *bdthzz_bkg = new TH1F("bdthzz_bkg","",100,-1.0,1.0);
  TH2F *bdthzziso_sig = new TH2F("bdthzziso_sig","",100,-1.0,1.0,100,0.0,5.0);
  TH2F *bdthzziso_bkg = new TH2F("bdthzziso_bkg","",100,-1.0,1.0,100,0.0,5.0);

  cout << "cut sig = " << cutSig << " bkg = " << cutBkg << endl;
  
  treeSig->Project("bdthww_sig","bdthww[0]",cutSig);
  treeBkg->Project("bdthww_bkg","bdthww[0]",cutBkg);
  treeSig->Project("bdthwwiso_sig","bdthww[1]",cutSig);
  treeBkg->Project("bdthwwiso_bkg","bdthww[1]",cutBkg);
  treeSig->Project("bdthzz_sig","newbdthww[3]",cutSig);
  treeBkg->Project("bdthzz_bkg","newbdthww[3]",cutBkg);
  treeSig->Project("bdthzziso_sig","combPFIsoHZZ/pt:newbdthww[3]",cutSig);
  treeBkg->Project("bdthzziso_bkg","combPFIsoHZZ/pt:newbdthww[3]",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,0.5);
  roc.addSignal("H #rightarrow WW BDT (2011)", bdthww_sig);
  roc.addBackgrounds(bdthww_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow WW BDT (with iso)", bdthwwiso_sig);
  roc.addBackgrounds(bdthwwiso_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow WW BDT (e#gamma)", bdthzz_sig);
  roc.addBackgrounds(bdthzz_bkg);
  roc.setCutDirection(">");
  roc.addSignal("H #rightarrow WW BDT (e#gamma + iso. cut)", bdthzziso_sig);
  roc.addBackgrounds(bdthzziso_bkg);
  roc.setCutDirectionXY(">:<");

  roc.drawResults(namefile.Data(),1);
}


void makeIsolationCurve(TTree *treeSig, TTree* treeBkg, 
			TString cutSig, TString cutBkg, TString namefile) {

  TH1F *isohww_sig = new TH1F("isohww_sig","",100,0.0,2.0);
  TH1F *isohww_bkg = new TH1F("isohww_bkg","",100,0.0,2.0);
  TH1F *isohzz_sig = new TH1F("isohzz_sig","",100,-1.0,2.0);
  TH1F *isohzz_bkg = new TH1F("isohzz_bkg","",100,-1.0,2.0);
  
  treeSig->Project("isohww_sig","chaPFIso[3]/pt",cutSig);
  treeBkg->Project("isohww_bkg","chaPFIso[3]/pt",cutBkg);
  treeSig->Project("isohzz_sig","combPFIsoHZZ/pt",cutSig);
  treeBkg->Project("isohzz_bkg","combPFIsoHZZ/pt",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("charged PF iso", isohww_sig);
  roc.addBackgrounds(isohww_bkg);
  roc.setCutDirection("<");
  roc.addSignal("combined PF iso (EA)", isohzz_sig);
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

  //  treeSig->Project("isohww_sig","chaPFIso/pt",cutSig);
  //  treeBkg->Project("isohww_bkg","chaPFIso/pt",cutBkg);
  treeSig->Project("isohzz_sig","combPFIsoHZZ/pt",cutSig);
  treeBkg->Project("isohzz_bkg","combPFIsoHZZ/pt",cutBkg);
  treeSig->Project("isohzznoEA_sig","combPFIsoHZZNoEA/pt",cutSig);
  treeBkg->Project("isohzznoEA_bkg","combPFIsoHZZNoEA/pt",cutBkg);

  FiguresOfMeritEvaluator roc;
  roc.setRange(0.6,1,0,1.0);
  roc.addSignal("PF iso (old vetoes)", isohww_sig);
  roc.addBackgrounds(isohww_bkg);
  roc.setCutDirection("<");
  roc.addSignal("full PF iso", isohzznoEA_sig);
  roc.addBackgrounds(isohzznoEA_bkg);
  roc.setCutDirection("<");
  roc.addSignal("full. PF iso, EA", isohzz_sig);
  roc.addBackgrounds(isohzz_bkg);
  roc.setCutDirection("<");
  roc.drawResults(namefile.Data(),1);
}


int main(int argc, char* argv[]) {

  // switch to decide if run on triggering or non-triggering eles
  bool trg = true;

  TFile *fileSig, *fileBkg;
  TTree *treeSig, *treeBkg;
  fileSig = fileBkg = 0;
  treeSig = treeBkg = 0;

  if(trg) {
    fileSig = TFile::Open("macro/results_data/electrons_zeemc.root");
    fileBkg = TFile::Open("macro/results_data/fakes.root");
  } else {
    fileSig = TFile::Open("macro/results_data/electrons_zeemc.root");
    fileBkg = TFile::Open("macro/results_data/fakes-unbiased-wlnu.root");
  }
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
  if(!trg) {
    cutBase.push_back(TString("abs(eta)<1.0                   && pt<10"));
    cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt<10"));
    cutBase.push_back(TString("abs(eta)>1.479                 && pt<10"));
  }
  cutBase.push_back(TString("abs(eta)<1.0                   && pt>10 && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt>10 && pt<20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt>10 && pt<20"));
  cutBase.push_back(TString("abs(eta)<1.0                   && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.0 && abs(eta)<1.479 && pt>20"));
  cutBase.push_back(TString("abs(eta)>1.479                 && pt>20"));

  if(trg) {
    for(int i=0;i<(int)cutBase.size();++i) cutBase[i]=cutBase[i]+TString("&& DenomFakeSmurf"); 
  }

  vector<TString> cutSignal;
  for(int i=0;i<(int)cutBase.size();++i) 
    cutSignal.push_back(cutBase[i]+TString("&& mcmatch && pt<35"));

  vector<TString> cutBackground;
  for(int i=0;i<(int)cutBase.size();++i)
    cutBackground.push_back(cutBase[i]+TString("&& pt<35 && !(run<=173692 && event%2==0)")); // HWW used the even events to train

  vector<TString> id;
  if(!trg) {
    id.push_back(TString("ROC_IdOnly_Data_inEB_VeryLowPt.pdf"));
    id.push_back(TString("ROC_IdOnly_Data_outEB_VeryLowPt.pdf"));
    id.push_back(TString("ROC_IdOnly_Data_EE_VeryLowPt.pdf"));
  }
  id.push_back(TString("ROC_IdOnly_Data_inEB_LowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_outEB_LowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_EE_LowPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_inEB_HighPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_outEB_HighPt.pdf"));
  id.push_back(TString("ROC_IdOnly_Data_EE_HighPt.pdf"));

  // HZZ isolations and new id bits are in friend trees
  treeSig->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/electrons_zeemc_hzzisoFriend.root");
  if(trg) treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/fakes_hzzisoFriend.root");
  else treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/fakes-unbiased-wlnu_hzzisoFriend.root");
  treeSig->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/electrons_zeemc_hzzidbitsFriend.root");
  if(trg) treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/fakes_hzzidbitsFriend.root");
  else treeBkg->AddFriend("eleIDdir/isoT1 = eleIDdir/T1","macro/results_data/fakes-unbiased-wlnu_hzzidbitsFriend.root");

  for(int i=0;i<(int)cutBase.size();++i) {
    if(trg) //makeIdCurveHWW(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
      make3IdCurves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
    else makeIdCurveHZZ(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
  }

  for(int i=0;i<(int)cutBase.size();++i) {
    id[i].ReplaceAll("IdOnly","IsoOnly");
    // makeIsolationCurve(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
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
    if(trg) //makeIdCurveHWW(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
      make3IdCurves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
    else makeIdCurveHZZ(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
  }

  for(int i=0;i<(int)cutBase.size();++i) {
    cutSignalHiPU.push_back(cutSignal[i]+TString(" && vertices>=8"));
    cutBackgroundHiPU.push_back(cutBackground[i]+TString(" && vertices>=8"));
    id[i].ReplaceAll("IdOnlyVtx1To7","IsoOnlyVtx8ToInf");
    make3IsolationCurves(treeSig,treeBkg,cutSignalHiPU[i],cutBackgroundHiPU[i],id[i]);
    id[i].ReplaceAll("IsoOnlyVtx8ToInf","IdOnlyVtx8ToInf");
    if(trg) // makeIdCurveHWW(treeSig,treeBkg,cutSignalHiPU[i],cutBackgroundHiPU[i],id[i]);
      make3IdCurves(treeSig,treeBkg,cutSignal[i],cutBackground[i],id[i]);
    else makeIdCurveHZZ(treeSig,treeBkg,cutSignalHiPU[i],cutBackgroundHiPU[i],id[i]);
  }


  return 0;
}
