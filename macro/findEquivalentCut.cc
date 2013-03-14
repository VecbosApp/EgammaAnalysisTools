#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <iostream>
#include <sstream>

using namespace std;

float getCut(TString presel, TString cut1, TString var2, float min2, float max2, TString cutDir) {
  
  TFile *file = 0;
  TTree *tree = 0;
  
  file = TFile::Open("results_data/fakes.root");
  if( file ) {
    file->cd();
    tree = (TTree*)file->Get("eleIDdir/T1");
  } else {
    cout << "File results_data/fakes.root not existing !" << endl;
    return 9999.;
  }
  if(!tree) {
    cout << "Tree eleIDdir/T1 not existing inside file!" << endl;
    return 9999.;
  }

  tree->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", "results_data/fakes_hzzisoFriend.root" );
  tree->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", "results_data/fakes_hzzidbitsFriend.root" );

  // first, find the efficiency on the original variable
  TString fullCut1 = presel + TString(" && ") + cut1;
  TH1F *h1 = new TH1F("h1","",10,-1000,1000);
  tree->Project("h1","ecaldriven",presel);
  float den1 = h1->Integral();
  tree->Project("h1","ecaldriven",fullCut1);
  float num1 = h1->Integral();
  float eff1 = num1/den1;

  // and now find the equivalent cut scanning var2
  TH1F *h2 = new TH1F("h2","",500,min2,max2);
  TAxis *axis2 = h2->GetXaxis();
  tree->Project("h2",var2,presel);
  int nBins2 = axis2->GetNbins();

  double Integral2 = h2->Integral(0,nBins2+1);
  double tmpIntegral2=0.0;

  float runningeff, runningeffprev;
  runningeff = runningeffprev = -1;
  
  for ( int ibin=0; ibin<=nBins2+1; ibin++) {
    
    if( strcmp(cutDir,"<")==0 )
      tmpIntegral2 = h2->Integral(0,ibin);
    else if( strcmp(cutDir,">")==0 )
      tmpIntegral2 = h2->Integral(ibin,nBins2+1);
    else {
      std::cout << "CONFIGURATION ERROR! direction of the cut not set." << std::endl
		<< "Please use: \">\" for var>x0 or  \"<\" for var<x0" << std::endl;
      return 0;
    }
    
    runningeffprev=runningeff;
    runningeff=tmpIntegral2/Integral2;

    if((eff1>runningeffprev && eff1<=runningeff && runningeffprev!=-1) || (eff1>=runningeff && eff1<runningeffprev)) {
      float cut = axis2->GetBinCenter(ibin);
      cout << "Equivalent cut for preselection: " << presel << " is " << var2 << cutDir << cut << endl;
      // cout << "eff(orig cut) = " << eff1 << "\teff(new cut)" << runningeff << endl;
      return cut;
    }
  }

  return 9999.;

}

float getCutFixedEff(float eff, TString presel, TString var2, float min2, float max2, TString cutDir, bool doFR=false) {
  
  TFile *file = 0;
  TTree *tree = 0;
  
  string filestr = "results_data/electrons_zeemc.root";
  string fileiso = "results_data/electrons_zeemc_hzzisoFriend.root";
  string fileidbits = "results_data/electrons_zeemc_hzzidbitsFriend.root";
  if(doFR) {
    filestr = "results_data/fakes-zll1e.root";
    fileiso = "results_data/fakes-zll1e_hzzisoFriend.root";
    fileidbits = "results_data/fakes-zll1e_hzzidbitsFriend.root";
  }

  file = TFile::Open(filestr.c_str());
  if( file ) {
    file->cd();
    tree = (TTree*)file->Get("eleIDdir/T1");
  } else {
    cout << "File " << file << " not existing !" << endl;
    return 9999.;
  }
  if(!tree) {
    cout << "Tree eleIDdir/T1 not existing inside file!" << endl;
    return 9999.;
  }

  tree->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", fileiso.c_str() );
  tree->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", fileidbits.c_str() );

  // and now find the equivalent cut scanning var2
  TH1F *h2 = new TH1F("h2","",500,min2,max2);
  TAxis *axis2 = h2->GetXaxis();
  tree->Project("h2",var2,presel);
  int nBins2 = axis2->GetNbins();

  double Integral2 = h2->Integral(0,nBins2+1);
  double tmpIntegral2=0.0;

  float runningeff, runningeffprev;
  runningeff = runningeffprev = -1;
  
  for ( int ibin=0; ibin<=nBins2+1; ibin++) {
    
    if( strcmp(cutDir,"<")==0 )
      tmpIntegral2 = h2->Integral(0,ibin);
    else if( strcmp(cutDir,">")==0 )
      tmpIntegral2 = h2->Integral(ibin,nBins2+1);
    else {
      std::cout << "CONFIGURATION ERROR! direction of the cut not set." << std::endl
		<< "Please use: \">\" for var>x0 or  \"<\" for var<x0" << std::endl;
      return 0;
    }
    
    runningeffprev=runningeff;
    runningeff=tmpIntegral2/Integral2;

    if(runningeff>eff-0.005 && runningeff<eff+0.005) {
      float cut = axis2->GetBinCenter(ibin);
      cout << "Equivalent cut for preselection: " << presel << " is " << var2 << cutDir << cut << endl;
      // cout << "eff(orig cut) = " << eff1 << "\teff(new cut)" << runningeff << endl;
      return cut;
    }
  }

  return 9999.;

}

float getEff(TString presel, TString var2, float cut, float min2, float max2, TString cutDir, bool doFR=false) {
  
  TFile *file = 0;
  TTree *tree = 0;
  
  string filestr = "results_data/electrons_zeemc.root";
  string fileiso = "results_data/electrons_zeemc_hzzisoFriend.root";
  string fileidbits = "results_data/electrons_zeemc_hzzidbitsFriend.root";
  if(doFR) {
    filestr = "results_data/fakes-zll1e.root";
    fileiso = "results_data/fakes-zll1e_hzzisoFriend.root";
    fileidbits = "results_data/fakes-zll1e_hzzidbitsFriend.root";
  }

  file = TFile::Open(filestr.c_str());
  if( file ) {
    file->cd();
    tree = (TTree*)file->Get("eleIDdir/T1");
  } else {
    cout << "File " << file << " not existing !" << endl;
    return 9999.;
  }
  if(!tree) {
    cout << "Tree eleIDdir/T1 not existing inside file!" << endl;
    return 9999.;
  }

  tree->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", fileiso.c_str() );
  tree->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", fileidbits.c_str() );


  // and now find the equivalent cut scanning var2
  TH1F *h1 = new TH1F("h1","",500,min2,max2);
  TH1F *h2 = new TH1F("h2","",500,min2,max2);
  TAxis *axis2 = h2->GetXaxis();
  tree->Project("h2",var2,presel);
  int nBins2 = axis2->GetNbins();

  double Integral2 = h2->Integral(0,nBins2+1);


  stringstream cutNom;
  cutNom << presel << " && " << var2.Data() << cutDir.Data() << cut;
  tree->Project("h1",var2,cutNom.str().c_str());
  double Integral1 = h1->Integral(0,nBins2+1);
  float eff = Integral1/Integral2;

  return eff;

}

void lookForCutsIso() {
  
  // pT > 20 GeV (<35 GeV is to keep EWK contamination engligible)
  getCut("DenomFakeSmurf && abs(eta) < 1.479 && pt>20 && pt<35","combPFIsoHWW/pt<0.13",
	 "combPFIsoHZZ/pt",-0.1,2.0,"<");
  getCut("DenomFakeSmurf && abs(eta) > 1.479 && pt>20 && pt<35","combPFIsoHWW/pt<0.09",
	 "combPFIsoHZZ/pt",-0.1,2.0,"<");

  // pT < 20 GeV 
  getCut("DenomFakeSmurf && abs(eta) < 1.479 && pt<20","combPFIsoHWW/pt<0.13",
	 "combPFIsoHZZ/pt",-0.1,2.0,"<");
  getCut("DenomFakeSmurf && abs(eta) > 1.479 && pt<20","combPFIsoHWW/pt<0.09",
	 "combPFIsoHZZ/pt",-0.1,2.0,"<");

}


void lookForCutsNoEAIso() {
  
  // barrel (pt<35 GeV is to keep EWK contamination engligible)
  getCut("DenomFakeSmurf && abs(eta) < 1.479 && pt<35","combPFIsoHWW/pt<0.13",
	 "combPFIsoHZZ/pt",0.0,2.0,"<");
  // endcap (pt<35 GeV is to keep EWK contamination engligible)
  getCut("DenomFakeSmurf && abs(eta) >= 1.479 && pt<35","combPFIsoHWW/pt<0.09",
	 "combPFIsoHZZ/pt",0.0,2.0,"<");

}

void lookForCutsId() {
  getCut("DenomFakeSmurf && mcmatch && combPFIsoHWW/pt< 0.10 && pt<35 && pt < 20 && abs(eta)<1.0","bdthww[0]>0.139","newbdthww[3]",-1.0,1.0,">");
  getCut("DenomFakeSmurf && mcmatch && combPFIsoHWW/pt< 0.10 && pt<35 && pt < 20 && abs(eta)>=1.0 && abs(eta)<1.479","bdthww[0]>0.525","newbdthww[3]",-1.0,1.0,">");
  getCut("DenomFakeSmurf && mcmatch && combPFIsoHWW/pt< 0.10 && pt<35 && pt < 20 && abs(eta)>=1.479 && abs(eta)<2.5","bdthww[0]>0.543","newbdthww[3]",-1.0,1.0,">");

  getCut("DenomFakeSmurf && mcmatch && combPFIsoHWW/pt< 0.10 && pt<35 && pt >= 20 && abs(eta)<1.0","bdthww[0]>0.947","newbdthww[3]",-1.0,1.0,">");
  getCut("DenomFakeSmurf && mcmatch && combPFIsoHWW/pt< 0.10 && pt<35 && pt >= 20 && abs(eta)>=1.0 && abs(eta)<1.479","bdthww[0]>0.950","newbdthww[3]",-1.0,1.0,">");
  getCut("DenomFakeSmurf && mcmatch && combPFIsoHWW/pt< 0.10 && pt<35 && pt >= 20 && abs(eta)>=1.479 && abs(eta)<2.5","bdthww[0]>0.884","newbdthww[3]",-1.0,1.0,">");
}

void lookForCutsCiCId() {
  getCut("abs(eta)<0.8 && pt<10 && missHits<=1","cicid[3]==1","bdthzz[3]",-1.0,1.0,">");
  getCut("abs(eta)>=0.8 && abs(eta)<1.479 && pt<10 && missHits<=1","cicid[3]==1","bdthzz[3]",-1.0,1.0,">");
  getCut("abs(eta)>1.479 && pt<10 && missHits<=1","cicid[3]==1","bdthzz[3]",-1.0,1.0,">");

  getCut("abs(eta)<0.8 && pt>=10 && missHits<=1","cicid[3]==1","bdthzz[3]",-1.0,1.0,">");
  getCut("abs(eta)>=0.8 && abs(eta)<1.479 && pt>=10 && missHits<=1","cicid[3]==1","bdthzz[3]",-1.0,1.0,">");
  getCut("abs(eta)>1.479 && pt>=10 && missHits<=1","cicid[3]==1","bdthzz[3]",-1.0,1.0,">");

  // Equivalent cut for preselection: abs(eta)<1.0 && pt<10 is bdthzz[3]>0.45
  // Equivalent cut for preselection: abs(eta)>=1.0 && abs(eta)<1.479 && pt<10 is bdthzz[3]>-0.074
  // Equivalent cut for preselection: abs(eta)>1.479 && pt<10 is bdthzz[3]>0.47
  // Equivalent cut for preselection: abs(eta)<1.0 && pt>=10 is bdthzz[3]>0.546
  // Equivalent cut for preselection: abs(eta)>=1.0 && abs(eta)<1.479 && pt>=10 is bdthzz[3]>0.118
  // Equivalent cut for preselection: abs(eta)>1.479 && pt>=10 is bdthzz[3]>0.802

}

void lookForCutsCiCIso() {
  getCut("abs(eta)<1.0 && pt<10","cicmediumiso==1","combPFIsoHZZ/pt",-1.0,1.0,"<");
  getCut("abs(eta)>=1.0 && abs(eta)<1.479 && pt<10","cicmediumiso==1","combPFIsoHZZ/pt",-1.0,1.0,"<");
  getCut("abs(eta)>1.479 && pt<10","cicmediumiso==1","combPFIsoHZZ/pt",-1.0,1.0,"<");

  getCut("abs(eta)<1.0 && pt>=10","cicmediumiso==1","combPFIsoHZZ/pt",-1.0,1.0,"<");
  getCut("abs(eta)>=1.0 && abs(eta)<1.479 && pt>=10","cicmediumiso==1","combPFIsoHZZ/pt",-1.0,1.0,"<");
  getCut("abs(eta)>1.479 && pt>=10","cicmediumiso==1","combPFIsoHZZ/pt",-1.0,1.0,"<");

  // Equivalent cut for preselection: abs(eta)<1.0 && pt<10 is combPFIsoHZZ/pt<0.106
  // Equivalent cut for preselection: abs(eta)>=1.0 && abs(eta)<1.479 && pt<10 is combPFIsoHZZ/pt<0.106
  // Equivalent cut for preselection: abs(eta)>1.479 && pt<10 is combPFIsoHZZ/pt<0.002
  // Equivalent cut for preselection: abs(eta)<1.0 && pt>=10 is combPFIsoHZZ/pt<0.222
  // Equivalent cut for preselection: abs(eta)>=1.0 && abs(eta)<1.479 && pt>=10 is combPFIsoHZZ/pt<0.198
  // Equivalent cut for preselection: abs(eta)>1.479 && pt>=10 is combPFIsoHZZ/pt<0.118

}


void lookForCutsFixedEff() {
  getCutFixedEff(0.9,"abs(eta)<0.8 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
  getCutFixedEff(0.9,"abs(eta)>=0.8 && abs(eta)<1.479 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
  getCutFixedEff(0.9,"abs(eta)>1.479 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");

  getCutFixedEff(0.97,"abs(eta)<0.8 && pt>=10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
  getCutFixedEff(0.97,"abs(eta)>=0.8 && abs(eta)<1.479 && pt>=10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
  getCutFixedEff(0.97,"abs(eta)>1.479 && pt>=10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");

  // Equivalent cut for preselection: abs(eta)<1.0 && pt<10 is bdthzz[3]>0.45
  // Equivalent cut for preselection: abs(eta)>=1.0 && abs(eta)<1.479 && pt<10 is bdthzz[3]>-0.074
  // Equivalent cut for preselection: abs(eta)>1.479 && pt<10 is bdthzz[3]>0.47
  // Equivalent cut for preselection: abs(eta)<1.0 && pt>=10 is bdthzz[3]>0.546
  // Equivalent cut for preselection: abs(eta)>=1.0 && abs(eta)<1.479 && pt>=10 is bdthzz[3]>0.118
  // Equivalent cut for preselection: abs(eta)>1.479 && pt>=10 is bdthzz[3]>0.802

}

void scanFixedEff() {

  /*
  // eta1
  for(float eff=0.95; eff>0.75; eff-=0.02) {
    cout << "probing eff = " << eff << endl;
    float cut = getCutFixedEff(eff,"abs(eta)<0.8 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
    float fr = getEff("abs(eta)<0.8 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4","bdthzz[3]",cut,-1.0,1.0,">",true);
    cout << "Fake rate = " << fr << endl;
  }
  */

  // eta2
  for(float eff=0.95; eff>0.75; eff-=0.02) {
    cout << "probing eff = " << eff << endl;
    float cut = getCutFixedEff(eff,"abs(eta)>0.8 && abs(eta)<1.479 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
    float fr = getEff("abs(eta)>0.8 && abs(eta)<1.479 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4","bdthzz[3]",cut,-1.0,1.0,">",true);
    cout << "Fake rate = " << fr << endl;
  }

  // eta3
  for(float eff=0.95; eff>0.75; eff-=0.02) {
    cout << "probing eff = " << eff << endl;
    float cut = getCutFixedEff(eff,"abs(eta)>1.479 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4 && mcmatch","bdthzz[3]",-1.0,1.0,">");
    float fr = getEff("abs(eta)>1.479 && pt<10 && missHits<=1 && combPFIsoHZZ/pt<0.4","bdthzz[3]",cut,-1.0,1.0,">",true);
    cout << "Fake rate = " << fr << endl;
  }

}
