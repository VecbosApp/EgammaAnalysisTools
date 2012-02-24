#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

float getCut(TString presel, TString cut1, TString var2, float min2, float max2, TString cutDir) {
  
  TFile *file = 0;
  TTree *tree = 0;
  
  file = TFile::Open("results_data_fakes/merged.root");
  if( file ) {
    file->cd();
    tree = (TTree*)file->Get("eleIDdir/T1");
  } else {
    cout << "File results_data_fakes/merged.root not existing !" << endl;
    return 9999.;
  }

  if(!tree) {
    cout << "Tree eleIDdir/T1 not existing inside file!" << endl;
    return 9999.;
  }

  // first, find the efficiency on the original variable
  TString fullCut1 = presel + TString(" && ") + cut1;
  TH1F *h1 = new TH1F("h1","",10,-1000,1000);
  tree->Project("h1","ecaldriven",presel);
  float den1 = h1->Integral();
  tree->Project("h1","ecaldriven",fullCut1);
  float num1 = h1->Integral();
  float eff1 = num1/den1;

  // and now find the equivalent cut scanning var2
  TH1F *h2 = new TH1F("h2","",200,min2,max2);
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

    if((eff1>runningeffprev && eff1<runningeff && runningeffprev!=-1) || (eff1>runningeff && eff1<runningeffprev)) {
      float cut = axis2->GetBinCenter(ibin);
      cout << "Equivalent cut for preselection: " << presel << " is " << var2 << cutDir << cut << endl;
      //      cout << "eff(orig cut) = " << eff1 << "\teff(new cut)" << runningeff << endl;
      return cut;
    }
  }

  return 9999.;

}


void lookForCutsIso() {
  
  getCut("abs(eta) <  1.0","combPFIsoHWW/pt<0.13",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.18 * rho))/pt",-0.1,2.0,"<");
  getCut("abs(eta)>1.0 && abs(eta)<1.479","combPFIsoHWW/pt<0.13",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.19 * rho))/pt",-0.1,2.0,"<");
  getCut("abs(eta)>1.479 && abs(eta)<2.0","combPFIsoHWW/pt<0.09",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.21 * rho))/pt",-0.1,2.0,"<");
  getCut("abs(eta)>2.0 && abs(eta)<2.2","combPFIsoHWW/pt<0.09",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.38 * rho))/pt",-0.1,2.0,"<");
  getCut("abs(eta)>2.2 && abs(eta)<2.3","combPFIsoHWW/pt<0.09",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.61 * rho))/pt",-0.1,2.0,"<");
  getCut("abs(eta)>2.3 && abs(eta)<2.4","combPFIsoHWW/pt<0.09",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.73 * rho))/pt",-0.1,2.0,"<");
  getCut("abs(eta)>2.4","combPFIsoHWW/pt<0.09",
	 "(chaPFIso + (neuPFIso+phoPFIso - 0.90 * rho))/pt",-0.1,2.0,"<");


  // on 2011 fake data:
  // Equivalent cut for preselection: abs(eta) <  1.0 is (chaPFIso + (neuPFIso+phoPFIso - 0.18 * rho))/pt<0.2045
  // Equivalent cut for preselection: abs(eta)>1.0 && abs(eta)<1.479 is (chaPFIso + (neuPFIso+phoPFIso - 0.19 * rho))/pt<0.2045
  // Equivalent cut for preselection: abs(eta)>1.479 && abs(eta)<2.0 is (chaPFIso + (neuPFIso+phoPFIso - 0.21 * rho))/pt<0.1205
  // Equivalent cut for preselection: abs(eta)>2.0 && abs(eta)<2.2 is (chaPFIso + (neuPFIso+phoPFIso - 0.38 * rho))/pt<0.1205
  // Equivalent cut for preselection: abs(eta)>2.2 && abs(eta)<2.3 is (chaPFIso + (neuPFIso+phoPFIso - 0.61 * rho))/pt<0.0575
  // Equivalent cut for preselection: abs(eta)>2.3 && abs(eta)<2.4 is (chaPFIso + (neuPFIso+phoPFIso - 0.73 * rho))/pt<0.0785
  // Equivalent cut for preselection: abs(eta)>2.4 is (chaPFIso + (neuPFIso+phoPFIso - 0.90 * rho))/pt<0.0155

}


void lookForCutsHZZId() {

  getCut("pt < 20 && abs(eta)<1.0","bdthww>0.139","bdthzz",-0.3,0.3,">");
  getCut("pt < 20 && abs(eta)>=1.0 && abs(eta)<1.479","bdthww>0.525","bdthzz",-0.3,0.3,">");
  getCut("pt < 20 && abs(eta)>=1.479 && abs(eta)<2.5","bdthww>0.543","bdthzz",-0.3,0.3,">");

  getCut("pt >= 20 && abs(eta)<1.0","bdthww>0.947","bdthzz",-0.3,0.3,">");
  getCut("pt >= 20 && abs(eta)>=1.0 && abs(eta)<1.479","bdthww>0.950","bdthzz",-0.3,0.3,">");
  getCut("pt >= 20 && abs(eta)>=1.479 && abs(eta)<2.5","bdthww>0.884","bdthzz",-0.3,0.3,">");

}
