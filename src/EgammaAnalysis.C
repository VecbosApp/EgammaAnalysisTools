//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
// Original code:
//    CutAnaHiggs_2e2nu.cpp
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdlib>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

// Offline analysis includes
#include "EgammaAnalysisTools/include/Application.hh"
#include "CommonTools/include/TriggerMask.hh"
// #if Application == 1
// #include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"
// #endif
#if Application == 2
#include "EgammaAnalysisTools/include/LHPdfsProducer.hh"
#endif
#if Application == 3
#include "EgammaAnalysisTools/include/sPlotsPdfsComparison.hh"
#endif
#if Application == 4
#include "EgammaAnalysisTools/include/SuperClusterWSelection.hh"
#endif
#if Application == 5
#include "EgammaAnalysisTools/include/sPlotsPdfsComparison.h"
#endif

int main(int argc, char* argv[]) {

  char inputFileName[300];
  char outputFileName[300];

#if Application != 5
  if ( argc < 4 ){
    std::cout << "missing argument: insert at least inputFile with list of root files" << std::endl; 
    std::cout << "EgammaAnalysis inputFile outputFile 0" << std::endl;
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  strcpy(outputFileName,argv[2]);
  int signal = atoi(argv[3]);

  // -------------------------
  // loading file:
  TChain *theChain = 0;

#if Application != 3
  theChain = new TChain("ntp1");
#else 
  theChain = new TChain(outputFileName);
#endif

  char Buffer[500];
  char MyRootFile[2000];  
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);

  // get the tree with the conditions from the first file
  TTree *treeCond = new TTree();

  int nfiles=1;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	theChain->Add(MyRootFile);
	if ( nfiles==1 ) {
	  TFile *firstfile = TFile::Open(MyRootFile);
	  treeCond = (TTree*)firstfile->Get("Conditions");
	}
	std::cout << "chaining " << MyRootFile << std::endl;
	nfiles++;
      }
  }
  inputFile->close();
  delete inputFile;
#endif

// #if Application == 1

//   LikelihoodAnalysis analysis(theChain);
//   analysis.reproduceEgammaCutID();
//   //  analysis.findEquivalentLHCut( 0.79935 );        // tight eleID
//   //  analysis.findEquivalentLHCut( 0.957 );       // loose eleID
//   //  analysis.estimateIDEfficiency();
//   //  analysis.estimateFakeRate();
  
// #endif

#if Application == 2

  char title[1000];
  
  LHPdfsProducer producer(theChain);
  
  // TriggerMask maskSignal(treeCond);
  // maskSignal.requireTrigger("HLT_Ele15_SW_L1R");
  // maskSignal.requireTrigger("HLT_Ele20_SW_L1R");
  // std::vector<int> requiredSignalTriggers = maskSignal.getBits();
  // producer.requireSignalTrigger(requiredSignalTriggers);

  std::vector<std::string> mask;
  mask.push_back("HLT_Jet15U");   
  producer.setRequiredTriggers(mask);
  
  //  sprintf(title,"%s_zTandP_tree.root",outputFileName);  
  // producer.LoopZTagAndProbe(title);
  // sprintf(title,"%s_zTandP_histos.root",outputFileName);    
  // producer.saveHistos(title);
  
  // sprintf(title,"%s_zMC_tree.root",outputFileName);  
  // producer.LoopZTagAndProbeForMcTruth(title);
  // sprintf(title,"%s_zMC_histos.root",outputFileName);    
  // producer.saveHistos(title);

  sprintf(title,"%s_qcdTandP_tree.root",outputFileName);  
  producer.LoopQCDTagAndProbe(title);
  sprintf(title,"%s_qcdTandP_counters.root",outputFileName);
  producer.displayEfficiencies(title);
  sprintf(title,"%s_qcdTandP_histos.root",outputFileName);    
  producer.saveHistos(title);

#endif

#if Application == 3

  sPlotsPdfsComparison p(theChain);
  p.RunOverMC(false);
  p.Loop();

#endif

#if Application == 4

  SuperClusterWSelection p(theChain);
  p.setPrefix(outputFileName);
  p.setSignal(signal);
  p.Loop();
  p.displayEfficiencies();

#endif

#if Application == 5

  int doSignalsPlots = 0;
  int isMC = 0;
  
  TFile *fileData, *fileMC;

  if(doSignalsPlots) {
    fileData = TFile::Open((std::string("results_data/sPlots/Wenu_tree.root")).c_str());
    fileMC = TFile::Open((std::string("results/treesW/WJetsMADGRAPH_Wenu.root")).c_str());
  } else {
    fileData = TFile::Open((std::string("results_data/sPlots/Wenu_bkgFit_tree.root")).c_str());
    fileMC = TFile::Open((std::string("results/treesW/QCD_Wenu.root")).c_str());    
  }

  TTree *tree = 0;
  if(isMC) tree = (TTree*) fileMC->Get("T1");
  else tree = (TTree*) fileData->Get("dataset");

  sPlotsPdfsComparison p;
  p.Init(tree, isMC);
  p.doSignalsPlots(doSignalsPlots);
  p.Loop();

#endif

  return 0;

}
