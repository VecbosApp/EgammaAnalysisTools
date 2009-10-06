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

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

// Offline analysis includes
#include "EgammaAnalysisTools/include/Application.hh"
#include "CommonTools/include/TriggerMask.hh"
#if Application == 1
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"
#endif
#if Application == 2
#include "EgammaAnalysisTools/include/LHPdfsProducer.hh"
#endif

int main(int argc, char* argv[]) {

  char inputFileName[150];
  char outputFileName[150];

  if ( argc < 2 ){
    std::cout << "missing argument: insert at least inputFile with list of root files" << std::endl; 
    std::cout << "EgammaAnalysis inputFile [outputFile]" << std::endl;
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  if ( argc < 3 ) sprintf(outputFileName,"pdfs.root");
  else strcpy(outputFileName,argv[2]);

  // -------------------------
  // loading file:
  TChain *theChain = new TChain("ntp1");
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

#if Application == 1

  LikelihoodAnalysis analysis(theChain);
  analysis.reproduceEgammaCutID();
  //  analysis.findEquivalentLHCut( 0.79935 );        // tight eleID
  //  analysis.findEquivalentLHCut( 0.957 );       // loose eleID
  //  analysis.estimateIDEfficiency();
  //  analysis.estimateFakeRate();
  
#endif

#if Application == 2
  
  LHPdfsProducer producer(theChain);

  TriggerMask maskSignal(treeCond);
  maskSignal.requireTrigger("HLT_Ele15_SW_L1R");
  maskSignal.requireTrigger("HLT_Ele15_SW_EleId_L1R");
  maskSignal.requireTrigger("HLT_Ele15_SW_LooseTrackIso_L1R");
  maskSignal.requireTrigger("HLT_Ele15_SiStrip_L1R");
  std::vector<int> requiredSignalTriggers = maskSignal.getBits();
  producer.requireSignalTrigger(requiredSignalTriggers);

  TriggerMask maskBackground(treeCond);
  maskBackground.requireTrigger("HLT_Ele15_SW_L1R");   // to be changed
  std::vector<int> requiredBackgroundTriggers = maskBackground.getBits();
  producer.requireBackgroundTrigger(requiredBackgroundTriggers);

  // producer.LoopZTagAndProbe("fromZ_TandP_");
  // producer.saveHistos("fromZ_TandP_histos");

  // producer.LoopZ("fromZ_MC_");
  // producer.saveHistos("fromZ_MC_histos");

  producer.LoopQCDTagAndProbe("fromQCD_TandP_");
  producer.saveHistos("fromQCD_TandP_histos");

#endif

  return 0;

}
