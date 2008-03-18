#include <iostream>

#include "TVector3.h"
#include "TH1F.h"

#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"


LikelihoodAnalysis::LikelihoodAnalysis(TTree *tree)
  : EgammaBase(tree) {
}



LikelihoodAnalysis::~LikelihoodAnalysis() {
}



float LikelihoodAnalysis::findEquivalentLHCut(float wantEfficiency) {

  int nbins = 500;
  TH1F *LHBinnedHisto = new TH1F("LHBinnedHisto", "LHBinnedHisto", nbins, 0.0, 1.0);


  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    for(int iele=0; iele<nEle; iele++) {
      if (eleClassEle[iele]==40) LHBinnedHisto->Fill(1.0);
      else {
	LHBinnedHisto->Fill(eleLikelihoodEle[iele]);
      }
    }

  }


  // scan the likelihood to find the cut
  float nEntries = LHBinnedHisto->GetEntries();
  float efficiency = 0.0;
  int bin = nbins+1;

  while (efficiency < wantEfficiency) {

    float integral = LHBinnedHisto->Integral(bin,nbins+1);
    efficiency = integral / nEntries;
    std::cout << "integral = " << integral 
	      << "\tefficiency = " << efficiency
	      << "bin = " << bin << std::endl;
    
    bin--;

  }

  float equivalentLHCut = LHBinnedHisto->GetBinLowEdge(bin);

  std::cout << "Equivalent cut on LH is LH > " << equivalentLHCut << std::endl
	    << "efficiency is: " << efficiency << std::endl
	    << "while wanted is: " << wantEfficiency << std::endl;

}




void LikelihoodAnalysis::reproduceEgammaCutID() {

  CutBasedEleIDSelector EgammaLooseCutBasedID;
  EgammaLooseCutBasedID.Configure("config/"); 

  int nEvtGoodElectron = 0;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    int iSelected = -1;
    for(int iele=0; iele<nEle; iele++) {

      TVector3 pTrkAtOuter(pxAtOuterEle[iele],pyAtOuterEle[iele],pzAtOuterEle[iele]);

      EgammaLooseCutBasedID.SetHOverE( eleHoEEle[iele] );
      EgammaLooseCutBasedID.SetS9S25( s9s25Ele[iele] );
      EgammaLooseCutBasedID.SetDEta( eleDeltaEtaAtVtxEle[iele] );
      EgammaLooseCutBasedID.SetDPhiIn( eleDeltaPhiAtVtxEle[iele] );
      EgammaLooseCutBasedID.SetDPhiOut( eleDeltaPhiAtCaloEle[iele] );
      EgammaLooseCutBasedID.SetInvEminusInvP( 1./eleCaloCorrEEle[iele]-1./eleTrackerPEle[iele] );
      EgammaLooseCutBasedID.SetBremFraction( fabs(eleTrackerPEle[iele]-pTrkAtOuter.Mag())/eleTrackerPEle[iele] );
      EgammaLooseCutBasedID.SetSigmaEtaEta( sqrt(covEtaEtaEle[iele]) );
      EgammaLooseCutBasedID.SetSigmaPhiPhi( sqrt(covPhiPhiEle[iele]) );
      EgammaLooseCutBasedID.SetEOverPout( eleCorrEoPoutEle[iele] );
      EgammaLooseCutBasedID.SetEOverPin( eleCorrEoPEle[iele] );
      EgammaLooseCutBasedID.SetElectronClass ( eleClassEle[iele] );
      EgammaLooseCutBasedID.SetEgammaCutBasedID ( eleIdCutBasedEle[iele] );
      
      bool isEleIDCutBased = EgammaLooseCutBasedID.output();

      if( isEleIDCutBased ) iSelected=iele;
    }

    nEvtGoodElectron++;

  }

  EgammaLooseCutBasedID.diplayEfficiencies();


}



void LikelihoodAnalysis::estimateIDEfficiency() {


}



void LikelihoodAnalysis::estimateFakeRate(const char *definition) {

}



bool LikelihoodAnalysis::getCustomEleID(int eleIndex, const char *fileCuts, const char *fileSwitches) {

}
