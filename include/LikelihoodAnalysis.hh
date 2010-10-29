//---------------------------------------
// Description:
// Class to estimate the performances of Likelihood Electron ID
//---------------------------------------

#ifndef LIKELIHOODANALYSIS_H
#define LIKELIHOODANALYSIS_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/CiCBasedEleSelector.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"

class LikelihoodAnalysis : public Egamma {

public:
  //! constructor
  LikelihoodAnalysis(TTree *tree=0);
  //! destructor
  virtual ~LikelihoodAnalysis();
  //! find the cut on LH that equals the electron ID with cuts (standard egamma)
  float findEquivalentLHCut();
  //! find the cut on LH that gives a wanted efficiency
  float findEquivalentLHCut(float wantEfficiency);
  //! reproduce the egamma cut based ID (for debugging)
  void reproduceEgammaCutID();
  //! produce the ID efficiency eta/pT distributions
  void estimateIDEfficiency(const char *outname="job0");
  //! produce the mis-ID eta/pT distributions
  void estimateFakeRate(const char *outname="job0");
  //! produce the mis-ID eta/pT distributions from QCD di-jets
  void estimateFakeRateQCD(const char *outname="job0");

private:

  //! apply the custom offline electron ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! apply the custom offline electron ID
  void isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);
  
  float SigmaiEiE(int electron);
  float SigmaiPiP(int electron);
  
  float likelihoodRatio(int eleIndex, ElectronLikelihood &lh);

  std::vector<CutBasedEleIDSelector> EgammaCutBasedID;
  std::vector<CiCBasedEleSelector> EgammaCiCBasedID;
  std::vector<CutBasedEleIDSelector> EgammaLHBasedID;

  std::vector<std::string> EgammaCutBasedIDWPs;
  std::vector<std::string> EgammaCiCBasedIDWPs;
  std::vector<std::string> EgammaLHBasedIDWPs;

  ElectronLikelihood *LH;

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;
  
  bool _isData;

};

#endif
