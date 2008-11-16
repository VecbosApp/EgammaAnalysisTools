//---------------------------------------
// Description:
// Class to estimate the performances of Likelihood Electron ID
//---------------------------------------

#ifndef LIKELIHOODANALYSIS_H
#define LIKELIHOODANALYSIS_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/EgammaBase.h"

class LikelihoodAnalysis : public EgammaBase {

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
  void estimateIDEfficiency(const char *outname);
  //! produce the mis-ID eta/pT distributions
  void estimateFakeRate(const char *outname="job0");

private:

  //! apply the custom offline electron ID
  bool getCustomEleID(int eleIndex, const char *fileCuts, const char *fileSwitches);

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;

};

#endif
