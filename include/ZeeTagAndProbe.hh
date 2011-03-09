///---------------------------------------
// Description:
// Class to estimate the performances of Likelihood Electron ID
//---------------------------------------

#ifndef ZTAGANDPROBE_H
#define ZTAGANDPROBE_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/CiCBasedEleSelector.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include <TLorentzVector.h>

class ZeeTagAndProbe : public Egamma {

public:
  //! constructor
  ZeeTagAndProbe(TTree *tree=0);
  //! destructor
  virtual ~ZeeTagAndProbe();
  //! do the loop over the events
  void Loop(const char *treefilesuffix);
  
private:

  //! configurations
  void configSelection(Selection* selection);

  //! apply the custom offline electron ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! apply the custom offline electron ID
  void isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);
  
  float SigmaiEiE(int electron);
  float SigmaiPiP(int electron);
  
  //! class members
  std::vector<CutBasedEleIDSelector> EgammaCutBasedID;
  std::vector<CiCBasedEleSelector> EgammaCiCBasedID;
  std::vector<CutBasedEleIDSelector> EgammaLHBasedID;

  std::vector<std::string> EgammaCutBasedIDWPs;
  std::vector<std::string> EgammaCiCBasedIDWPs;
  std::vector<std::string> EgammaLHBasedIDWPs;

  ElectronLikelihood *LH;

  bool isData_;

  //! the configurable selection
  Selection *m_selection;

  // the tag and the probe
  int electrons[2];

};

#endif
