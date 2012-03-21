///---------------------------------------
// Description:
// Class to fill a tree of fake electron candidates to study ele ID
//---------------------------------------

#ifndef FAKEELECTRONSELECTORWMUNUPLUSONEJET_H
#define FAKEELECTRONSELECTORWMUNUPLUSONEJET_H

#include <vector>

#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/FakeTree.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "CommonTools/include/Counters.hh"

#include <TLorentzVector.h>
#include <TTree.h>


class FakeElectronSelectorWmunuPlusOneJet : public Egamma {

public:
  //! constructor
  FakeElectronSelectorWmunuPlusOneJet(TTree *tree=0);
  //! destructor
  virtual ~FakeElectronSelectorWmunuPlusOneJet();
  //! do the loop over the events
  void Loop(const char *treefilesuffix);
  
private:

  std::pair<int,int> getBestGoodMuonPair(std::vector<int> goodMuons);

  bool _isData;

  // the tag and the probe
  int electrons[2];

  // counters
  Counters myCounter;

  /// MVA for electron ID. To be created and initialized from the children classes
  ElectronLikelihood *LH;
  ElectronIDMVA *fMVAHWW, *fMVAHWWNoIP;
  ElectronIDMVAHZZ *fMVAHZZDanV0, *fMVAHZZSiV0, *fMVAHZZSiV1, *fMVAHZZSiDanV0;
  ElectronIDMVAHZZ *fMVAHWWDanV0, *fMVAHWWSiV0, *fMVAHWWSiV1, *fMVAHWWSiDanV0;

  // outputs for the kinematics and the id
  FakeTree *myOutKineTree;
  RedEleIDTree *myOutIDTree;

};

#endif
