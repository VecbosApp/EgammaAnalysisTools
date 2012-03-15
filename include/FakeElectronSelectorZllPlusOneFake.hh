///---------------------------------------
// Description:
// Class to fill a tree of fake electron candidates to study ele ID
//---------------------------------------

#ifndef FAKEELECTRONSELECTORZLLPLUSONEFAKE_H
#define FAKEELECTRONSELECTORZLLPLUSONEFAKE_H

#include <vector>

#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/FakeTree.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "CommonTools/include/Counters.hh"

#include <TLorentzVector.h>
#include <TTree.h>


class FakeElectronSelectorZllPlusOneFake : public Egamma {

public:
  //! constructor
  FakeElectronSelectorZllPlusOneFake(TTree *tree=0);
  //! destructor
  virtual ~FakeElectronSelectorZllPlusOneFake();
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
  ElectronIDMVAHZZ *fMVAHZZ, *fMVAHZZNoIP, *fMVAHZZMC;

  // outputs for the kinematics and the id
  FakeTree *myOutKineTree;
  RedEleIDTree *myOutIDTree;

};

#endif
