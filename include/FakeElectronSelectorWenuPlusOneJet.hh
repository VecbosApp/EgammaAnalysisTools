///---------------------------------------
// Description:
// Class to fill a tree of fake electron candidates to study ele ID
//---------------------------------------

#ifndef FAKEELECTRONSELECTORWENUPLUSONEJET_H
#define FAKEELECTRONSELECTORWENUPLUSONEJET_H

#include <vector>

#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/FakeTree.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "CommonTools/include/Counters.hh"

#include <TLorentzVector.h>
#include <TTree.h>


class FakeElectronSelectorWenuPlusOneJet : public Egamma {

public:
  //! constructor
  FakeElectronSelectorWenuPlusOneJet(TTree *tree=0);
  //! destructor
  virtual ~FakeElectronSelectorWenuPlusOneJet();
  //! do the loop over the events
  void Loop(const char *treefilesuffix);
  
private:

  bool _isData;

  // the tag and the probe
  int electrons[2];

  // counters
  Counters myCounter;

  /// MVA for electron ID. To be created and initialized from the children classes
  ElectronLikelihood *LH;
  ElectronIDMVA *fMVAHWW, *fMVAHWWWithIso;
  ElectronIDMVAHZZ *fMVAHZZDanV0, *fMVAHZZSiV0, *fMVAHZZSiV1, *fMVAHZZSiDanV2;
  ElectronIDMVAHZZ *fMVAHWWDanV0, *fMVAHWWSiV0, *fMVAHWWSiV1, *fMVAHWWSiDanV2;

  // outputs for the kinematics and the id
  FakeTree *myOutKineTree;
  RedEleIDTree *myOutIDTree;

};

#endif
