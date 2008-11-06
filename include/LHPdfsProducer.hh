//---------------------------------------
// Description:
// Create the Likelihood Pdfs from Z->ee Tag and probe
//---------------------------------------

#ifndef LHPDFSPRODUCER_H
#define LHPDFSPRODUCER_H

#include <vector>

#include "TH1F.h"
#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/EgammaBase.h"

class LHPdfsProducer : public EgammaBase {

public:
  //! constructor
  LHPdfsProducer(TTree *tree=0);
  //! destructor
  virtual ~LHPdfsProducer();
  //! loop over events
  void Loop();
  //! set the list of the required triggers
  void requireTrigger(std::vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  //! save the pdfs in a ROOT file
  void saveHistos(const char *filename);

private:

  //! book the histograms of the PDFs
  void bookHistos();

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;

  //! the tag and the probe
  int electrons[2];

  //! the required trigger bits
  std::vector<int> m_requiredTriggers;

  //! the configurable selection
  Selection *m_selection;

  // ---------- monitoring histograms ------------
  TH1F *m_Zmass;

  /// Electrons: not splitted
  /// histo[ecalsubdet][ptbin]
  TH1F *dPhiCaloUnsplitEle[2][2];
  TH1F *dPhiVtxUnsplitEle[2][2];
  TH1F *dEtaUnsplitEle[2][2];
  TH1F *EoPoutUnsplitEle[2][2];
  TH1F *HoEUnsplitEle[2][2];  
  TH1F *shapeFisherUnsplitEle[2][2];  
  TH1F *sigmaEtaEtaUnsplitEle[2][2];
  TH1F *sigmaEtaPhiUnsplitEle[2][2];
  TH1F *sigmaPhiPhiUnsplitEle[2][2];
  TH1F *s1s9UnsplitEle[2][2];
  TH1F *s9s25UnsplitEle[2][2];
  TH1F *LATUnsplitEle[2][2];
  TH1F *etaLATUnsplitEle[2][2];
  TH1F *phiLATUnsplitEle[2][2];
  TH1F *a20UnsplitEle[2][2];
  TH1F *a42UnsplitEle[2][2];
  
  /// Electrons class-splitted
  /// histo[ecalsubdet][ptbin][class]
  TH1F *dPhiCaloClassEle[2][2][2];
  TH1F *dPhiVtxClassEle[2][2][2];
  TH1F *dEtaClassEle[2][2][2];
  TH1F *EoPoutClassEle[2][2][2];
  TH1F *HoEClassEle[2][2][2];
  TH1F *shapeFisherClassEle[2][2][2];
  TH1F *sigmaEtaEtaClassEle[2][2][2];
  TH1F *sigmaEtaPhiClassEle[2][2][2];
  TH1F *sigmaPhiPhiClassEle[2][2][2];
  TH1F *s1s9ClassEle[2][2][2];
  TH1F *s9s25ClassEle[2][2][2];
  TH1F *LATClassEle[2][2][2];
  TH1F *etaLATClassEle[2][2][2];
  TH1F *phiLATClassEle[2][2][2];
  TH1F *a20ClassEle[2][2][2];
  TH1F *a42ClassEle[2][2][2];

  /// Electrons fullclass-splitted
  /// histo[ecalsubdet][ptbin][fullclass]
  TH1F *dPhiCaloFullclassEle[2][2][4];
  TH1F *dPhiVtxFullclassEle[2][2][4];
  TH1F *dEtaFullclassEle[2][2][4];
  TH1F *EoPoutFullclassEle[2][2][4];
  TH1F *HoEFullclassEle[2][2][4];
  TH1F *shapeFisherFullclassEle[2][2][4];
  TH1F *sigmaEtaEtaFullclassEle[2][2][4];
  TH1F *sigmaEtaPhiFullclassEle[2][2][4];
  TH1F *sigmaPhiPhiFullclassEle[2][2][4];
  TH1F *s1s9FullclassEle[2][2][4];
  TH1F *s9s25FullclassEle[2][2][4];
  TH1F *LATFullclassEle[2][2][4];
  TH1F *etaLATFullclassEle[2][2][4];
  TH1F *phiLATFullclassEle[2][2][4];
  TH1F *a20FullclassEle[2][2][4];
  TH1F *a42FullclassEle[2][2][4];

};

#endif
