/// The Egamma class is an auxiliary class which contains basic
/// functionality useful for any analysis of Vecbos+jets events.
/// It derives from EgammaBase.

#ifndef Egamma_h
#define Egamma_h

#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "EgammaAnalysisTools/include/EgammaBase.h"
// ROOT includes
#include <TLorentzVector.h>
#include <TVector3.h>
// std includes
#include <string>
#include <vector>
#include <map>

class Egamma : public EgammaBase{

public:
  typedef std::pair<unsigned int,unsigned int> aLSSegment;
  typedef std::vector< std::pair<unsigned int,unsigned int> > LSSegments;
  typedef unsigned int aRun;
  typedef std::map< aRun, LSSegments > runsLSSegmentsMap;
  typedef std::pair < aRun, LSSegments > aRunsLSSegmentsMapElement;

  /// Class Constructor
  Egamma(TTree *tree=0);
  /// Class Destructor
  virtual ~Egamma();

  /// Fill RunLSMap according to json file
  void fillRunLSMap();
  /// Set Good Run LS
  void setJsonGoodRunList(const std::string& jsonFilePath);
  /// check if Run/LS is a good one
  bool isGoodRunLS();
  /// reload TriggerMask if necessary (data file is changed). Should be called for each event inside the event loop
  bool reloadTriggerMask(bool newVersion=false);
  /// set the list of required trigger to produce the bitmask
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers);
  //check if the event passed HLT. To be called per event
  bool hasPassedHLT();
  //get the value of the requested bits
  std::vector<int> getHLTOutput();

  /// Get pt given x/y coordinates
  float GetPt(float px, float py) { return TMath::Sqrt(px*px + py*py); }

  // useful electron functions
  /// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
  float SigmaiEiE(int electron);
  /// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
  float SigmaiPiP(int electron);
  // get the likelihood electron ID
  float likelihoodRatio(int eleIndex, ElectronLikelihood &lh);

protected:
  ///goodRUN/LS list
  runsLSSegmentsMap goodRunLS; 
  std::string jsonFile;

  std::string lastFile;
  std::vector<std::string> requiredTriggers;

  //! the list of required triggers
  std::vector<int> m_requiredTriggers;

  /// calculate transverse mass
  /// definitions in http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=104213
  float mT3(TLorentzVector pl1, TLorentzVector pl2, TVector3 met);

};

#endif
