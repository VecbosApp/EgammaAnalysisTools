/// The Egamma class is an auxiliary class which contains basic
/// functionality useful for any analysis of Vecbos+jets events.
/// It derives from EgammaBase.

#ifndef Egamma_h
#define Egamma_h

#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "EgammaAnalysisTools/include/EgammaBase.h"
#include "EgammaAnalysisTools/include/ElectronIDMVA.h"
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
  //check for matching HLT object                          
  bool triggerMatch(float eta, float phi, float Dr);
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
  /// return the value of electron BDT
  float eleBDT(ElectronIDMVA *mva, int iele);
  /// apply the BDT cut
  bool passEleBDT(float pt, float eta, float bdtoutput);

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
  /// dxy, dz and dsz parameters with respect to PV for tracks
  double trackDxyPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom);
  double trackDzPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom);
  double trackDszPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom);
  /// dxy, dz and dsz parameters with respect to PV for electrons
  double eleDxyPV(int iele, int iPV);
  double eleDzPV(int iele, int iPV);
  double eleDszPV(int iele, int iPV);

  /// some useful template  
  template <typename T>
    T DeltaR(T eta1, T phi1, T eta2, T phi2); //< Delta R in between two pairs (eta,phi)
  
  template <typename T>
    T DeltaPhi(T phi1, T phi2); //< Delta phi in radians in between two angles.  
};

template <typename T>
T Egamma::DeltaPhi(T phi1, T phi2) {
  T result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

template <typename T>
T Egamma::DeltaR(T eta1, T phi1, T eta2, T phi2) {
  T dphi = DeltaPhi(phi1,phi2);
  T result = sqrt((eta1-eta2)*(eta1-eta2)+dphi*dphi);
  return result;
}

#endif
