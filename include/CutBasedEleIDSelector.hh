#ifndef CutBasedEleIDSelector_h
#define CutBasedEleIDSelector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class CutBasedEleIDSelector {

public:

  //! constructor
  CutBasedEleIDSelector();
  //! destructor
  virtual ~CutBasedEleIDSelector();   
  //! configure from files (class dependent)
  void Configure(const char *configDir);
  //! configure from files (class independent)
  void ConfigureNoClass(const char *configDir);
  //! set event by event observables
  void SetHOverE(float HOverE) { m_HOverE = HOverE; m_useHOverE = true; }
  void SetS9S25(float S9S25) { m_S9S25 = S9S25; m_useS9S25 = true; }
  void SetDEta(float DEta) { m_DEta = DEta; m_useDEta = true; }
  void SetDPhiIn(float DPhiIn) { m_DPhiIn = DPhiIn; m_useDPhiIn; }
  void SetDPhiOut(float DPhiOut) { m_DPhiOut = DPhiOut; m_useDPhiOut = true; }
  void SetInvEminusInvP(float InvEminusInvP) { m_InvEminusInvP = InvEminusInvP; m_useInvEminusInvP = true; }
  void SetBremFraction(float BremFraction) { m_BremFraction = BremFraction; m_useBremFraction = true; }
  void SetSigmaEtaEta(float SigmaEtaEta) { m_SigmaEtaEta = SigmaEtaEta; m_useSigmaEtaEta = true; }
  void SetSigmaPhiPhi(float SigmaPhiPhi) { m_SigmaPhiPhi = SigmaPhiPhi; m_useSigmaPhiPhi = true; }
  void SetEOverPout(float EOverPout) { m_EOverPout = EOverPout; m_useEOverPout = true; }
  void SetEOverPin(float EOverPin) { m_EOverPin = EOverPin; m_useEOverPin = true; }
  void SetElectronClass(int electronClass) { m_electronClass = electronClass; m_electronClassInitialised = true; }
  void SetLikelihood(float Likelihood) { m_Likelihood = Likelihood; m_useLikelihood = true; }
  void SetEcalFiducialRegion(int word) { m_fiducialflag = word; }
  void SetRecoFlag(int word) { m_recoflag = word; }
  //! do not apply electron ID if the electron is tracker driven
  void applyElectronIDOnPFlowElectrons(bool what) { m_applyIDOnPFlow = what; }
  //! set event by event output of egamma cut-based electron ID
  void SetEgammaCutBasedID(bool egammaCutBased) { m_egammaCutBased = egammaCutBased; m_egammaCutBasedInitialised = true; }
  //! get output of the selector (class dependent)
  bool output();
  //! get output of the selector (class independent)
  bool outputNoClass();
  //! display the electron efficiency
  void diplayEfficiencies();

private:

  bool m_useCrackElectrons;
  bool m_useHOverE;
  bool m_useS9S25;
  bool m_useDEta;
  bool m_useDPhiIn;
  bool m_useDPhiOut;
  bool m_useInvEminusInvP;
  bool m_useBremFraction;
  bool m_useSigmaEtaEta;
  bool m_useSigmaPhiPhi;
  bool m_useEOverPout;
  bool m_useEOverPin;
  bool m_useLikelihood;
  bool m_electronClassInitialised;
  bool m_egammaCutBasedInitialised;

  float m_HOverE;
  float m_S9S25;
  float m_DEta;
  float m_DPhiIn;
  float m_DPhiOut;
  float m_InvEminusInvP;
  float m_BremFraction;
  float m_SigmaEtaEta;
  float m_SigmaPhiPhi;
  float m_EOverPout;
  float m_EOverPin;
  float m_Likelihood;
  int m_electronClass;
  int m_fiducialflag;
  int m_recoflag;
  bool m_egammaCutBased;
  
  bool m_applyIDOnPFlow;

  bool m_classDep;

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;

  //! counters for the efficiencies display, based on electron candidates
  Counters m_electronCounter;

};

#endif
