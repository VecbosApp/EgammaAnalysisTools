#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include <iostream>
#include <math.h>

CutBasedEleIDSelector::CutBasedEleIDSelector() {

  // if the variable value is not initialized, it is not used
  m_useHOverE = false;
  m_useS9S25 = false;
  m_useDEta = false;
  m_useDPhiIn = false;
  m_useDPhiOut = false;
  m_useInvEminusInvP = false;
  m_useBremFraction = false;
  m_useSigmaEtaEta = false;
  m_useSigmaPhiPhi = false;
  m_useEOverPout = false;
  m_useEOverPin = false;
  m_useLikelihood = false;
  m_electronClassInitialised = false;
  m_egammaCutBasedInitialised = false;

}

CutBasedEleIDSelector::~CutBasedEleIDSelector() {}

void CutBasedEleIDSelector::Configure(const char *configDir) {

  std::string fileSwitches = std::string(configDir) + "electronIDSwitches.txt";

  std::string goldenCutsEB = std::string(configDir) + "electronIDGoldenCutsEB.txt";
  std::string bigbremCutsEB = std::string(configDir) + "electronIDBigBremCutsEB.txt";
  std::string narrowCutsEB = std::string(configDir) + "electronIDNarrowCutsEB.txt";
  std::string showeringCutsEB = std::string(configDir) + "electronIDShoweringCutsEB.txt";

  std::string goldenCutsEE = std::string(configDir) + "electronIDGoldenCutsEE.txt";
  std::string bigbremCutsEE = std::string(configDir) + "electronIDBigBremCutsEE.txt";
  std::string narrowCutsEE = std::string(configDir) + "electronIDNarrowCutsEE.txt";
  std::string showeringCutsEE = std::string(configDir) + "electronIDShoweringCutsEE.txt";
  
  Selection *goldenSelectionEB = new Selection(goldenCutsEB, fileSwitches);
  Selection *bigbremSelectionEB = new Selection(bigbremCutsEB, fileSwitches);
  Selection *narrowSelectionEB = new Selection(narrowCutsEB, fileSwitches);
  Selection *showeringSelectionEB = new Selection(showeringCutsEB, fileSwitches);

  Selection *goldenSelectionEE = new Selection(goldenCutsEE, fileSwitches);
  Selection *bigbremSelectionEE = new Selection(bigbremCutsEE, fileSwitches);
  Selection *narrowSelectionEE = new Selection(narrowCutsEE, fileSwitches);
  Selection *showeringSelectionEE = new Selection(showeringCutsEE, fileSwitches);

  m_EgammaCutIDSelection.push_back(goldenSelectionEB);
  m_EgammaCutIDSelection.push_back(bigbremSelectionEB);
  m_EgammaCutIDSelection.push_back(narrowSelectionEB);
  m_EgammaCutIDSelection.push_back(showeringSelectionEB);
  m_EgammaCutIDSelection.push_back(goldenSelectionEE);
  m_EgammaCutIDSelection.push_back(bigbremSelectionEE);
  m_EgammaCutIDSelection.push_back(narrowSelectionEE);
  m_EgammaCutIDSelection.push_back(showeringSelectionEE);

  std::vector<Selection*>::iterator selectionItr;
  for(selectionItr=m_EgammaCutIDSelection.begin();selectionItr!=m_EgammaCutIDSelection.end();selectionItr++) {
    Selection *eleSelection = *selectionItr;
    eleSelection->addSwitch("egammaCutBased");
    eleSelection->addCut("hOverE");
    eleSelection->addCut("s9s25");
    eleSelection->addCut("deta");
    eleSelection->addCut("dphiIn");
    eleSelection->addCut("dphiOut");
    eleSelection->addCut("invEminusInvP");
    eleSelection->addCut("bremFraction"); 
    eleSelection->addCut("sigmaEtaEta");
    eleSelection->addCut("sigmaPhiPhi");
    eleSelection->addCut("eOverPout");
    eleSelection->addCut("eOverPin");
    eleSelection->addCut("likelihood");
    eleSelection->summary();
  }


  m_electronCounter.SetTitle("SINGLE ELECTRON COUNTER");
  m_electronCounter.AddVar("electrons");
  m_electronCounter.AddVar("egammaCutBased");
  m_electronCounter.AddVar("hOverE");
  m_electronCounter.AddVar("s9s25");
  m_electronCounter.AddVar("deta");
  m_electronCounter.AddVar("dphiIn");
  m_electronCounter.AddVar("dphiOut");
  m_electronCounter.AddVar("invEminusInvP");
  m_electronCounter.AddVar("bremFraction");
  m_electronCounter.AddVar("sigmaEtaEta");
  m_electronCounter.AddVar("sigmaPhiPhi");
  m_electronCounter.AddVar("eOverPout");
  m_electronCounter.AddVar("eOverPin");
  m_electronCounter.AddVar("finalCustomEleID");
  m_electronCounter.AddVar("likelihood");


}



bool CutBasedEleIDSelector::output() {

  int GsfClass = m_electronClass;
  int offset=0;
  if(GsfClass>=100) {
    GsfClass-=100;
    offset=4;
  }

  Selection *selection;
  if(GsfClass==0) 
    selection=m_EgammaCutIDSelection[offset];
  else if(GsfClass==10)
    selection=m_EgammaCutIDSelection[offset+1];
  else if(GsfClass==20)
    selection=m_EgammaCutIDSelection[offset+2];
  else if(GsfClass>=30 && GsfClass<=40) {
    selection=m_EgammaCutIDSelection[offset+3];
  }

  m_electronCounter.IncrVar("electrons");

  if(selection->getSwitch("egammaCutBased") && 
     m_egammaCutBased ) m_electronCounter.IncrVar("egammaCutBased");

  if(selection->getSwitch("likelihood") && 
     selection->passCut("likelihood", fabs(m_Likelihood) )) {  
    m_electronCounter.IncrVar("likelihood");
  }

  if(selection->getSwitch("hOverE") && 
     !selection->passCut("hOverE", fabs(m_HOverE) )) return false;  
  m_electronCounter.IncrVar("hOverE");
  
  if(selection->getSwitch("s9s25") && 
     !selection->passCut("s9s25", m_S9S25)) return false; 
  m_electronCounter.IncrVar("s9s25");
  
  if(selection->getSwitch("deta") && 
     !selection->passCut("deta", fabs(m_DEta) )) return false; 
  m_electronCounter.IncrVar("deta");

  if(selection->getSwitch("dphiIn") && 
     !selection->passCut("dphiIn", fabs(m_DPhiIn) )) return false; 
  m_electronCounter.IncrVar("dphiIn");

  if(selection->getSwitch("dphiOut") && 
     !selection->passCut("dphiOut", fabs(m_DPhiOut) )) return false; 
  m_electronCounter.IncrVar("dphiOut");
    
  if(selection->getSwitch("invEminusInvP") && 
     !selection->passCut("invEminusInvP", m_InvEminusInvP)) return false; 
  m_electronCounter.IncrVar("invEminusInvP");
    
  if(selection->getSwitch("bremFraction") && 
     !selection->passCut("bremFraction", m_BremFraction)) return false; 
  m_electronCounter.IncrVar("bremFraction");
  
  if(selection->getSwitch("sigmaEtaEta") && 
     !selection->passCut("sigmaEtaEta", m_SigmaEtaEta)) return false; 
  m_electronCounter.IncrVar("sigmaEtaEta");

  if(selection->getSwitch("sigmaPhiPhi") && 
     !selection->passCut("sigmaPhiPhi", m_SigmaEtaEta)) return false; 
  m_electronCounter.IncrVar("sigmaPhiPhi");

  if(selection->getSwitch("eOverPout") && 
     !selection->passCut("eOverPout", m_EOverPout)) return false; 
  m_electronCounter.IncrVar("eOverPout");

  if(selection->getSwitch("eOverPin") && 
     !selection->passCut("eOverPin", m_EOverPin)) return false; 
  m_electronCounter.IncrVar("eOverPin");
  
  m_electronCounter.IncrVar("finalCustomEleID");

  return true;


}


void CutBasedEleIDSelector::diplayEfficiencies() {

  m_electronCounter.Draw();
  m_electronCounter.Draw("egammaCutBased", "electrons");
  m_electronCounter.Draw("hOverE","electrons");
  m_electronCounter.Draw("s9s25","hOverE");
  m_electronCounter.Draw("deta","s9s25");
  m_electronCounter.Draw("dphiIn","deta");
  m_electronCounter.Draw("dphiOut","dphiIn");
  m_electronCounter.Draw("invEminusInvP", "dphiOut");
  m_electronCounter.Draw("bremFraction", "invEminusInvP");
  m_electronCounter.Draw("sigmaEtaEta","bremFraction");
  m_electronCounter.Draw("sigmaPhiPhi","sigmaEtaEta");
  m_electronCounter.Draw("eOverPout", "sigmaPhiPhi");
  m_electronCounter.Draw("eOverPin", "eOverPout");
  m_electronCounter.Draw("finalCustomEleID","electrons");
  m_electronCounter.Draw("likelihood","electrons");

}
