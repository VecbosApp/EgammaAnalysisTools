#include <iostream>
#include "EgammaAnalysisTools/include/ElectronCaloIsolation.hh"

ElectronCaloIsolation::ElectronCaloIsolation (){}

ElectronCaloIsolation::ElectronCaloIsolation (TLorentzVector electronSuperCluster) {
  m_electronSuperCluster = electronSuperCluster;

  //! default configuration
  m_extRadius = 0.40; // deltaR
  m_intRadius = 0.00; // deltaR
  m_relative = true;
  m_electronAtVertex.SetPxPyPzE(0,0,0,0);
}

float ElectronCaloIsolation::getSumEt() {

  float dummyPt = 0 ;

  float eleSC_pt  = m_electronSuperCluster.Pt();
  float ele_ptAtVtx = m_electronAtVertex.Pt();

  std::vector<TLorentzVector>::const_iterator tower;
  for(tower=m_towers.begin(); tower!=m_towers.end(); tower++) { 
    
    float tower_pt = tower->Pt();

    double dr = m_electronSuperCluster.DeltaR(*tower);
    if ( fabs(dr) < m_extRadius && fabs(dr) > m_intRadius ){ 
      if ( m_relative && ele_ptAtVtx == 0 ) {
	std::cout << "relative calo isolation required. But electron momentum at vertex not set."
		  << "\nPlease use setElectronMomentumAtVtx(TLorentzVector electronAtVertex)!"
		  << "\nreturning sumEt=10000" << std::endl;
	return 10000.;
      }
      else if ( m_relative && ele_ptAtVtx!=0 ) {
	dummyPt += tower_pt/ele_ptAtVtx;
      }
      else {
	dummyPt += tower_pt;
      }
    } 
    
  } //end loop over towers		       
  
  return dummyPt;
}
