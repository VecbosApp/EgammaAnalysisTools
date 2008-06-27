#ifndef ELECTRONCALOISOLATION_H
#define ELECTRONCALOISOLATION_H

#include <vector>
#include "TLorentzVector.h"

class ElectronCaloIsolation {
 public:

  //! constructors
  ElectronCaloIsolation();
  ElectronCaloIsolation( TLorentzVector electronSuperCluster );

  //! destructor
  ~ElectronCaloIsolation() {}

  //! HCAL/ECAL towers of the event
  void setTowers( std::vector<TLorentzVector> towers) { m_towers = towers; }

  //! internal and external cones
  void setIntRadius (float intRadius) { m_intRadius = intRadius; }
  void setExtRadius (float extRadius) { m_extRadius = extRadius; }

  //! if true, returns sumPt / ptEle
  void setRelative(bool relative) { m_relative = relative; }
  //! if relative, it is necessary to set the electron Pt at vertex
  void setElectronMomentumAtVtx(TLorentzVector electronAtVertex) { m_electronAtVertex = electronAtVertex; }

  //! get the Et sums in a cone. If relative => give sumEt/pT electron
  float getSumEt();

 private:
  
  TLorentzVector m_electronSuperCluster;
  TLorentzVector m_electronAtVertex;
  std::vector<TLorentzVector> m_towers;
  
  float m_extRadius;
  float m_intRadius;

  bool m_relative;

};

#endif

