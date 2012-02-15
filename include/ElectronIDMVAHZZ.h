//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVAHZZ
//
// Helper Class for applying MVA electron ID selection
//
// Authors: D.Benedetti, impl. E.Di Marco
//--------------------------------------------------------------------------------------------------

#ifndef HIGGSANALYSIS_ElectronIDMVAHZZ_H
#define HIGGSANALYSIS_ElectronIDMVAHZZ_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class ElectronIDMVAHZZ {
  public:
    ElectronIDMVAHZZ();
    ~ElectronIDMVAHZZ(); 

    enum MVAType {
      kBDTSimpleCat = 0,      // the BDT used in H->ZZ
      kBDTSimpleCatData      // the BDT used in H->ZZ, trained on DATA
    };

    void   Initialize(std::string methodName,
                      std::string weights , 
                      ElectronIDMVAHZZ::MVAType type );
    Bool_t IsInitialized() const { return fIsInitialized; }
    
    double MVAValue(double ElePt , double EleSCEta,
                    double EleFBrem,
                    double EleDEtaIn,
                    double EleDPhiIn,
                    double EleDEtaEleOut,
                    double EleSigmaIEtaIEta,
                    double EleHoverE,
                    double EleSuperClusterEOverP,
                    double EleE1x5E5x5,
                    double EleEOverPout,
                    double EleKFChi2,
                    double EleKFHits,
                    double EleMissHits,
                    double EleDistConv,
                    double EleDcotConv,
                    double NVtx,
                    double EleEcalSeeded,
                    double EleEtaWidth,
                    double ElePhiWidth,
                    double EleD0,
                    double EleIP3d,
                    double EleIP3dSig);



  protected:
    TMVA::Reader             *fTMVAReader[1];
    std::string               fMethodname;
    MVAType                   fMVAType;
    
    Bool_t                    fIsInitialized;
    Float_t                   fMVAVar_EleSCEta;
    Float_t                   fMVAVar_ElePt;
    Float_t                   fMVAVar_EleFBrem;     
    Float_t                   fMVAVar_EleDEtaIn; 
    Float_t                   fMVAVar_EleDPhiIn; 
    Float_t                   fMVAVar_EleDEtaEleOut;
    Float_t                   fMVAVar_EleSigmaIEtaIEta; 
    Float_t                   fMVAVar_EleHoverE; 
    Float_t                   fMVAVar_EleSuperClusterEOverP;
    Float_t                   fMVAVar_EleE1x5E5x5;
    Float_t                   fMVAVar_EleEOverPout;
    Float_t                   fMVAVar_EleKFChi2;
    Float_t                   fMVAVar_EleKFHits;
    Float_t                   fMVAVar_EleMissHits;
    Float_t                   fMVAVar_EleDistConv;
    Float_t                   fMVAVar_EleDcotConv;
    Float_t                   fMVAVar_NVtx;
    Float_t                   fMVAVar_EleEcalSeeded;
    Float_t                   fMVAVar_EleEtaWidth;
    Float_t                   fMVAVar_ElePhiWidth;
    Float_t                   fMVAVar_EleD0;
    Float_t                   fMVAVar_EleIP3d;
    Float_t                   fMVAVar_EleIP3dSig;
    
};

#endif
