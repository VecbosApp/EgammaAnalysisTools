#include <TFile.h>
#include "EgammaAnalysisTools/include/ElectronIDMVAHZZ.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <math.h>
#include <assert.h>

//--------------------------------------------------------------------------------------------------
ElectronIDMVAHZZ::ElectronIDMVAHZZ() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<1; ++i) {
    fTMVAReader[i] = 0;
  }
}



//--------------------------------------------------------------------------------------------------
ElectronIDMVAHZZ::~ElectronIDMVAHZZ()
{
  for(UInt_t i=0; i<1; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVAHZZ::Initialize( std::string methodName,
                                   std::string weights , 
                                   ElectronIDMVAHZZ::MVAType type) {

  fIsInitialized = kTRUE;
  fMVAType = type;

  fMethodname = methodName;
    
  for(UInt_t i=0; i<1; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kBDTSimpleCat) {

      fTMVAReader[i]->AddVariable("fbrem",       &fMVAVar_EleFBrem);
      fTMVAReader[i]->AddVariable("detain",      &fMVAVar_EleDEtaIn);
      fTMVAReader[i]->AddVariable("dphiin",      &fMVAVar_EleDPhiIn);
      fTMVAReader[i]->AddVariable("sieie",       &fMVAVar_EleSigmaIEtaIEta);
      fTMVAReader[i]->AddVariable("hoe",         &fMVAVar_EleHoverE);
      fTMVAReader[i]->AddVariable("eop",         &fMVAVar_EleSuperClusterEOverP);
      fTMVAReader[i]->AddVariable("e1x5e5x5",    &fMVAVar_EleE1x5E5x5);
      fTMVAReader[i]->AddVariable("eleopout",    &fMVAVar_EleEOverPout);
      fTMVAReader[i]->AddVariable("detaeleout",  &fMVAVar_EleDEtaEleOut);
      fTMVAReader[i]->AddVariable("kfchi2",      &fMVAVar_EleKFChi2);
      fTMVAReader[i]->AddVariable("kfhits",      &fMVAVar_EleKFHits);
      fTMVAReader[i]->AddVariable("mishits",     &fMVAVar_EleMissHits);
      fTMVAReader[i]->AddVariable("dist",        &fMVAVar_EleDistConv);
      fTMVAReader[i]->AddVariable("dcot",        &fMVAVar_EleDcotConv);
      fTMVAReader[i]->AddVariable("nvtx",        &fMVAVar_NVtx);  // for the new weight file

      fTMVAReader[i]->AddSpectator("eta",        &fMVAVar_EleSCEta);
      fTMVAReader[i]->AddSpectator("pt",         &fMVAVar_ElePt);
      fTMVAReader[i]->AddSpectator("ecalseed",   &fMVAVar_EleEcalSeeded);

    }
    
    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , weights );

  }

  std::cout << "Electron ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << weights << std::endl;

}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVAHZZ::MVAValue(Double_t ElePt , Double_t EleSCEta,
                                    Double_t EleFBrem,
                                    Double_t EleDEtaIn,
                                    Double_t EleDPhiIn,
                                    Double_t EleDEtaEleOut,
                                    Double_t EleSigmaIEtaIEta,
                                    Double_t EleHoverE,
                                    double EleSuperClusterEOverP,
                                    Double_t EleE1x5E5x5,
                                    Double_t EleEOverPout,
                                    Double_t EleKFChi2,
                                    Double_t EleKFHits,
                                    Double_t EleMissHits,
                                    Double_t EleDistConv,
                                    Double_t EleDcotConv,
                                    Double_t NVtx,
                                    Double_t EleEcalSeeded
                                    ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVAHZZ not properly initialized.\n"; 
    return -9999;
  }

  //set all input variables
  fMVAVar_EleSCEta = EleSCEta;
  fMVAVar_ElePt = ElePt;
  fMVAVar_EleFBrem = EleFBrem;     
  fMVAVar_EleDEtaIn = EleDEtaIn; 
  fMVAVar_EleDPhiIn = EleDPhiIn; 
  fMVAVar_EleDEtaEleOut = EleDEtaEleOut;
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta; 
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleSuperClusterEOverP = EleSuperClusterEOverP;
  fMVAVar_EleE1x5E5x5 = EleE1x5E5x5;
  fMVAVar_EleEOverPout = EleEOverPout;
  fMVAVar_EleKFChi2 = EleKFChi2;
  fMVAVar_EleKFHits = EleMissHits;
  fMVAVar_EleMissHits = EleMissHits;
  fMVAVar_EleDistConv = EleDistConv;
  fMVAVar_EleDcotConv = EleDcotConv;
  fMVAVar_NVtx = NVtx;
  fMVAVar_EleEcalSeeded = EleEcalSeeded;

  // apply the boundaries used in the training
  if(fMVAVar_EleFBrem < -1.)
    fMVAVar_EleFBrem = -1.;

  fMVAVar_EleDEtaIn = fabs(fMVAVar_EleDEtaIn);
  if(fMVAVar_EleDEtaIn > 0.06)
    fMVAVar_EleDEtaIn = 0.06;

  fMVAVar_EleDPhiIn = fabs(fMVAVar_EleDPhiIn);
  if(fMVAVar_EleDPhiIn > 0.6)
    fMVAVar_EleDPhiIn = 0.6;

  if(fMVAVar_EleEOverPout > 20.)
    fMVAVar_EleEOverPout = 20;

  if(fMVAVar_EleSuperClusterEOverP > 20.)
    fMVAVar_EleSuperClusterEOverP = 20.;

  fMVAVar_EleDEtaEleOut = fabs(fMVAVar_EleDEtaEleOut);
  if(fMVAVar_EleDEtaEleOut > 0.2)
    fMVAVar_EleDEtaEleOut = 0.2;

  if(fMVAVar_EleKFChi2 < 0.)
    fMVAVar_EleKFChi2 = 0.;

  if(fMVAVar_EleKFChi2 > 15.)
    fMVAVar_EleKFChi2 = 15.;


  if(fMVAVar_EleE1x5E5x5 < -1.)
    fMVAVar_EleE1x5E5x5 = -1;

  if(fMVAVar_EleE1x5E5x5 > 2.)
    fMVAVar_EleE1x5E5x5 = 2.;

  if(fMVAVar_EleDistConv > 15.)
    fMVAVar_EleDistConv = 15.;
  if(fMVAVar_EleDistConv < -15.)
    fMVAVar_EleDistConv = -15.;
  fMVAVar_EleDistConv = fabs(fMVAVar_EleDistConv);

  if(fMVAVar_EleDcotConv > 3.)
    fMVAVar_EleDcotConv = 3.;
  if(fMVAVar_EleDcotConv < -3.)
    fMVAVar_EleDcotConv = -3.;
  fMVAVar_EleDcotConv = fabs(fMVAVar_EleDcotConv);

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = 0;
  // assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
  
                                                
  mva = reader->EvaluateMVA( fMethodname );

//   std::cout << "********************************\n";
//   std::cout << "Electron MVA\n";
//   std::cout << ElePt << " no eleeta " << " non ho phi " <<
//     " : " << EleSCEta << " : " << MVABin << std::endl;
//   std::cout << fMVAVar_EleFBrem << endl     
//   << fMVAVar_EleDEtaIn << endl 
//   << fMVAVar_EleDPhiIn << endl 
//   << fMVAVar_EleDEtaEleOut << endl
//   << fMVAVar_EleSigmaIEtaIEta << endl 
//   << fMVAVar_EleHoverE << endl 
//   << fMVAVar_EleE1x5E5x5 << endl
//   << fMVAVar_EleEOverPout << endl
//   << fMVAVar_EleKFChi2 << endl
//   << fMVAVar_EleKFHits << endl
//   << fMVAVar_EleMissHits << endl
//   << fMVAVar_EleDistConv << endl
//   << fMVAVar_EleDcotConv << endl
//   << fMVAVar_NVtx << endl
//   << fMVAVar_EleEcalSeeded << endl
//   << std::endl;
//   std::cout << "MVA: " << mva << std::endl;
//   std::cout << "********************************\n";

  return mva;
}
