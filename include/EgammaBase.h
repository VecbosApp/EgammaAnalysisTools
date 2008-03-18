//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  4 11:53:25 2008 by ROOT version 5.14/00f
// from TTree ntp1/ntp1
// found on file: default.root
//////////////////////////////////////////////////////////

#ifndef EgammaBase_h
#define EgammaBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Utils;

class EgammaBase {
public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           nMc;
   Float_t         pMc[101];   //[nMc]
   Float_t         massMc[101];   //[nMc]
   Float_t         thetaMc[101];   //[nMc]
   Float_t         etaMc[101];   //[nMc]
   Float_t         phiMc[101];   //[nMc]
   Int_t           idMc[101];   //[nMc]
   Int_t           mothMc[101];   //[nMc]
   Int_t           nDauMc[101];   //[nMc]
   Float_t         xMc[101];   //[nMc]
   Float_t         yMc[101];   //[nMc]
   Float_t         zMc[101];   //[nMc]
   Int_t           nTrg;
   Bool_t          firedTrg[90];   //[nTrg]
   Bool_t          evtPresel;
   Int_t           nEle;
   Int_t           chargeEle[50];   //[nEle]
   Float_t         energyEle[50];   //[nEle]
   Float_t         etEle[50];   //[nEle]
   Float_t         momentumEle[50];   //[nEle]
   Float_t         thetaEle[50];   //[nEle]
   Float_t         etaEle[50];   //[nEle]
   Float_t         phiEle[50];   //[nEle]
   Float_t         pxEle[50];   //[nEle]
   Float_t         pyEle[50];   //[nEle]
   Float_t         pzEle[50];   //[nEle]
   Float_t         vertexXEle[50];   //[nEle]
   Float_t         vertexYEle[50];   //[nEle]
   Float_t         vertexZEle[50];   //[nEle]
   Float_t         massEle[50];   //[nEle]
   Float_t         mtEle[50];   //[nEle]
   Int_t           pdgIdEle[50];   //[nEle]
   Int_t           nDauEle[50];   //[nEle]
   Int_t           d1IndexEle[50];   //[nEle]
   Int_t           d2IndexEle[50];   //[nEle]
   Int_t           d1pdgIdEle[50];   //[nEle]
   Int_t           d2pdgIdEle[50];   //[nEle]
   Float_t         ecalEle[50];   //[nEle]
   Int_t           nCluEle[50];   //[nEle]
   Int_t           nCryEle[50];   //[nEle]
   Float_t         e3x3Ele[50];   //[nEle]
   Float_t         e5x5Ele[50];   //[nEle]
   Float_t         eMaxEle[50];   //[nEle]
/*    Float_t         latEle[50];   //[nEle] */
/*    Float_t         phiLatEle[50];   //[nEle] */
/*    Float_t         etaLatEle[50];   //[nEle] */
   Float_t         erawEle[50];   //[nEle]
   Float_t         caloEtaEle[50];   //[nEle]
   Float_t         caloPhiEle[50];   //[nEle]
   Float_t         e2x2Ele[50];   //[nEle]
   Float_t         e2ndEle[50];   //[nEle]
   Float_t         s1s9Ele[50];   //[nEle]
   Float_t         s9s25Ele[50];   //[nEle]
   Float_t         covEtaEtaEle[50];   //[nEle]
   Float_t         covEtaPhiEle[50];   //[nEle]
   Float_t         covPhiPhiEle[50];   //[nEle]
/*    Float_t         a20Ele[50];   //[nEle] */
/*    Float_t         a42Ele[50];   //[nEle] */
   Float_t         pxAtOuterEle[50];   //[nEle]
   Float_t         pyAtOuterEle[50];   //[nEle]
   Float_t         pzAtOuterEle[50];   //[nEle]
   Float_t         xAtOuterEle[50];   //[nEle]
   Float_t         yAtOuterEle[50];   //[nEle]
   Float_t         zAtOuterEle[50];   //[nEle]
   Float_t         pxAtInnerEle[50];   //[nEle]
   Float_t         pyAtInnerEle[50];   //[nEle]
   Float_t         pzAtInnerEle[50];   //[nEle]
   Float_t         xAtInnerEle[50];   //[nEle]
   Float_t         yAtInnerEle[50];   //[nEle]
   Float_t         zAtInnerEle[50];   //[nEle]
   Float_t         eleTrackNormalizedChi2Ele[50];   //[nEle]
   Float_t         eleTrackDxyEle[50];   //[nEle]
   Float_t         eleTrackD0Ele[50];   //[nEle]
   Float_t         eleTrackDszEle[50];   //[nEle]
   Float_t         eleTrackDzEle[50];   //[nEle]
   Float_t         eleTrackDxyErrorEle[50];   //[nEle]
   Float_t         eleTrackD0ErrorEle[50];   //[nEle]
   Float_t         eleTrackDszErrorEle[50];   //[nEle]
   Float_t         eleTrackDzErrorEle[50];   //[nEle]
   Float_t         eleTrackValidHitsEle[50];   //[nEle]
   Float_t         eleTrackLostHitsEle[50];   //[nEle]
   Float_t         eleTrackVxEle[50];   //[nEle]
   Float_t         eleTrackVyEle[50];   //[nEle]
   Float_t         eleTrackVzEle[50];   //[nEle]
   Float_t         eleFullCorrEEle[50];   //[nEle]
   Float_t         eleCaloCorrEEle[50];   //[nEle]
   Float_t         eleNxtalCorrEEle[50];   //[nEle]
   Float_t         eleRawEEle[50];   //[nEle]
   Float_t         eleTrackerPEle[50];   //[nEle]
   Int_t           eleClassEle[50];   //[nEle]
   Float_t         eleHoEEle[50];   //[nEle]
   Float_t         eleCorrEoPEle[50];   //[nEle]
   Float_t         eleNotCorrEoPEle[50];   //[nEle]
   Float_t         eleCorrEoPoutEle[50];   //[nEle]
   Float_t         eleNotCorrEoPoutEle[50];   //[nEle]
   Float_t         eleDeltaEtaAtVtxEle[50];   //[nEle]
   Float_t         eleDeltaPhiAtVtxEle[50];   //[nEle]
   Float_t         eleDeltaEtaAtCaloEle[50];   //[nEle]
   Float_t         eleDeltaPhiAtCaloEle[50];   //[nEle]
   Float_t         eleTrackerIso_sumPtEle[50];   //[nEle]
   Float_t         eleCaloIso_sumPtEle[50];   //[nEle]
   Bool_t          eleIdCutBasedEle[50];   //[nEle]
   Float_t         eleLikelihoodEle[50];   //[nEle]
   Float_t         eleTipEle[50];   //[nEle]
   Int_t           nMuon;
   Int_t           chargeMuon[50];   //[nMuon]
   Float_t         energyMuon[50];   //[nMuon]
   Float_t         etMuon[50];   //[nMuon]
   Float_t         momentumMuon[50];   //[nMuon]
   Float_t         thetaMuon[50];   //[nMuon]
   Float_t         etaMuon[50];   //[nMuon]
   Float_t         phiMuon[50];   //[nMuon]
   Float_t         pxMuon[50];   //[nMuon]
   Float_t         pyMuon[50];   //[nMuon]
   Float_t         pzMuon[50];   //[nMuon]
   Float_t         vertexXMuon[50];   //[nMuon]
   Float_t         vertexYMuon[50];   //[nMuon]
   Float_t         vertexZMuon[50];   //[nMuon]
   Float_t         massMuon[50];   //[nMuon]
   Float_t         mtMuon[50];   //[nMuon]
   Int_t           pdgIdMuon[50];   //[nMuon]
   Int_t           nDauMuon[50];   //[nMuon]
   Int_t           d1IndexMuon[50];   //[nMuon]
   Int_t           d2IndexMuon[50];   //[nMuon]
   Int_t           d1pdgIdMuon[50];   //[nMuon]
   Int_t           d2pdgIdMuon[50];   //[nMuon]
   Int_t           nMet;
   Int_t           chargeMet[1];   //[nMet]
   Float_t         energyMet[1];   //[nMet]
   Float_t         etMet[1];   //[nMet]
   Float_t         momentumMet[1];   //[nMet]
   Float_t         thetaMet[1];   //[nMet]
   Float_t         etaMet[1];   //[nMet]
   Float_t         phiMet[1];   //[nMet]
   Float_t         pxMet[1];   //[nMet]
   Float_t         pyMet[1];   //[nMet]
   Float_t         pzMet[1];   //[nMet]
   Float_t         vertexXMet[1];   //[nMet]
   Float_t         vertexYMet[1];   //[nMet]
   Float_t         vertexZMet[1];   //[nMet]
   Float_t         massMet[1];   //[nMet]
   Float_t         mtMet[1];   //[nMet]
   Int_t           pdgIdMet[1];   //[nMet]
   Int_t           nDauMet[1];   //[nMet]
   Int_t           d1IndexMet[1];   //[nMet]
   Int_t           d2IndexMet[1];   //[nMet]
   Int_t           d1pdgIdMet[1];   //[nMet]
   Int_t           d2pdgIdMet[1];   //[nMet]
   Int_t           nGenMet;
   Int_t           chargeGenMet[1];   //[nGenMet]
   Float_t         energyGenMet[1];   //[nGenMet]
   Float_t         etGenMet[1];   //[nGenMet]
   Float_t         momentumGenMet[1];   //[nGenMet]
   Float_t         thetaGenMet[1];   //[nGenMet]
   Float_t         etaGenMet[1];   //[nGenMet]
   Float_t         phiGenMet[1];   //[nGenMet]
   Float_t         pxGenMet[1];   //[nGenMet]
   Float_t         pyGenMet[1];   //[nGenMet]
   Float_t         pzGenMet[1];   //[nGenMet]
   Float_t         vertexXGenMet[1];   //[nGenMet]
   Float_t         vertexYGenMet[1];   //[nGenMet]
   Float_t         vertexZGenMet[1];   //[nGenMet]
   Float_t         massGenMet[1];   //[nGenMet]
   Float_t         mtGenMet[1];   //[nGenMet]
   Int_t           pdgIdGenMet[1];   //[nGenMet]
   Int_t           nDauGenMet[1];   //[nGenMet]
   Int_t           d1IndexGenMet[1];   //[nGenMet]
   Int_t           d2IndexGenMet[1];   //[nGenMet]
   Int_t           d1pdgIdGenMet[1];   //[nGenMet]
   Int_t           d2pdgIdGenMet[1];   //[nGenMet]
   Int_t           nJet;
   Int_t           chargeJet[40];   //[nJet]
   Float_t         energyJet[40];   //[nJet]
   Float_t         etJet[40];   //[nJet]
   Float_t         momentumJet[40];   //[nJet]
   Float_t         thetaJet[40];   //[nJet]
   Float_t         etaJet[40];   //[nJet]
   Float_t         phiJet[40];   //[nJet]
   Float_t         pxJet[40];   //[nJet]
   Float_t         pyJet[40];   //[nJet]
   Float_t         pzJet[40];   //[nJet]
   Float_t         vertexXJet[40];   //[nJet]
   Float_t         vertexYJet[40];   //[nJet]
   Float_t         vertexZJet[40];   //[nJet]
   Float_t         massJet[40];   //[nJet]
   Float_t         mtJet[40];   //[nJet]
   Int_t           pdgIdJet[40];   //[nJet]
   Int_t           nDauJet[40];   //[nJet]
   Int_t           d1IndexJet[40];   //[nJet]
   Int_t           d2IndexJet[40];   //[nJet]
   Int_t           d1pdgIdJet[40];   //[nJet]
   Int_t           d2pdgIdJet[40];   //[nJet]
   Float_t         alphaJet[40];   //[nJet]
   Float_t         emFracJet[40];   //[nJet]
   Float_t         hadFracJet[40];   //[nJet]
   Int_t           nGenJet;
   Int_t           chargeGenJet[68];   //[nGenJet]
   Float_t         energyGenJet[68];   //[nGenJet]
   Float_t         etGenJet[68];   //[nGenJet]
   Float_t         momentumGenJet[68];   //[nGenJet]
   Float_t         thetaGenJet[68];   //[nGenJet]
   Float_t         etaGenJet[68];   //[nGenJet]
   Float_t         phiGenJet[68];   //[nGenJet]
   Float_t         pxGenJet[68];   //[nGenJet]
   Float_t         pyGenJet[68];   //[nGenJet]
   Float_t         pzGenJet[68];   //[nGenJet]
   Float_t         vertexXGenJet[68];   //[nGenJet]
   Float_t         vertexYGenJet[68];   //[nGenJet]
   Float_t         vertexZGenJet[68];   //[nGenJet]
   Float_t         massGenJet[68];   //[nGenJet]
   Float_t         mtGenJet[68];   //[nGenJet]
   Int_t           pdgIdGenJet[68];   //[nGenJet]
   Int_t           nDauGenJet[68];   //[nGenJet]
   Int_t           d1IndexGenJet[68];   //[nGenJet]
   Int_t           d2IndexGenJet[68];   //[nGenJet]
   Int_t           d1pdgIdGenJet[68];   //[nGenJet]
   Int_t           d2pdgIdGenJet[68];   //[nGenJet]

   // List of branches
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_massMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_nDauMc;   //!
   TBranch        *b_xMc;   //!
   TBranch        *b_yMc;   //!
   TBranch        *b_zMc;   //!
   TBranch        *b_nTrg;   //!
   TBranch        *b_firedTrg;   //!
   TBranch        *b_evtPresel;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_etEle;   //!
   TBranch        *b_momentumEle;   //!
   TBranch        *b_thetaEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_pxEle;   //!
   TBranch        *b_pyEle;   //!
   TBranch        *b_pzEle;   //!
   TBranch        *b_vertexXEle;   //!
   TBranch        *b_vertexYEle;   //!
   TBranch        *b_vertexZEle;   //!
   TBranch        *b_massEle;   //!
   TBranch        *b_mtEle;   //!
   TBranch        *b_pdgIdEle;   //!
   TBranch        *b_nDauEle;   //!
   TBranch        *b_d1IndexEle;   //!
   TBranch        *b_d2IndexEle;   //!
   TBranch        *b_d1pdgIdEle;   //!
   TBranch        *b_d2pdgIdEle;   //!
   TBranch        *b_ecalEle;   //!
   TBranch        *b_nCluEle;   //!
   TBranch        *b_nCryEle;   //!
   TBranch        *b_e3x3Ele;   //!
   TBranch        *b_e5x5Ele;   //!
   TBranch        *b_eMaxEle;   //!
/*    TBranch        *b_latEle;   //! */
/*    TBranch        *b_phiLatEle;   //! */
/*    TBranch        *b_etaLatEle;   //! */
   TBranch        *b_erawEle;   //!
   TBranch        *b_caloEtaEle;   //!
   TBranch        *b_caloPhiEle;   //!
   TBranch        *b_e2x2Ele;   //!
   TBranch        *b_e2ndEle;   //!
   TBranch        *b_s1s9Ele;   //!
   TBranch        *b_s9s25Ele;   //!
   TBranch        *b_covEtaEtaEle;   //!
   TBranch        *b_covEtaPhiEle;   //!
   TBranch        *b_covPhiPhiEle;   //!
/*    TBranch        *b_a20Ele;   //! */
/*    TBranch        *b_a42Ele;   //! */
   TBranch        *b_pxAtOuterEle;   //!
   TBranch        *b_pyAtOuterEle;   //!
   TBranch        *b_pzAtOuterEle;   //!
   TBranch        *b_xAtOuterEle;   //!
   TBranch        *b_yAtOuterEle;   //!
   TBranch        *b_zAtOuterEle;   //!
   TBranch        *b_pxAtInnerEle;   //!
   TBranch        *b_pyAtInnerEle;   //!
   TBranch        *b_pzAtInnerEle;   //!
   TBranch        *b_xAtInnerEle;   //!
   TBranch        *b_yAtInnerEle;   //!
   TBranch        *b_zAtInnerEle;   //!
   TBranch        *b_eleTrackNormalizedChi2Ele;   //!
   TBranch        *b_eleTrackDxyEle;   //!
   TBranch        *b_eleTrackD0Ele;   //!
   TBranch        *b_eleTrackDszEle;   //!
   TBranch        *b_eleTrackDzEle;   //!
   TBranch        *b_eleTrackDxyErrorEle;   //!
   TBranch        *b_eleTrackD0ErrorEle;   //!
   TBranch        *b_eleTrackDszErrorEle;   //!
   TBranch        *b_eleTrackDzErrorEle;   //!
   TBranch        *b_eleTrackValidHitsEle;   //!
   TBranch        *b_eleTrackLostHitsEle;   //!
   TBranch        *b_eleTrackVxEle;   //!
   TBranch        *b_eleTrackVyEle;   //!
   TBranch        *b_eleTrackVzEle;   //!
   TBranch        *b_eleFullCorrEEle;   //!
   TBranch        *b_eleCaloCorrEEle;   //!
   TBranch        *b_eleNxtalCorrEEle;   //!
   TBranch        *b_eleRawEEle;   //!
   TBranch        *b_eleTrackerPEle;   //!
   TBranch        *b_eleClassEle;   //!
   TBranch        *b_eleHoEEle;   //!
   TBranch        *b_eleCorrEoPEle;   //!
   TBranch        *b_eleNotCorrEoPEle;   //!
   TBranch        *b_eleCorrEoPoutEle;   //!
   TBranch        *b_eleNotCorrEoPoutEle;   //!
   TBranch        *b_eleDeltaEtaAtVtxEle;   //!
   TBranch        *b_eleDeltaPhiAtVtxEle;   //!
   TBranch        *b_eleDeltaEtaAtCaloEle;   //!
   TBranch        *b_eleDeltaPhiAtCaloEle;   //!
   TBranch        *b_eleTrackerIso_sumPtEle;   //!
   TBranch        *b_eleCaloIso_sumPtEle;   //!
   TBranch        *b_eleIdCutBasedEle;   //!
   TBranch        *b_eleLikelihoodEle;   //!
   TBranch        *b_eleTipEle;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_etMuon;   //!
   TBranch        *b_momentumMuon;   //!
   TBranch        *b_thetaMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_pxMuon;   //!
   TBranch        *b_pyMuon;   //!
   TBranch        *b_pzMuon;   //!
   TBranch        *b_vertexXMuon;   //!
   TBranch        *b_vertexYMuon;   //!
   TBranch        *b_vertexZMuon;   //!
   TBranch        *b_massMuon;   //!
   TBranch        *b_mtMuon;   //!
   TBranch        *b_pdgIdMuon;   //!
   TBranch        *b_nDauMuon;   //!
   TBranch        *b_d1IndexMuon;   //!
   TBranch        *b_d2IndexMuon;   //!
   TBranch        *b_d1pdgIdMuon;   //!
   TBranch        *b_d2pdgIdMuon;   //!
   TBranch        *b_nMet;   //!
   TBranch        *b_chargeMet;   //!
   TBranch        *b_energyMet;   //!
   TBranch        *b_etMet;   //!
   TBranch        *b_momentumMet;   //!
   TBranch        *b_thetaMet;   //!
   TBranch        *b_etaMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_pxMet;   //!
   TBranch        *b_pyMet;   //!
   TBranch        *b_pzMet;   //!
   TBranch        *b_vertexXMet;   //!
   TBranch        *b_vertexYMet;   //!
   TBranch        *b_vertexZMet;   //!
   TBranch        *b_massMet;   //!
   TBranch        *b_mtMet;   //!
   TBranch        *b_pdgIdMet;   //!
   TBranch        *b_nDauMet;   //!
   TBranch        *b_d1IndexMet;   //!
   TBranch        *b_d2IndexMet;   //!
   TBranch        *b_d1pdgIdMet;   //!
   TBranch        *b_d2pdgIdMet;   //!
   TBranch        *b_nGenMet;   //!
   TBranch        *b_chargeGenMet;   //!
   TBranch        *b_energyGenMet;   //!
   TBranch        *b_etGenMet;   //!
   TBranch        *b_momentumGenMet;   //!
   TBranch        *b_thetaGenMet;   //!
   TBranch        *b_etaGenMet;   //!
   TBranch        *b_phiGenMet;   //!
   TBranch        *b_pxGenMet;   //!
   TBranch        *b_pyGenMet;   //!
   TBranch        *b_pzGenMet;   //!
   TBranch        *b_vertexXGenMet;   //!
   TBranch        *b_vertexYGenMet;   //!
   TBranch        *b_vertexZGenMet;   //!
   TBranch        *b_massGenMet;   //!
   TBranch        *b_mtGenMet;   //!
   TBranch        *b_pdgIdGenMet;   //!
   TBranch        *b_nDauGenMet;   //!
   TBranch        *b_d1IndexGenMet;   //!
   TBranch        *b_d2IndexGenMet;   //!
   TBranch        *b_d1pdgIdGenMet;   //!
   TBranch        *b_d2pdgIdGenMet;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_chargeJet;   //!
   TBranch        *b_energyJet;   //!
   TBranch        *b_etJet;   //!
   TBranch        *b_momentumJet;   //!
   TBranch        *b_thetaJet;   //!
   TBranch        *b_etaJet;   //!
   TBranch        *b_phiJet;   //!
   TBranch        *b_pxJet;   //!
   TBranch        *b_pyJet;   //!
   TBranch        *b_pzJet;   //!
   TBranch        *b_vertexXJet;   //!
   TBranch        *b_vertexYJet;   //!
   TBranch        *b_vertexZJet;   //!
   TBranch        *b_massJet;   //!
   TBranch        *b_mtJet;   //!
   TBranch        *b_pdgIdJet;   //!
   TBranch        *b_nDauJet;   //!
   TBranch        *b_d1IndexJet;   //!
   TBranch        *b_d2IndexJet;   //!
   TBranch        *b_d1pdgIdJet;   //!
   TBranch        *b_d2pdgIdJet;   //!
   TBranch        *b_alphaJet;   //!
   TBranch        *b_emFracJet;   //!
   TBranch        *b_hadFracJet;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_chargeGenJet;   //!
   TBranch        *b_energyGenJet;   //!
   TBranch        *b_etGenJet;   //!
   TBranch        *b_momentumGenJet;   //!
   TBranch        *b_thetaGenJet;   //!
   TBranch        *b_etaGenJet;   //!
   TBranch        *b_phiGenJet;   //!
   TBranch        *b_pxGenJet;   //!
   TBranch        *b_pyGenJet;   //!
   TBranch        *b_pzGenJet;   //!
   TBranch        *b_vertexXGenJet;   //!
   TBranch        *b_vertexYGenJet;   //!
   TBranch        *b_vertexZGenJet;   //!
   TBranch        *b_massGenJet;   //!
   TBranch        *b_mtGenJet;   //!
   TBranch        *b_pdgIdGenJet;   //!
   TBranch        *b_nDauGenJet;   //!
   TBranch        *b_d1IndexGenJet;   //!
   TBranch        *b_d2IndexGenJet;   //!
   TBranch        *b_d1pdgIdGenJet;   //!
   TBranch        *b_d2pdgIdGenJet;   //!

   EgammaBase(TTree *tree=0);
   virtual ~EgammaBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EgammaBase_cxx
EgammaBase::EgammaBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("default.root");
      if (!f) {
         f = new TFile("default.root");
      }
      tree = (TTree*)gDirectory->Get("ntp1");

   }
   Init(tree);
}

EgammaBase::~EgammaBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EgammaBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EgammaBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EgammaBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
   fChain->SetBranchAddress("pMc", pMc, &b_pMc);
   fChain->SetBranchAddress("massMc", massMc, &b_massMc);
   fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
   fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
   fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("nDauMc", nDauMc, &b_nDauMc);
   fChain->SetBranchAddress("xMc", xMc, &b_xMc);
   fChain->SetBranchAddress("yMc", yMc, &b_yMc);
   fChain->SetBranchAddress("zMc", zMc, &b_zMc);
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("evtPresel", &evtPresel, &b_evtPresel);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("etEle", etEle, &b_etEle);
   fChain->SetBranchAddress("momentumEle", momentumEle, &b_momentumEle);
   fChain->SetBranchAddress("thetaEle", thetaEle, &b_thetaEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("pxEle", pxEle, &b_pxEle);
   fChain->SetBranchAddress("pyEle", pyEle, &b_pyEle);
   fChain->SetBranchAddress("pzEle", pzEle, &b_pzEle);
   fChain->SetBranchAddress("vertexXEle", vertexXEle, &b_vertexXEle);
   fChain->SetBranchAddress("vertexYEle", vertexYEle, &b_vertexYEle);
   fChain->SetBranchAddress("vertexZEle", vertexZEle, &b_vertexZEle);
   fChain->SetBranchAddress("massEle", massEle, &b_massEle);
   fChain->SetBranchAddress("mtEle", mtEle, &b_mtEle);
   fChain->SetBranchAddress("pdgIdEle", pdgIdEle, &b_pdgIdEle);
   fChain->SetBranchAddress("nDauEle", nDauEle, &b_nDauEle);
   fChain->SetBranchAddress("d1IndexEle", d1IndexEle, &b_d1IndexEle);
   fChain->SetBranchAddress("d2IndexEle", d2IndexEle, &b_d2IndexEle);
   fChain->SetBranchAddress("d1pdgIdEle", d1pdgIdEle, &b_d1pdgIdEle);
   fChain->SetBranchAddress("d2pdgIdEle", d2pdgIdEle, &b_d2pdgIdEle);
   fChain->SetBranchAddress("ecalEle", ecalEle, &b_ecalEle);
   fChain->SetBranchAddress("nCluEle", nCluEle, &b_nCluEle);
   fChain->SetBranchAddress("nCryEle", nCryEle, &b_nCryEle);
   fChain->SetBranchAddress("e3x3Ele", e3x3Ele, &b_e3x3Ele);
   fChain->SetBranchAddress("e5x5Ele", e5x5Ele, &b_e5x5Ele);
   fChain->SetBranchAddress("eMaxEle", eMaxEle, &b_eMaxEle);
/*    fChain->SetBranchAddress("latEle", latEle, &b_latEle); */
/*    fChain->SetBranchAddress("phiLatEle", phiLatEle, &b_phiLatEle); */
/*    fChain->SetBranchAddress("etaLatEle", etaLatEle, &b_etaLatEle); */
   fChain->SetBranchAddress("erawEle", erawEle, &b_erawEle);
   fChain->SetBranchAddress("caloEtaEle", caloEtaEle, &b_caloEtaEle);
   fChain->SetBranchAddress("caloPhiEle", caloPhiEle, &b_caloPhiEle);
   fChain->SetBranchAddress("e2x2Ele", e2x2Ele, &b_e2x2Ele);
   fChain->SetBranchAddress("e2ndEle", e2ndEle, &b_e2ndEle);
   fChain->SetBranchAddress("s1s9Ele", s1s9Ele, &b_s1s9Ele);
   fChain->SetBranchAddress("s9s25Ele", s9s25Ele, &b_s9s25Ele);
   fChain->SetBranchAddress("covEtaEtaEle", covEtaEtaEle, &b_covEtaEtaEle);
   fChain->SetBranchAddress("covEtaPhiEle", covEtaPhiEle, &b_covEtaPhiEle);
   fChain->SetBranchAddress("covPhiPhiEle", covPhiPhiEle, &b_covPhiPhiEle);
/*    fChain->SetBranchAddress("a20Ele", a20Ele, &b_a20Ele); */
/*    fChain->SetBranchAddress("a42Ele", a42Ele, &b_a42Ele); */
   fChain->SetBranchAddress("pxAtOuterEle", pxAtOuterEle, &b_pxAtOuterEle);
   fChain->SetBranchAddress("pyAtOuterEle", pyAtOuterEle, &b_pyAtOuterEle);
   fChain->SetBranchAddress("pzAtOuterEle", pzAtOuterEle, &b_pzAtOuterEle);
   fChain->SetBranchAddress("xAtOuterEle", xAtOuterEle, &b_xAtOuterEle);
   fChain->SetBranchAddress("yAtOuterEle", yAtOuterEle, &b_yAtOuterEle);
   fChain->SetBranchAddress("zAtOuterEle", zAtOuterEle, &b_zAtOuterEle);
   fChain->SetBranchAddress("pxAtInnerEle", pxAtInnerEle, &b_pxAtInnerEle);
   fChain->SetBranchAddress("pyAtInnerEle", pyAtInnerEle, &b_pyAtInnerEle);
   fChain->SetBranchAddress("pzAtInnerEle", pzAtInnerEle, &b_pzAtInnerEle);
   fChain->SetBranchAddress("xAtInnerEle", xAtInnerEle, &b_xAtInnerEle);
   fChain->SetBranchAddress("yAtInnerEle", yAtInnerEle, &b_yAtInnerEle);
   fChain->SetBranchAddress("zAtInnerEle", zAtInnerEle, &b_zAtInnerEle);
   fChain->SetBranchAddress("eleTrackNormalizedChi2Ele", eleTrackNormalizedChi2Ele, &b_eleTrackNormalizedChi2Ele);
   fChain->SetBranchAddress("eleTrackDxyEle", eleTrackDxyEle, &b_eleTrackDxyEle);
   fChain->SetBranchAddress("eleTrackD0Ele", eleTrackD0Ele, &b_eleTrackD0Ele);
   fChain->SetBranchAddress("eleTrackDszEle", eleTrackDszEle, &b_eleTrackDszEle);
   fChain->SetBranchAddress("eleTrackDzEle", eleTrackDzEle, &b_eleTrackDzEle);
   fChain->SetBranchAddress("eleTrackDxyErrorEle", eleTrackDxyErrorEle, &b_eleTrackDxyErrorEle);
   fChain->SetBranchAddress("eleTrackD0ErrorEle", eleTrackD0ErrorEle, &b_eleTrackD0ErrorEle);
   fChain->SetBranchAddress("eleTrackDszErrorEle", eleTrackDszErrorEle, &b_eleTrackDszErrorEle);
   fChain->SetBranchAddress("eleTrackDzErrorEle", eleTrackDzErrorEle, &b_eleTrackDzErrorEle);
   fChain->SetBranchAddress("eleTrackValidHitsEle", eleTrackValidHitsEle, &b_eleTrackValidHitsEle);
   fChain->SetBranchAddress("eleTrackLostHitsEle", eleTrackLostHitsEle, &b_eleTrackLostHitsEle);
   fChain->SetBranchAddress("eleTrackVxEle", eleTrackVxEle, &b_eleTrackVxEle);
   fChain->SetBranchAddress("eleTrackVyEle", eleTrackVyEle, &b_eleTrackVyEle);
   fChain->SetBranchAddress("eleTrackVzEle", eleTrackVzEle, &b_eleTrackVzEle);
   fChain->SetBranchAddress("eleFullCorrEEle", eleFullCorrEEle, &b_eleFullCorrEEle);
   fChain->SetBranchAddress("eleCaloCorrEEle", eleCaloCorrEEle, &b_eleCaloCorrEEle);
   fChain->SetBranchAddress("eleNxtalCorrEEle", eleNxtalCorrEEle, &b_eleNxtalCorrEEle);
   fChain->SetBranchAddress("eleRawEEle", eleRawEEle, &b_eleRawEEle);
   fChain->SetBranchAddress("eleTrackerPEle", eleTrackerPEle, &b_eleTrackerPEle);
   fChain->SetBranchAddress("eleClassEle", eleClassEle, &b_eleClassEle);
   fChain->SetBranchAddress("eleHoEEle", eleHoEEle, &b_eleHoEEle);
   fChain->SetBranchAddress("eleCorrEoPEle", eleCorrEoPEle, &b_eleCorrEoPEle);
   fChain->SetBranchAddress("eleNotCorrEoPEle", eleNotCorrEoPEle, &b_eleNotCorrEoPEle);
   fChain->SetBranchAddress("eleCorrEoPoutEle", eleCorrEoPoutEle, &b_eleCorrEoPoutEle);
   fChain->SetBranchAddress("eleNotCorrEoPoutEle", eleNotCorrEoPoutEle, &b_eleNotCorrEoPoutEle);
   fChain->SetBranchAddress("eleDeltaEtaAtVtxEle", eleDeltaEtaAtVtxEle, &b_eleDeltaEtaAtVtxEle);
   fChain->SetBranchAddress("eleDeltaPhiAtVtxEle", eleDeltaPhiAtVtxEle, &b_eleDeltaPhiAtVtxEle);
   fChain->SetBranchAddress("eleDeltaEtaAtCaloEle", eleDeltaEtaAtCaloEle, &b_eleDeltaEtaAtCaloEle);
   fChain->SetBranchAddress("eleDeltaPhiAtCaloEle", eleDeltaPhiAtCaloEle, &b_eleDeltaPhiAtCaloEle);
   fChain->SetBranchAddress("eleTrackerIso_sumPtEle", eleTrackerIso_sumPtEle, &b_eleTrackerIso_sumPtEle);
   fChain->SetBranchAddress("eleCaloIso_sumPtEle", eleCaloIso_sumPtEle, &b_eleCaloIso_sumPtEle);
   fChain->SetBranchAddress("eleIdCutBasedEle", eleIdCutBasedEle, &b_eleIdCutBasedEle);
   fChain->SetBranchAddress("eleLikelihoodEle", eleLikelihoodEle, &b_eleLikelihoodEle);
   fChain->SetBranchAddress("eleTipEle", eleTipEle, &b_eleTipEle);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("chargeMuon", &chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("energyMuon", &energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("etMuon", &etMuon, &b_etMuon);
   fChain->SetBranchAddress("momentumMuon", &momentumMuon, &b_momentumMuon);
   fChain->SetBranchAddress("thetaMuon", &thetaMuon, &b_thetaMuon);
   fChain->SetBranchAddress("etaMuon", &etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", &phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("pxMuon", &pxMuon, &b_pxMuon);
   fChain->SetBranchAddress("pyMuon", &pyMuon, &b_pyMuon);
   fChain->SetBranchAddress("pzMuon", &pzMuon, &b_pzMuon);
   fChain->SetBranchAddress("vertexXMuon", &vertexXMuon, &b_vertexXMuon);
   fChain->SetBranchAddress("vertexYMuon", &vertexYMuon, &b_vertexYMuon);
   fChain->SetBranchAddress("vertexZMuon", &vertexZMuon, &b_vertexZMuon);
   fChain->SetBranchAddress("massMuon", &massMuon, &b_massMuon);
   fChain->SetBranchAddress("mtMuon", &mtMuon, &b_mtMuon);
   fChain->SetBranchAddress("pdgIdMuon", &pdgIdMuon, &b_pdgIdMuon);
   fChain->SetBranchAddress("nDauMuon", &nDauMuon, &b_nDauMuon);
   fChain->SetBranchAddress("d1IndexMuon", &d1IndexMuon, &b_d1IndexMuon);
   fChain->SetBranchAddress("d2IndexMuon", &d2IndexMuon, &b_d2IndexMuon);
   fChain->SetBranchAddress("d1pdgIdMuon", &d1pdgIdMuon, &b_d1pdgIdMuon);
   fChain->SetBranchAddress("d2pdgIdMuon", &d2pdgIdMuon, &b_d2pdgIdMuon);
   fChain->SetBranchAddress("nMet", &nMet, &b_nMet);
   fChain->SetBranchAddress("chargeMet", chargeMet, &b_chargeMet);
   fChain->SetBranchAddress("energyMet", energyMet, &b_energyMet);
   fChain->SetBranchAddress("etMet", etMet, &b_etMet);
   fChain->SetBranchAddress("momentumMet", momentumMet, &b_momentumMet);
   fChain->SetBranchAddress("thetaMet", thetaMet, &b_thetaMet);
   fChain->SetBranchAddress("etaMet", etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMet", pxMet, &b_pxMet);
   fChain->SetBranchAddress("pyMet", pyMet, &b_pyMet);
   fChain->SetBranchAddress("pzMet", pzMet, &b_pzMet);
   fChain->SetBranchAddress("vertexXMet", vertexXMet, &b_vertexXMet);
   fChain->SetBranchAddress("vertexYMet", vertexYMet, &b_vertexYMet);
   fChain->SetBranchAddress("vertexZMet", vertexZMet, &b_vertexZMet);
   fChain->SetBranchAddress("massMet", massMet, &b_massMet);
   fChain->SetBranchAddress("mtMet", mtMet, &b_mtMet);
   fChain->SetBranchAddress("pdgIdMet", pdgIdMet, &b_pdgIdMet);
   fChain->SetBranchAddress("nDauMet", nDauMet, &b_nDauMet);
   fChain->SetBranchAddress("d1IndexMet", d1IndexMet, &b_d1IndexMet);
   fChain->SetBranchAddress("d2IndexMet", d2IndexMet, &b_d2IndexMet);
   fChain->SetBranchAddress("d1pdgIdMet", d1pdgIdMet, &b_d1pdgIdMet);
   fChain->SetBranchAddress("d2pdgIdMet", d2pdgIdMet, &b_d2pdgIdMet);
   fChain->SetBranchAddress("nGenMet", &nGenMet, &b_nGenMet);
   fChain->SetBranchAddress("chargeGenMet", chargeGenMet, &b_chargeGenMet);
   fChain->SetBranchAddress("energyGenMet", energyGenMet, &b_energyGenMet);
   fChain->SetBranchAddress("etGenMet", etGenMet, &b_etGenMet);
   fChain->SetBranchAddress("momentumGenMet", momentumGenMet, &b_momentumGenMet);
   fChain->SetBranchAddress("thetaGenMet", thetaGenMet, &b_thetaGenMet);
   fChain->SetBranchAddress("etaGenMet", etaGenMet, &b_etaGenMet);
   fChain->SetBranchAddress("phiGenMet", phiGenMet, &b_phiGenMet);
   fChain->SetBranchAddress("pxGenMet", pxGenMet, &b_pxGenMet);
   fChain->SetBranchAddress("pyGenMet", pyGenMet, &b_pyGenMet);
   fChain->SetBranchAddress("pzGenMet", pzGenMet, &b_pzGenMet);
   fChain->SetBranchAddress("vertexXGenMet", vertexXGenMet, &b_vertexXGenMet);
   fChain->SetBranchAddress("vertexYGenMet", vertexYGenMet, &b_vertexYGenMet);
   fChain->SetBranchAddress("vertexZGenMet", vertexZGenMet, &b_vertexZGenMet);
   fChain->SetBranchAddress("massGenMet", massGenMet, &b_massGenMet);
   fChain->SetBranchAddress("mtGenMet", mtGenMet, &b_mtGenMet);
   fChain->SetBranchAddress("pdgIdGenMet", pdgIdGenMet, &b_pdgIdGenMet);
   fChain->SetBranchAddress("nDauGenMet", nDauGenMet, &b_nDauGenMet);
   fChain->SetBranchAddress("d1IndexGenMet", d1IndexGenMet, &b_d1IndexGenMet);
   fChain->SetBranchAddress("d2IndexGenMet", d2IndexGenMet, &b_d2IndexGenMet);
   fChain->SetBranchAddress("d1pdgIdGenMet", d1pdgIdGenMet, &b_d1pdgIdGenMet);
   fChain->SetBranchAddress("d2pdgIdGenMet", d2pdgIdGenMet, &b_d2pdgIdGenMet);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("chargeJet", chargeJet, &b_chargeJet);
   fChain->SetBranchAddress("energyJet", energyJet, &b_energyJet);
   fChain->SetBranchAddress("etJet", etJet, &b_etJet);
   fChain->SetBranchAddress("momentumJet", momentumJet, &b_momentumJet);
   fChain->SetBranchAddress("thetaJet", thetaJet, &b_thetaJet);
   fChain->SetBranchAddress("etaJet", etaJet, &b_etaJet);
   fChain->SetBranchAddress("phiJet", phiJet, &b_phiJet);
   fChain->SetBranchAddress("pxJet", pxJet, &b_pxJet);
   fChain->SetBranchAddress("pyJet", pyJet, &b_pyJet);
   fChain->SetBranchAddress("pzJet", pzJet, &b_pzJet);
   fChain->SetBranchAddress("vertexXJet", vertexXJet, &b_vertexXJet);
   fChain->SetBranchAddress("vertexYJet", vertexYJet, &b_vertexYJet);
   fChain->SetBranchAddress("vertexZJet", vertexZJet, &b_vertexZJet);
   fChain->SetBranchAddress("massJet", massJet, &b_massJet);
   fChain->SetBranchAddress("mtJet", mtJet, &b_mtJet);
   fChain->SetBranchAddress("pdgIdJet", pdgIdJet, &b_pdgIdJet);
   fChain->SetBranchAddress("nDauJet", nDauJet, &b_nDauJet);
   fChain->SetBranchAddress("d1IndexJet", d1IndexJet, &b_d1IndexJet);
   fChain->SetBranchAddress("d2IndexJet", d2IndexJet, &b_d2IndexJet);
   fChain->SetBranchAddress("d1pdgIdJet", d1pdgIdJet, &b_d1pdgIdJet);
   fChain->SetBranchAddress("d2pdgIdJet", d2pdgIdJet, &b_d2pdgIdJet);
   fChain->SetBranchAddress("alphaJet", alphaJet, &b_alphaJet);
   fChain->SetBranchAddress("emFracJet", emFracJet, &b_emFracJet);
   fChain->SetBranchAddress("hadFracJet", hadFracJet, &b_hadFracJet);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("chargeGenJet", chargeGenJet, &b_chargeGenJet);
   fChain->SetBranchAddress("energyGenJet", energyGenJet, &b_energyGenJet);
   fChain->SetBranchAddress("etGenJet", etGenJet, &b_etGenJet);
   fChain->SetBranchAddress("momentumGenJet", momentumGenJet, &b_momentumGenJet);
   fChain->SetBranchAddress("thetaGenJet", thetaGenJet, &b_thetaGenJet);
   fChain->SetBranchAddress("etaGenJet", etaGenJet, &b_etaGenJet);
   fChain->SetBranchAddress("phiGenJet", phiGenJet, &b_phiGenJet);
   fChain->SetBranchAddress("pxGenJet", pxGenJet, &b_pxGenJet);
   fChain->SetBranchAddress("pyGenJet", pyGenJet, &b_pyGenJet);
   fChain->SetBranchAddress("pzGenJet", pzGenJet, &b_pzGenJet);
   fChain->SetBranchAddress("vertexXGenJet", vertexXGenJet, &b_vertexXGenJet);
   fChain->SetBranchAddress("vertexYGenJet", vertexYGenJet, &b_vertexYGenJet);
   fChain->SetBranchAddress("vertexZGenJet", vertexZGenJet, &b_vertexZGenJet);
   fChain->SetBranchAddress("massGenJet", massGenJet, &b_massGenJet);
   fChain->SetBranchAddress("mtGenJet", mtGenJet, &b_mtGenJet);
   fChain->SetBranchAddress("pdgIdGenJet", pdgIdGenJet, &b_pdgIdGenJet);
   fChain->SetBranchAddress("nDauGenJet", nDauGenJet, &b_nDauGenJet);
   fChain->SetBranchAddress("d1IndexGenJet", d1IndexGenJet, &b_d1IndexGenJet);
   fChain->SetBranchAddress("d2IndexGenJet", d2IndexGenJet, &b_d2IndexGenJet);
   fChain->SetBranchAddress("d1pdgIdGenJet", d1pdgIdGenJet, &b_d1pdgIdGenJet);
   fChain->SetBranchAddress("d2pdgIdGenJet", d2pdgIdGenJet, &b_d2pdgIdGenJet);
   Notify();
}

Bool_t EgammaBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EgammaBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EgammaBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EgammaBase_cxx
