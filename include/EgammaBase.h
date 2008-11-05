//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov  5 11:08:38 2008 by ROOT version 5.18/00a
// from TTree ntp1/ntp1
// found on file: default.root
//////////////////////////////////////////////////////////

#ifndef EgammaBase_h
#define EgammaBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class EgammaBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nMc;
   Float_t         pMc[101];   //[nMc]
   Float_t         massMc[101];   //[nMc]
   Float_t         thetaMc[101];   //[nMc]
   Float_t         etaMc[101];   //[nMc]
   Float_t         phiMc[101];   //[nMc]
   Float_t         energyMc[101];   //[nMc]
   Int_t           idMc[101];   //[nMc]
   Int_t           mothMc[101];   //[nMc]
   Int_t           nDauMc[101];   //[nMc]
   Int_t           statusMc[101];   //[nMc]
   Float_t         xMc[101];   //[nMc]
   Float_t         yMc[101];   //[nMc]
   Float_t         zMc[101];   //[nMc]
   Int_t           nTrg;
   UChar_t         firedTrg[160];   //[nTrg]
   UChar_t         evtPresel;
   Double_t        evtKfactor;
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
   Float_t         latEle[50];   //[nEle]
   Float_t         phiLatEle[50];   //[nEle]
   Float_t         etaLatEle[50];   //[nEle]
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
   Float_t         a20Ele[50];   //[nEle]
   Float_t         a42Ele[50];   //[nEle]
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
   Float_t         pxAtInnerEle[50];   //[nEle]
   Float_t         pyAtInnerEle[50];   //[nEle]
   Float_t         pzAtInnerEle[50];   //[nEle]
   Float_t         xAtInnerEle[50];   //[nEle]
   Float_t         yAtInnerEle[50];   //[nEle]
   Float_t         zAtInnerEle[50];   //[nEle]
   Float_t         pxAtOuterEle[50];   //[nEle]
   Float_t         pyAtOuterEle[50];   //[nEle]
   Float_t         pzAtOuterEle[50];   //[nEle]
   Float_t         xAtOuterEle[50];   //[nEle]
   Float_t         yAtOuterEle[50];   //[nEle]
   Float_t         zAtOuterEle[50];   //[nEle]
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
   Float_t         eleMinDR03Ele[50];   //[nEle]
   Float_t         eleMinDRveto03Ele[50];   //[nEle]
   Float_t         eleSumPt03Ele[50];   //[nEle]
   Float_t         eleSumPtSquared03Ele[50];   //[nEle]
   Float_t         eleSumN03Ele[50];   //[nEle]
   Float_t         eleSumPt04Ele[50];   //[nEle]
   Float_t         eleSumPt05Ele[50];   //[nEle]
   Float_t         eleSumPtPreselectionEle[50];   //[nEle]
   Float_t         eleSumHadEt04Ele[50];   //[nEle]
   Float_t         eleSumEmEt04Ele[50];   //[nEle]
   Float_t         eleSumHadEt05Ele[50];   //[nEle]
   Float_t         eleSumEmEt05Ele[50];   //[nEle]
   Float_t         eleIsoFromDepsTkEle[50];   //[nEle]
   Float_t         eleIsoFromDepsEcalEle[50];   //[nEle]
   Float_t         eleIsoFromDepsHcalEle[50];   //[nEle]
   UChar_t         eleIdCutBasedEle[50];   //[nEle]
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
   Int_t           nPFMet;
   Int_t           chargePFMet[1];   //[nPFMet]
   Float_t         energyPFMet[1];   //[nPFMet]
   Float_t         etPFMet[1];   //[nPFMet]
   Float_t         momentumPFMet[1];   //[nPFMet]
   Float_t         thetaPFMet[1];   //[nPFMet]
   Float_t         etaPFMet[1];   //[nPFMet]
   Float_t         phiPFMet[1];   //[nPFMet]
   Float_t         pxPFMet[1];   //[nPFMet]
   Float_t         pyPFMet[1];   //[nPFMet]
   Float_t         pzPFMet[1];   //[nPFMet]
   Float_t         vertexXPFMet[1];   //[nPFMet]
   Float_t         vertexYPFMet[1];   //[nPFMet]
   Float_t         vertexZPFMet[1];   //[nPFMet]
   Float_t         massPFMet[1];   //[nPFMet]
   Float_t         mtPFMet[1];   //[nPFMet]
   Int_t           pdgIdPFMet[1];   //[nPFMet]
   Int_t           nDauPFMet[1];   //[nPFMet]
   Int_t           d1IndexPFMet[1];   //[nPFMet]
   Int_t           d2IndexPFMet[1];   //[nPFMet]
   Int_t           d1pdgIdPFMet[1];   //[nPFMet]
   Int_t           d2pdgIdPFMet[1];   //[nPFMet]
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
   Int_t           nIterativeJet;
   Int_t           chargeIterativeJet[100];   //[nIterativeJet]
   Float_t         energyIterativeJet[100];   //[nIterativeJet]
   Float_t         etIterativeJet[100];   //[nIterativeJet]
   Float_t         momentumIterativeJet[100];   //[nIterativeJet]
   Float_t         thetaIterativeJet[100];   //[nIterativeJet]
   Float_t         etaIterativeJet[100];   //[nIterativeJet]
   Float_t         phiIterativeJet[100];   //[nIterativeJet]
   Float_t         pxIterativeJet[100];   //[nIterativeJet]
   Float_t         pyIterativeJet[100];   //[nIterativeJet]
   Float_t         pzIterativeJet[100];   //[nIterativeJet]
   Float_t         vertexXIterativeJet[100];   //[nIterativeJet]
   Float_t         vertexYIterativeJet[100];   //[nIterativeJet]
   Float_t         vertexZIterativeJet[100];   //[nIterativeJet]
   Float_t         massIterativeJet[100];   //[nIterativeJet]
   Float_t         mtIterativeJet[100];   //[nIterativeJet]
   Int_t           pdgIdIterativeJet[100];   //[nIterativeJet]
   Int_t           nDauIterativeJet[100];   //[nIterativeJet]
   Int_t           d1IndexIterativeJet[100];   //[nIterativeJet]
   Int_t           d2IndexIterativeJet[100];   //[nIterativeJet]
   Int_t           d1pdgIdIterativeJet[100];   //[nIterativeJet]
   Int_t           d2pdgIdIterativeJet[100];   //[nIterativeJet]
   Float_t         alphaIterativeJet[100];   //[nIterativeJet]
   Float_t         emFracIterativeJet[100];   //[nIterativeJet]
   Float_t         hadFracIterativeJet[100];   //[nIterativeJet]
   Int_t           nSisConeJet;
   Int_t           chargeSisConeJet[100];   //[nSisConeJet]
   Float_t         energySisConeJet[100];   //[nSisConeJet]
   Float_t         etSisConeJet[100];   //[nSisConeJet]
   Float_t         momentumSisConeJet[100];   //[nSisConeJet]
   Float_t         thetaSisConeJet[100];   //[nSisConeJet]
   Float_t         etaSisConeJet[100];   //[nSisConeJet]
   Float_t         phiSisConeJet[100];   //[nSisConeJet]
   Float_t         pxSisConeJet[100];   //[nSisConeJet]
   Float_t         pySisConeJet[100];   //[nSisConeJet]
   Float_t         pzSisConeJet[100];   //[nSisConeJet]
   Float_t         vertexXSisConeJet[100];   //[nSisConeJet]
   Float_t         vertexYSisConeJet[100];   //[nSisConeJet]
   Float_t         vertexZSisConeJet[100];   //[nSisConeJet]
   Float_t         massSisConeJet[100];   //[nSisConeJet]
   Float_t         mtSisConeJet[100];   //[nSisConeJet]
   Int_t           pdgIdSisConeJet[100];   //[nSisConeJet]
   Int_t           nDauSisConeJet[100];   //[nSisConeJet]
   Int_t           d1IndexSisConeJet[100];   //[nSisConeJet]
   Int_t           d2IndexSisConeJet[100];   //[nSisConeJet]
   Int_t           d1pdgIdSisConeJet[100];   //[nSisConeJet]
   Int_t           d2pdgIdSisConeJet[100];   //[nSisConeJet]
   Float_t         alphaSisConeJet[100];   //[nSisConeJet]
   Float_t         emFracSisConeJet[100];   //[nSisConeJet]
   Float_t         hadFracSisConeJet[100];   //[nSisConeJet]
   Int_t           nIterativePFJet;
   Int_t           chargeIterativePFJet[100];   //[nIterativePFJet]
   Float_t         energyIterativePFJet[100];   //[nIterativePFJet]
   Float_t         etIterativePFJet[100];   //[nIterativePFJet]
   Float_t         momentumIterativePFJet[100];   //[nIterativePFJet]
   Float_t         thetaIterativePFJet[100];   //[nIterativePFJet]
   Float_t         etaIterativePFJet[100];   //[nIterativePFJet]
   Float_t         phiIterativePFJet[100];   //[nIterativePFJet]
   Float_t         pxIterativePFJet[100];   //[nIterativePFJet]
   Float_t         pyIterativePFJet[100];   //[nIterativePFJet]
   Float_t         pzIterativePFJet[100];   //[nIterativePFJet]
   Float_t         vertexXIterativePFJet[100];   //[nIterativePFJet]
   Float_t         vertexYIterativePFJet[100];   //[nIterativePFJet]
   Float_t         vertexZIterativePFJet[100];   //[nIterativePFJet]
   Float_t         massIterativePFJet[100];   //[nIterativePFJet]
   Float_t         mtIterativePFJet[100];   //[nIterativePFJet]
   Int_t           pdgIdIterativePFJet[100];   //[nIterativePFJet]
   Int_t           nDauIterativePFJet[100];   //[nIterativePFJet]
   Int_t           d1IndexIterativePFJet[100];   //[nIterativePFJet]
   Int_t           d2IndexIterativePFJet[100];   //[nIterativePFJet]
   Int_t           d1pdgIdIterativePFJet[100];   //[nIterativePFJet]
   Int_t           d2pdgIdIterativePFJet[100];   //[nIterativePFJet]
   Float_t         chargedHadronEnergyIterativePFJet[100];   //[nIterativePFJet]
   Float_t         neutralHadronEnergyIterativePFJet[100];   //[nIterativePFJet]
   Float_t         chargedEmEnergyIterativePFJet[100];   //[nIterativePFJet]
   Float_t         neutralEmEnergyIterativePFJet[100];   //[nIterativePFJet]
   Float_t         neutralMultiplicityIterativePFJet[100];   //[nIterativePFJet]
   Float_t         chargedMultiplicityIterativePFJet[100];   //[nIterativePFJet]
   Float_t         muonMultiplicityIterativePFJet[100];   //[nIterativePFJet]
   Int_t           nSisConePFJet;
   Int_t           chargeSisConePFJet[100];   //[nSisConePFJet]
   Float_t         energySisConePFJet[100];   //[nSisConePFJet]
   Float_t         etSisConePFJet[100];   //[nSisConePFJet]
   Float_t         momentumSisConePFJet[100];   //[nSisConePFJet]
   Float_t         thetaSisConePFJet[100];   //[nSisConePFJet]
   Float_t         etaSisConePFJet[100];   //[nSisConePFJet]
   Float_t         phiSisConePFJet[100];   //[nSisConePFJet]
   Float_t         pxSisConePFJet[100];   //[nSisConePFJet]
   Float_t         pySisConePFJet[100];   //[nSisConePFJet]
   Float_t         pzSisConePFJet[100];   //[nSisConePFJet]
   Float_t         vertexXSisConePFJet[100];   //[nSisConePFJet]
   Float_t         vertexYSisConePFJet[100];   //[nSisConePFJet]
   Float_t         vertexZSisConePFJet[100];   //[nSisConePFJet]
   Float_t         massSisConePFJet[100];   //[nSisConePFJet]
   Float_t         mtSisConePFJet[100];   //[nSisConePFJet]
   Int_t           pdgIdSisConePFJet[100];   //[nSisConePFJet]
   Int_t           nDauSisConePFJet[100];   //[nSisConePFJet]
   Int_t           d1IndexSisConePFJet[100];   //[nSisConePFJet]
   Int_t           d2IndexSisConePFJet[100];   //[nSisConePFJet]
   Int_t           d1pdgIdSisConePFJet[100];   //[nSisConePFJet]
   Int_t           d2pdgIdSisConePFJet[100];   //[nSisConePFJet]
   Float_t         chargedHadronEnergySisConePFJet[100];   //[nSisConePFJet]
   Float_t         neutralHadronEnergySisConePFJet[100];   //[nSisConePFJet]
   Float_t         chargedEmEnergySisConePFJet[100];   //[nSisConePFJet]
   Float_t         neutralEmEnergySisConePFJet[100];   //[nSisConePFJet]
   Float_t         neutralMultiplicitySisConePFJet[100];   //[nSisConePFJet]
   Float_t         chargedMultiplicitySisConePFJet[100];   //[nSisConePFJet]
   Float_t         muonMultiplicitySisConePFJet[100];   //[nSisConePFJet]
   Int_t           nIterativeGenJet;
   Int_t           chargeIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         energyIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         etIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         momentumIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         thetaIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         etaIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         phiIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         pxIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         pyIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         pzIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         vertexXIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         vertexYIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         vertexZIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         massIterativeGenJet[100];   //[nIterativeGenJet]
   Float_t         mtIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           pdgIdIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           nDauIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           d1IndexIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           d2IndexIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           d1pdgIdIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           d2pdgIdIterativeGenJet[100];   //[nIterativeGenJet]
   Int_t           nSisConeGenJet;
   Int_t           chargeSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         energySisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         etSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         momentumSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         thetaSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         etaSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         phiSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         pxSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         pySisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         pzSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         vertexXSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         vertexYSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         vertexZSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         massSisConeGenJet[100];   //[nSisConeGenJet]
   Float_t         mtSisConeGenJet[100];   //[nSisConeGenJet]
   Int_t           pdgIdSisConeGenJet[100];   //[nSisConeGenJet]
   Int_t           nDauSisConeGenJet[100];   //[nSisConeGenJet]
   Int_t           d1IndexSisConeGenJet[100];   //[nSisConeGenJet]
   Int_t           d2IndexSisConeGenJet[100];   //[nSisConeGenJet]
   Int_t           d1pdgIdSisConeGenJet[100];   //[nSisConeGenJet]
   Int_t           d2pdgIdSisConeGenJet[100];   //[nSisConeGenJet]

   // List of branches
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_massMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_energyMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_nDauMc;   //!
   TBranch        *b_statusMc;   //!
   TBranch        *b_xMc;   //!
   TBranch        *b_yMc;   //!
   TBranch        *b_zMc;   //!
   TBranch        *b_nTrg;   //!
   TBranch        *b_firedTrg;   //!
   TBranch        *b_evtPresel;   //!
   TBranch        *b_evtKfactor;   //!
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
   TBranch        *b_latEle;   //!
   TBranch        *b_phiLatEle;   //!
   TBranch        *b_etaLatEle;   //!
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
   TBranch        *b_a20Ele;   //!
   TBranch        *b_a42Ele;   //!
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
   TBranch        *b_pxAtInnerEle;   //!
   TBranch        *b_pyAtInnerEle;   //!
   TBranch        *b_pzAtInnerEle;   //!
   TBranch        *b_xAtInnerEle;   //!
   TBranch        *b_yAtInnerEle;   //!
   TBranch        *b_zAtInnerEle;   //!
   TBranch        *b_pxAtOuterEle;   //!
   TBranch        *b_pyAtOuterEle;   //!
   TBranch        *b_pzAtOuterEle;   //!
   TBranch        *b_xAtOuterEle;   //!
   TBranch        *b_yAtOuterEle;   //!
   TBranch        *b_zAtOuterEle;   //!
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
   TBranch        *b_eleMinDR03Ele;   //!
   TBranch        *b_eleMinDRveto03Ele;   //!
   TBranch        *b_eleSumPt03Ele;   //!
   TBranch        *b_eleSumPtSquared03Ele;   //!
   TBranch        *b_eleSumN03Ele;   //!
   TBranch        *b_eleSumPt04Ele;   //!
   TBranch        *b_eleSumPt05Ele;   //!
   TBranch        *b_eleSumPtPreselectionEle;   //!
   TBranch        *b_eleSumHadEt04Ele;   //!
   TBranch        *b_eleSumEmEt04Ele;   //!
   TBranch        *b_eleSumHadEt05Ele;   //!
   TBranch        *b_eleSumEmEt05Ele;   //!
   TBranch        *b_eleIsoFromDepsTkEle;   //!
   TBranch        *b_eleIsoFromDepsEcalEle;   //!
   TBranch        *b_eleIsoFromDepsHcalEle;   //!
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
   TBranch        *b_nPFMet;   //!
   TBranch        *b_chargePFMet;   //!
   TBranch        *b_energyPFMet;   //!
   TBranch        *b_etPFMet;   //!
   TBranch        *b_momentumPFMet;   //!
   TBranch        *b_thetaPFMet;   //!
   TBranch        *b_etaPFMet;   //!
   TBranch        *b_phiPFMet;   //!
   TBranch        *b_pxPFMet;   //!
   TBranch        *b_pyPFMet;   //!
   TBranch        *b_pzPFMet;   //!
   TBranch        *b_vertexXPFMet;   //!
   TBranch        *b_vertexYPFMet;   //!
   TBranch        *b_vertexZPFMet;   //!
   TBranch        *b_massPFMet;   //!
   TBranch        *b_mtPFMet;   //!
   TBranch        *b_pdgIdPFMet;   //!
   TBranch        *b_nDauPFMet;   //!
   TBranch        *b_d1IndexPFMet;   //!
   TBranch        *b_d2IndexPFMet;   //!
   TBranch        *b_d1pdgIdPFMet;   //!
   TBranch        *b_d2pdgIdPFMet;   //!
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
   TBranch        *b_nIterativeJet;   //!
   TBranch        *b_chargeIterativeJet;   //!
   TBranch        *b_energyIterativeJet;   //!
   TBranch        *b_etIterativeJet;   //!
   TBranch        *b_momentumIterativeJet;   //!
   TBranch        *b_thetaIterativeJet;   //!
   TBranch        *b_etaIterativeJet;   //!
   TBranch        *b_phiIterativeJet;   //!
   TBranch        *b_pxIterativeJet;   //!
   TBranch        *b_pyIterativeJet;   //!
   TBranch        *b_pzIterativeJet;   //!
   TBranch        *b_vertexXIterativeJet;   //!
   TBranch        *b_vertexYIterativeJet;   //!
   TBranch        *b_vertexZIterativeJet;   //!
   TBranch        *b_massIterativeJet;   //!
   TBranch        *b_mtIterativeJet;   //!
   TBranch        *b_pdgIdIterativeJet;   //!
   TBranch        *b_nDauIterativeJet;   //!
   TBranch        *b_d1IndexIterativeJet;   //!
   TBranch        *b_d2IndexIterativeJet;   //!
   TBranch        *b_d1pdgIdIterativeJet;   //!
   TBranch        *b_d2pdgIdIterativeJet;   //!
   TBranch        *b_alphaIterativeJet;   //!
   TBranch        *b_emFracIterativeJet;   //!
   TBranch        *b_hadFracIterativeJet;   //!
   TBranch        *b_nSisConeJet;   //!
   TBranch        *b_chargeSisConeJet;   //!
   TBranch        *b_energySisConeJet;   //!
   TBranch        *b_etSisConeJet;   //!
   TBranch        *b_momentumSisConeJet;   //!
   TBranch        *b_thetaSisConeJet;   //!
   TBranch        *b_etaSisConeJet;   //!
   TBranch        *b_phiSisConeJet;   //!
   TBranch        *b_pxSisConeJet;   //!
   TBranch        *b_pySisConeJet;   //!
   TBranch        *b_pzSisConeJet;   //!
   TBranch        *b_vertexXSisConeJet;   //!
   TBranch        *b_vertexYSisConeJet;   //!
   TBranch        *b_vertexZSisConeJet;   //!
   TBranch        *b_massSisConeJet;   //!
   TBranch        *b_mtSisConeJet;   //!
   TBranch        *b_pdgIdSisConeJet;   //!
   TBranch        *b_nDauSisConeJet;   //!
   TBranch        *b_d1IndexSisConeJet;   //!
   TBranch        *b_d2IndexSisConeJet;   //!
   TBranch        *b_d1pdgIdSisConeJet;   //!
   TBranch        *b_d2pdgIdSisConeJet;   //!
   TBranch        *b_alphaSisConeJet;   //!
   TBranch        *b_emFracSisConeJet;   //!
   TBranch        *b_hadFracSisConeJet;   //!
   TBranch        *b_nIterativePFJet;   //!
   TBranch        *b_chargeIterativePFJet;   //!
   TBranch        *b_energyIterativePFJet;   //!
   TBranch        *b_etIterativePFJet;   //!
   TBranch        *b_momentumIterativePFJet;   //!
   TBranch        *b_thetaIterativePFJet;   //!
   TBranch        *b_etaIterativePFJet;   //!
   TBranch        *b_phiIterativePFJet;   //!
   TBranch        *b_pxIterativePFJet;   //!
   TBranch        *b_pyIterativePFJet;   //!
   TBranch        *b_pzIterativePFJet;   //!
   TBranch        *b_vertexXIterativePFJet;   //!
   TBranch        *b_vertexYIterativePFJet;   //!
   TBranch        *b_vertexZIterativePFJet;   //!
   TBranch        *b_massIterativePFJet;   //!
   TBranch        *b_mtIterativePFJet;   //!
   TBranch        *b_pdgIdIterativePFJet;   //!
   TBranch        *b_nDauIterativePFJet;   //!
   TBranch        *b_d1IndexIterativePFJet;   //!
   TBranch        *b_d2IndexIterativePFJet;   //!
   TBranch        *b_d1pdgIdIterativePFJet;   //!
   TBranch        *b_d2pdgIdIterativePFJet;   //!
   TBranch        *b_chargedHadronEnergyIterativePFJet;   //!
   TBranch        *b_neutralHadronEnergyIterativePFJet;   //!
   TBranch        *b_chargedEmEnergyIterativePFJet;   //!
   TBranch        *b_neutralEmEnergyIterativePFJet;   //!
   TBranch        *b_neutralMultiplicityIterativePFJet;   //!
   TBranch        *b_chargedMultiplicityIterativePFJet;   //!
   TBranch        *b_muonMultiplicityIterativePFJet;   //!
   TBranch        *b_nSisConePFJet;   //!
   TBranch        *b_chargeSisConePFJet;   //!
   TBranch        *b_energySisConePFJet;   //!
   TBranch        *b_etSisConePFJet;   //!
   TBranch        *b_momentumSisConePFJet;   //!
   TBranch        *b_thetaSisConePFJet;   //!
   TBranch        *b_etaSisConePFJet;   //!
   TBranch        *b_phiSisConePFJet;   //!
   TBranch        *b_pxSisConePFJet;   //!
   TBranch        *b_pySisConePFJet;   //!
   TBranch        *b_pzSisConePFJet;   //!
   TBranch        *b_vertexXSisConePFJet;   //!
   TBranch        *b_vertexYSisConePFJet;   //!
   TBranch        *b_vertexZSisConePFJet;   //!
   TBranch        *b_massSisConePFJet;   //!
   TBranch        *b_mtSisConePFJet;   //!
   TBranch        *b_pdgIdSisConePFJet;   //!
   TBranch        *b_nDauSisConePFJet;   //!
   TBranch        *b_d1IndexSisConePFJet;   //!
   TBranch        *b_d2IndexSisConePFJet;   //!
   TBranch        *b_d1pdgIdSisConePFJet;   //!
   TBranch        *b_d2pdgIdSisConePFJet;   //!
   TBranch        *b_chargedHadronEnergySisConePFJet;   //!
   TBranch        *b_neutralHadronEnergySisConePFJet;   //!
   TBranch        *b_chargedEmEnergySisConePFJet;   //!
   TBranch        *b_neutralEmEnergySisConePFJet;   //!
   TBranch        *b_neutralMultiplicitySisConePFJet;   //!
   TBranch        *b_chargedMultiplicitySisConePFJet;   //!
   TBranch        *b_muonMultiplicitySisConePFJet;   //!
   TBranch        *b_nIterativeGenJet;   //!
   TBranch        *b_chargeIterativeGenJet;   //!
   TBranch        *b_energyIterativeGenJet;   //!
   TBranch        *b_etIterativeGenJet;   //!
   TBranch        *b_momentumIterativeGenJet;   //!
   TBranch        *b_thetaIterativeGenJet;   //!
   TBranch        *b_etaIterativeGenJet;   //!
   TBranch        *b_phiIterativeGenJet;   //!
   TBranch        *b_pxIterativeGenJet;   //!
   TBranch        *b_pyIterativeGenJet;   //!
   TBranch        *b_pzIterativeGenJet;   //!
   TBranch        *b_vertexXIterativeGenJet;   //!
   TBranch        *b_vertexYIterativeGenJet;   //!
   TBranch        *b_vertexZIterativeGenJet;   //!
   TBranch        *b_massIterativeGenJet;   //!
   TBranch        *b_mtIterativeGenJet;   //!
   TBranch        *b_pdgIdIterativeGenJet;   //!
   TBranch        *b_nDauIterativeGenJet;   //!
   TBranch        *b_d1IndexIterativeGenJet;   //!
   TBranch        *b_d2IndexIterativeGenJet;   //!
   TBranch        *b_d1pdgIdIterativeGenJet;   //!
   TBranch        *b_d2pdgIdIterativeGenJet;   //!
   TBranch        *b_nSisConeGenJet;   //!
   TBranch        *b_chargeSisConeGenJet;   //!
   TBranch        *b_energySisConeGenJet;   //!
   TBranch        *b_etSisConeGenJet;   //!
   TBranch        *b_momentumSisConeGenJet;   //!
   TBranch        *b_thetaSisConeGenJet;   //!
   TBranch        *b_etaSisConeGenJet;   //!
   TBranch        *b_phiSisConeGenJet;   //!
   TBranch        *b_pxSisConeGenJet;   //!
   TBranch        *b_pySisConeGenJet;   //!
   TBranch        *b_pzSisConeGenJet;   //!
   TBranch        *b_vertexXSisConeGenJet;   //!
   TBranch        *b_vertexYSisConeGenJet;   //!
   TBranch        *b_vertexZSisConeGenJet;   //!
   TBranch        *b_massSisConeGenJet;   //!
   TBranch        *b_mtSisConeGenJet;   //!
   TBranch        *b_pdgIdSisConeGenJet;   //!
   TBranch        *b_nDauSisConeGenJet;   //!
   TBranch        *b_d1IndexSisConeGenJet;   //!
   TBranch        *b_d2IndexSisConeGenJet;   //!
   TBranch        *b_d1pdgIdSisConeGenJet;   //!
   TBranch        *b_d2pdgIdSisConeGenJet;   //!

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
   // It is normally not necessary to make changes to the generated
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
   fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("nDauMc", nDauMc, &b_nDauMc);
   fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   fChain->SetBranchAddress("xMc", xMc, &b_xMc);
   fChain->SetBranchAddress("yMc", yMc, &b_yMc);
   fChain->SetBranchAddress("zMc", zMc, &b_zMc);
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("evtPresel", &evtPresel, &b_evtPresel);
   fChain->SetBranchAddress("evtKfactor", &evtKfactor, &b_evtKfactor);
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
   fChain->SetBranchAddress("latEle", latEle, &b_latEle);
   fChain->SetBranchAddress("phiLatEle", phiLatEle, &b_phiLatEle);
   fChain->SetBranchAddress("etaLatEle", etaLatEle, &b_etaLatEle);
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
   fChain->SetBranchAddress("a20Ele", a20Ele, &b_a20Ele);
   fChain->SetBranchAddress("a42Ele", a42Ele, &b_a42Ele);
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
   fChain->SetBranchAddress("pxAtInnerEle", pxAtInnerEle, &b_pxAtInnerEle);
   fChain->SetBranchAddress("pyAtInnerEle", pyAtInnerEle, &b_pyAtInnerEle);
   fChain->SetBranchAddress("pzAtInnerEle", pzAtInnerEle, &b_pzAtInnerEle);
   fChain->SetBranchAddress("xAtInnerEle", xAtInnerEle, &b_xAtInnerEle);
   fChain->SetBranchAddress("yAtInnerEle", yAtInnerEle, &b_yAtInnerEle);
   fChain->SetBranchAddress("zAtInnerEle", zAtInnerEle, &b_zAtInnerEle);
   fChain->SetBranchAddress("pxAtOuterEle", pxAtOuterEle, &b_pxAtOuterEle);
   fChain->SetBranchAddress("pyAtOuterEle", pyAtOuterEle, &b_pyAtOuterEle);
   fChain->SetBranchAddress("pzAtOuterEle", pzAtOuterEle, &b_pzAtOuterEle);
   fChain->SetBranchAddress("xAtOuterEle", xAtOuterEle, &b_xAtOuterEle);
   fChain->SetBranchAddress("yAtOuterEle", yAtOuterEle, &b_yAtOuterEle);
   fChain->SetBranchAddress("zAtOuterEle", zAtOuterEle, &b_zAtOuterEle);
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
   fChain->SetBranchAddress("eleMinDR03Ele", eleMinDR03Ele, &b_eleMinDR03Ele);
   fChain->SetBranchAddress("eleMinDRveto03Ele", eleMinDRveto03Ele, &b_eleMinDRveto03Ele);
   fChain->SetBranchAddress("eleSumPt03Ele", eleSumPt03Ele, &b_eleSumPt03Ele);
   fChain->SetBranchAddress("eleSumPtSquared03Ele", eleSumPtSquared03Ele, &b_eleSumPtSquared03Ele);
   fChain->SetBranchAddress("eleSumN03Ele", eleSumN03Ele, &b_eleSumN03Ele);
   fChain->SetBranchAddress("eleSumPt04Ele", eleSumPt04Ele, &b_eleSumPt04Ele);
   fChain->SetBranchAddress("eleSumPt05Ele", eleSumPt05Ele, &b_eleSumPt05Ele);
   fChain->SetBranchAddress("eleSumPtPreselectionEle", eleSumPtPreselectionEle, &b_eleSumPtPreselectionEle);
   fChain->SetBranchAddress("eleSumHadEt04Ele", eleSumHadEt04Ele, &b_eleSumHadEt04Ele);
   fChain->SetBranchAddress("eleSumEmEt04Ele", eleSumEmEt04Ele, &b_eleSumEmEt04Ele);
   fChain->SetBranchAddress("eleSumHadEt05Ele", eleSumHadEt05Ele, &b_eleSumHadEt05Ele);
   fChain->SetBranchAddress("eleSumEmEt05Ele", eleSumEmEt05Ele, &b_eleSumEmEt05Ele);
   fChain->SetBranchAddress("eleIsoFromDepsTkEle", eleIsoFromDepsTkEle, &b_eleIsoFromDepsTkEle);
   fChain->SetBranchAddress("eleIsoFromDepsEcalEle", eleIsoFromDepsEcalEle, &b_eleIsoFromDepsEcalEle);
   fChain->SetBranchAddress("eleIsoFromDepsHcalEle", eleIsoFromDepsHcalEle, &b_eleIsoFromDepsHcalEle);
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
   fChain->SetBranchAddress("nPFMet", &nPFMet, &b_nPFMet);
   fChain->SetBranchAddress("chargePFMet", chargePFMet, &b_chargePFMet);
   fChain->SetBranchAddress("energyPFMet", energyPFMet, &b_energyPFMet);
   fChain->SetBranchAddress("etPFMet", etPFMet, &b_etPFMet);
   fChain->SetBranchAddress("momentumPFMet", momentumPFMet, &b_momentumPFMet);
   fChain->SetBranchAddress("thetaPFMet", thetaPFMet, &b_thetaPFMet);
   fChain->SetBranchAddress("etaPFMet", etaPFMet, &b_etaPFMet);
   fChain->SetBranchAddress("phiPFMet", phiPFMet, &b_phiPFMet);
   fChain->SetBranchAddress("pxPFMet", pxPFMet, &b_pxPFMet);
   fChain->SetBranchAddress("pyPFMet", pyPFMet, &b_pyPFMet);
   fChain->SetBranchAddress("pzPFMet", pzPFMet, &b_pzPFMet);
   fChain->SetBranchAddress("vertexXPFMet", vertexXPFMet, &b_vertexXPFMet);
   fChain->SetBranchAddress("vertexYPFMet", vertexYPFMet, &b_vertexYPFMet);
   fChain->SetBranchAddress("vertexZPFMet", vertexZPFMet, &b_vertexZPFMet);
   fChain->SetBranchAddress("massPFMet", massPFMet, &b_massPFMet);
   fChain->SetBranchAddress("mtPFMet", mtPFMet, &b_mtPFMet);
   fChain->SetBranchAddress("pdgIdPFMet", pdgIdPFMet, &b_pdgIdPFMet);
   fChain->SetBranchAddress("nDauPFMet", nDauPFMet, &b_nDauPFMet);
   fChain->SetBranchAddress("d1IndexPFMet", d1IndexPFMet, &b_d1IndexPFMet);
   fChain->SetBranchAddress("d2IndexPFMet", d2IndexPFMet, &b_d2IndexPFMet);
   fChain->SetBranchAddress("d1pdgIdPFMet", d1pdgIdPFMet, &b_d1pdgIdPFMet);
   fChain->SetBranchAddress("d2pdgIdPFMet", d2pdgIdPFMet, &b_d2pdgIdPFMet);
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
   fChain->SetBranchAddress("nIterativeJet", &nIterativeJet, &b_nIterativeJet);
   fChain->SetBranchAddress("chargeIterativeJet", chargeIterativeJet, &b_chargeIterativeJet);
   fChain->SetBranchAddress("energyIterativeJet", energyIterativeJet, &b_energyIterativeJet);
   fChain->SetBranchAddress("etIterativeJet", etIterativeJet, &b_etIterativeJet);
   fChain->SetBranchAddress("momentumIterativeJet", momentumIterativeJet, &b_momentumIterativeJet);
   fChain->SetBranchAddress("thetaIterativeJet", thetaIterativeJet, &b_thetaIterativeJet);
   fChain->SetBranchAddress("etaIterativeJet", etaIterativeJet, &b_etaIterativeJet);
   fChain->SetBranchAddress("phiIterativeJet", phiIterativeJet, &b_phiIterativeJet);
   fChain->SetBranchAddress("pxIterativeJet", pxIterativeJet, &b_pxIterativeJet);
   fChain->SetBranchAddress("pyIterativeJet", pyIterativeJet, &b_pyIterativeJet);
   fChain->SetBranchAddress("pzIterativeJet", pzIterativeJet, &b_pzIterativeJet);
   fChain->SetBranchAddress("vertexXIterativeJet", vertexXIterativeJet, &b_vertexXIterativeJet);
   fChain->SetBranchAddress("vertexYIterativeJet", vertexYIterativeJet, &b_vertexYIterativeJet);
   fChain->SetBranchAddress("vertexZIterativeJet", vertexZIterativeJet, &b_vertexZIterativeJet);
   fChain->SetBranchAddress("massIterativeJet", massIterativeJet, &b_massIterativeJet);
   fChain->SetBranchAddress("mtIterativeJet", mtIterativeJet, &b_mtIterativeJet);
   fChain->SetBranchAddress("pdgIdIterativeJet", pdgIdIterativeJet, &b_pdgIdIterativeJet);
   fChain->SetBranchAddress("nDauIterativeJet", nDauIterativeJet, &b_nDauIterativeJet);
   fChain->SetBranchAddress("d1IndexIterativeJet", d1IndexIterativeJet, &b_d1IndexIterativeJet);
   fChain->SetBranchAddress("d2IndexIterativeJet", d2IndexIterativeJet, &b_d2IndexIterativeJet);
   fChain->SetBranchAddress("d1pdgIdIterativeJet", d1pdgIdIterativeJet, &b_d1pdgIdIterativeJet);
   fChain->SetBranchAddress("d2pdgIdIterativeJet", d2pdgIdIterativeJet, &b_d2pdgIdIterativeJet);
   fChain->SetBranchAddress("alphaIterativeJet", alphaIterativeJet, &b_alphaIterativeJet);
   fChain->SetBranchAddress("emFracIterativeJet", emFracIterativeJet, &b_emFracIterativeJet);
   fChain->SetBranchAddress("hadFracIterativeJet", hadFracIterativeJet, &b_hadFracIterativeJet);
   fChain->SetBranchAddress("nSisConeJet", &nSisConeJet, &b_nSisConeJet);
   fChain->SetBranchAddress("chargeSisConeJet", chargeSisConeJet, &b_chargeSisConeJet);
   fChain->SetBranchAddress("energySisConeJet", energySisConeJet, &b_energySisConeJet);
   fChain->SetBranchAddress("etSisConeJet", etSisConeJet, &b_etSisConeJet);
   fChain->SetBranchAddress("momentumSisConeJet", momentumSisConeJet, &b_momentumSisConeJet);
   fChain->SetBranchAddress("thetaSisConeJet", thetaSisConeJet, &b_thetaSisConeJet);
   fChain->SetBranchAddress("etaSisConeJet", etaSisConeJet, &b_etaSisConeJet);
   fChain->SetBranchAddress("phiSisConeJet", phiSisConeJet, &b_phiSisConeJet);
   fChain->SetBranchAddress("pxSisConeJet", pxSisConeJet, &b_pxSisConeJet);
   fChain->SetBranchAddress("pySisConeJet", pySisConeJet, &b_pySisConeJet);
   fChain->SetBranchAddress("pzSisConeJet", pzSisConeJet, &b_pzSisConeJet);
   fChain->SetBranchAddress("vertexXSisConeJet", vertexXSisConeJet, &b_vertexXSisConeJet);
   fChain->SetBranchAddress("vertexYSisConeJet", vertexYSisConeJet, &b_vertexYSisConeJet);
   fChain->SetBranchAddress("vertexZSisConeJet", vertexZSisConeJet, &b_vertexZSisConeJet);
   fChain->SetBranchAddress("massSisConeJet", massSisConeJet, &b_massSisConeJet);
   fChain->SetBranchAddress("mtSisConeJet", mtSisConeJet, &b_mtSisConeJet);
   fChain->SetBranchAddress("pdgIdSisConeJet", pdgIdSisConeJet, &b_pdgIdSisConeJet);
   fChain->SetBranchAddress("nDauSisConeJet", nDauSisConeJet, &b_nDauSisConeJet);
   fChain->SetBranchAddress("d1IndexSisConeJet", d1IndexSisConeJet, &b_d1IndexSisConeJet);
   fChain->SetBranchAddress("d2IndexSisConeJet", d2IndexSisConeJet, &b_d2IndexSisConeJet);
   fChain->SetBranchAddress("d1pdgIdSisConeJet", d1pdgIdSisConeJet, &b_d1pdgIdSisConeJet);
   fChain->SetBranchAddress("d2pdgIdSisConeJet", d2pdgIdSisConeJet, &b_d2pdgIdSisConeJet);
   fChain->SetBranchAddress("alphaSisConeJet", alphaSisConeJet, &b_alphaSisConeJet);
   fChain->SetBranchAddress("emFracSisConeJet", emFracSisConeJet, &b_emFracSisConeJet);
   fChain->SetBranchAddress("hadFracSisConeJet", hadFracSisConeJet, &b_hadFracSisConeJet);
   fChain->SetBranchAddress("nIterativePFJet", &nIterativePFJet, &b_nIterativePFJet);
   fChain->SetBranchAddress("chargeIterativePFJet", chargeIterativePFJet, &b_chargeIterativePFJet);
   fChain->SetBranchAddress("energyIterativePFJet", energyIterativePFJet, &b_energyIterativePFJet);
   fChain->SetBranchAddress("etIterativePFJet", etIterativePFJet, &b_etIterativePFJet);
   fChain->SetBranchAddress("momentumIterativePFJet", momentumIterativePFJet, &b_momentumIterativePFJet);
   fChain->SetBranchAddress("thetaIterativePFJet", thetaIterativePFJet, &b_thetaIterativePFJet);
   fChain->SetBranchAddress("etaIterativePFJet", etaIterativePFJet, &b_etaIterativePFJet);
   fChain->SetBranchAddress("phiIterativePFJet", phiIterativePFJet, &b_phiIterativePFJet);
   fChain->SetBranchAddress("pxIterativePFJet", pxIterativePFJet, &b_pxIterativePFJet);
   fChain->SetBranchAddress("pyIterativePFJet", pyIterativePFJet, &b_pyIterativePFJet);
   fChain->SetBranchAddress("pzIterativePFJet", pzIterativePFJet, &b_pzIterativePFJet);
   fChain->SetBranchAddress("vertexXIterativePFJet", vertexXIterativePFJet, &b_vertexXIterativePFJet);
   fChain->SetBranchAddress("vertexYIterativePFJet", vertexYIterativePFJet, &b_vertexYIterativePFJet);
   fChain->SetBranchAddress("vertexZIterativePFJet", vertexZIterativePFJet, &b_vertexZIterativePFJet);
   fChain->SetBranchAddress("massIterativePFJet", massIterativePFJet, &b_massIterativePFJet);
   fChain->SetBranchAddress("mtIterativePFJet", mtIterativePFJet, &b_mtIterativePFJet);
   fChain->SetBranchAddress("pdgIdIterativePFJet", pdgIdIterativePFJet, &b_pdgIdIterativePFJet);
   fChain->SetBranchAddress("nDauIterativePFJet", nDauIterativePFJet, &b_nDauIterativePFJet);
   fChain->SetBranchAddress("d1IndexIterativePFJet", d1IndexIterativePFJet, &b_d1IndexIterativePFJet);
   fChain->SetBranchAddress("d2IndexIterativePFJet", d2IndexIterativePFJet, &b_d2IndexIterativePFJet);
   fChain->SetBranchAddress("d1pdgIdIterativePFJet", d1pdgIdIterativePFJet, &b_d1pdgIdIterativePFJet);
   fChain->SetBranchAddress("d2pdgIdIterativePFJet", d2pdgIdIterativePFJet, &b_d2pdgIdIterativePFJet);
   fChain->SetBranchAddress("chargedHadronEnergyIterativePFJet", chargedHadronEnergyIterativePFJet, &b_chargedHadronEnergyIterativePFJet);
   fChain->SetBranchAddress("neutralHadronEnergyIterativePFJet", neutralHadronEnergyIterativePFJet, &b_neutralHadronEnergyIterativePFJet);
   fChain->SetBranchAddress("chargedEmEnergyIterativePFJet", chargedEmEnergyIterativePFJet, &b_chargedEmEnergyIterativePFJet);
   fChain->SetBranchAddress("neutralEmEnergyIterativePFJet", neutralEmEnergyIterativePFJet, &b_neutralEmEnergyIterativePFJet);
   fChain->SetBranchAddress("neutralMultiplicityIterativePFJet", neutralMultiplicityIterativePFJet, &b_neutralMultiplicityIterativePFJet);
   fChain->SetBranchAddress("chargedMultiplicityIterativePFJet", chargedMultiplicityIterativePFJet, &b_chargedMultiplicityIterativePFJet);
   fChain->SetBranchAddress("muonMultiplicityIterativePFJet", muonMultiplicityIterativePFJet, &b_muonMultiplicityIterativePFJet);
   fChain->SetBranchAddress("nSisConePFJet", &nSisConePFJet, &b_nSisConePFJet);
   fChain->SetBranchAddress("chargeSisConePFJet", chargeSisConePFJet, &b_chargeSisConePFJet);
   fChain->SetBranchAddress("energySisConePFJet", energySisConePFJet, &b_energySisConePFJet);
   fChain->SetBranchAddress("etSisConePFJet", etSisConePFJet, &b_etSisConePFJet);
   fChain->SetBranchAddress("momentumSisConePFJet", momentumSisConePFJet, &b_momentumSisConePFJet);
   fChain->SetBranchAddress("thetaSisConePFJet", thetaSisConePFJet, &b_thetaSisConePFJet);
   fChain->SetBranchAddress("etaSisConePFJet", etaSisConePFJet, &b_etaSisConePFJet);
   fChain->SetBranchAddress("phiSisConePFJet", phiSisConePFJet, &b_phiSisConePFJet);
   fChain->SetBranchAddress("pxSisConePFJet", pxSisConePFJet, &b_pxSisConePFJet);
   fChain->SetBranchAddress("pySisConePFJet", pySisConePFJet, &b_pySisConePFJet);
   fChain->SetBranchAddress("pzSisConePFJet", pzSisConePFJet, &b_pzSisConePFJet);
   fChain->SetBranchAddress("vertexXSisConePFJet", vertexXSisConePFJet, &b_vertexXSisConePFJet);
   fChain->SetBranchAddress("vertexYSisConePFJet", vertexYSisConePFJet, &b_vertexYSisConePFJet);
   fChain->SetBranchAddress("vertexZSisConePFJet", vertexZSisConePFJet, &b_vertexZSisConePFJet);
   fChain->SetBranchAddress("massSisConePFJet", massSisConePFJet, &b_massSisConePFJet);
   fChain->SetBranchAddress("mtSisConePFJet", mtSisConePFJet, &b_mtSisConePFJet);
   fChain->SetBranchAddress("pdgIdSisConePFJet", pdgIdSisConePFJet, &b_pdgIdSisConePFJet);
   fChain->SetBranchAddress("nDauSisConePFJet", nDauSisConePFJet, &b_nDauSisConePFJet);
   fChain->SetBranchAddress("d1IndexSisConePFJet", d1IndexSisConePFJet, &b_d1IndexSisConePFJet);
   fChain->SetBranchAddress("d2IndexSisConePFJet", d2IndexSisConePFJet, &b_d2IndexSisConePFJet);
   fChain->SetBranchAddress("d1pdgIdSisConePFJet", d1pdgIdSisConePFJet, &b_d1pdgIdSisConePFJet);
   fChain->SetBranchAddress("d2pdgIdSisConePFJet", d2pdgIdSisConePFJet, &b_d2pdgIdSisConePFJet);
   fChain->SetBranchAddress("chargedHadronEnergySisConePFJet", chargedHadronEnergySisConePFJet, &b_chargedHadronEnergySisConePFJet);
   fChain->SetBranchAddress("neutralHadronEnergySisConePFJet", neutralHadronEnergySisConePFJet, &b_neutralHadronEnergySisConePFJet);
   fChain->SetBranchAddress("chargedEmEnergySisConePFJet", chargedEmEnergySisConePFJet, &b_chargedEmEnergySisConePFJet);
   fChain->SetBranchAddress("neutralEmEnergySisConePFJet", neutralEmEnergySisConePFJet, &b_neutralEmEnergySisConePFJet);
   fChain->SetBranchAddress("neutralMultiplicitySisConePFJet", neutralMultiplicitySisConePFJet, &b_neutralMultiplicitySisConePFJet);
   fChain->SetBranchAddress("chargedMultiplicitySisConePFJet", chargedMultiplicitySisConePFJet, &b_chargedMultiplicitySisConePFJet);
   fChain->SetBranchAddress("muonMultiplicitySisConePFJet", muonMultiplicitySisConePFJet, &b_muonMultiplicitySisConePFJet);
   fChain->SetBranchAddress("nIterativeGenJet", &nIterativeGenJet, &b_nIterativeGenJet);
   fChain->SetBranchAddress("chargeIterativeGenJet", chargeIterativeGenJet, &b_chargeIterativeGenJet);
   fChain->SetBranchAddress("energyIterativeGenJet", energyIterativeGenJet, &b_energyIterativeGenJet);
   fChain->SetBranchAddress("etIterativeGenJet", etIterativeGenJet, &b_etIterativeGenJet);
   fChain->SetBranchAddress("momentumIterativeGenJet", momentumIterativeGenJet, &b_momentumIterativeGenJet);
   fChain->SetBranchAddress("thetaIterativeGenJet", thetaIterativeGenJet, &b_thetaIterativeGenJet);
   fChain->SetBranchAddress("etaIterativeGenJet", etaIterativeGenJet, &b_etaIterativeGenJet);
   fChain->SetBranchAddress("phiIterativeGenJet", phiIterativeGenJet, &b_phiIterativeGenJet);
   fChain->SetBranchAddress("pxIterativeGenJet", pxIterativeGenJet, &b_pxIterativeGenJet);
   fChain->SetBranchAddress("pyIterativeGenJet", pyIterativeGenJet, &b_pyIterativeGenJet);
   fChain->SetBranchAddress("pzIterativeGenJet", pzIterativeGenJet, &b_pzIterativeGenJet);
   fChain->SetBranchAddress("vertexXIterativeGenJet", vertexXIterativeGenJet, &b_vertexXIterativeGenJet);
   fChain->SetBranchAddress("vertexYIterativeGenJet", vertexYIterativeGenJet, &b_vertexYIterativeGenJet);
   fChain->SetBranchAddress("vertexZIterativeGenJet", vertexZIterativeGenJet, &b_vertexZIterativeGenJet);
   fChain->SetBranchAddress("massIterativeGenJet", massIterativeGenJet, &b_massIterativeGenJet);
   fChain->SetBranchAddress("mtIterativeGenJet", mtIterativeGenJet, &b_mtIterativeGenJet);
   fChain->SetBranchAddress("pdgIdIterativeGenJet", pdgIdIterativeGenJet, &b_pdgIdIterativeGenJet);
   fChain->SetBranchAddress("nDauIterativeGenJet", nDauIterativeGenJet, &b_nDauIterativeGenJet);
   fChain->SetBranchAddress("d1IndexIterativeGenJet", d1IndexIterativeGenJet, &b_d1IndexIterativeGenJet);
   fChain->SetBranchAddress("d2IndexIterativeGenJet", d2IndexIterativeGenJet, &b_d2IndexIterativeGenJet);
   fChain->SetBranchAddress("d1pdgIdIterativeGenJet", d1pdgIdIterativeGenJet, &b_d1pdgIdIterativeGenJet);
   fChain->SetBranchAddress("d2pdgIdIterativeGenJet", d2pdgIdIterativeGenJet, &b_d2pdgIdIterativeGenJet);
   fChain->SetBranchAddress("nSisConeGenJet", &nSisConeGenJet, &b_nSisConeGenJet);
   fChain->SetBranchAddress("chargeSisConeGenJet", chargeSisConeGenJet, &b_chargeSisConeGenJet);
   fChain->SetBranchAddress("energySisConeGenJet", energySisConeGenJet, &b_energySisConeGenJet);
   fChain->SetBranchAddress("etSisConeGenJet", etSisConeGenJet, &b_etSisConeGenJet);
   fChain->SetBranchAddress("momentumSisConeGenJet", momentumSisConeGenJet, &b_momentumSisConeGenJet);
   fChain->SetBranchAddress("thetaSisConeGenJet", thetaSisConeGenJet, &b_thetaSisConeGenJet);
   fChain->SetBranchAddress("etaSisConeGenJet", etaSisConeGenJet, &b_etaSisConeGenJet);
   fChain->SetBranchAddress("phiSisConeGenJet", phiSisConeGenJet, &b_phiSisConeGenJet);
   fChain->SetBranchAddress("pxSisConeGenJet", pxSisConeGenJet, &b_pxSisConeGenJet);
   fChain->SetBranchAddress("pySisConeGenJet", pySisConeGenJet, &b_pySisConeGenJet);
   fChain->SetBranchAddress("pzSisConeGenJet", pzSisConeGenJet, &b_pzSisConeGenJet);
   fChain->SetBranchAddress("vertexXSisConeGenJet", vertexXSisConeGenJet, &b_vertexXSisConeGenJet);
   fChain->SetBranchAddress("vertexYSisConeGenJet", vertexYSisConeGenJet, &b_vertexYSisConeGenJet);
   fChain->SetBranchAddress("vertexZSisConeGenJet", vertexZSisConeGenJet, &b_vertexZSisConeGenJet);
   fChain->SetBranchAddress("massSisConeGenJet", massSisConeGenJet, &b_massSisConeGenJet);
   fChain->SetBranchAddress("mtSisConeGenJet", mtSisConeGenJet, &b_mtSisConeGenJet);
   fChain->SetBranchAddress("pdgIdSisConeGenJet", pdgIdSisConeGenJet, &b_pdgIdSisConeGenJet);
   fChain->SetBranchAddress("nDauSisConeGenJet", nDauSisConeGenJet, &b_nDauSisConeGenJet);
   fChain->SetBranchAddress("d1IndexSisConeGenJet", d1IndexSisConeGenJet, &b_d1IndexSisConeGenJet);
   fChain->SetBranchAddress("d2IndexSisConeGenJet", d2IndexSisConeGenJet, &b_d2IndexSisConeGenJet);
   fChain->SetBranchAddress("d1pdgIdSisConeGenJet", d1pdgIdSisConeGenJet, &b_d1pdgIdSisConeGenJet);
   fChain->SetBranchAddress("d2pdgIdSisConeGenJet", d2pdgIdSisConeGenJet, &b_d2pdgIdSisConeGenJet);
   Notify();
}

Bool_t EgammaBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
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
