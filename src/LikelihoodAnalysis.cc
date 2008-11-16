#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"


LikelihoodAnalysis::LikelihoodAnalysis(TTree *tree)
  : EgammaBase(tree) {
}



LikelihoodAnalysis::~LikelihoodAnalysis() {
}



float LikelihoodAnalysis::findEquivalentLHCut(float wantEfficiency) {

  int nbins = 500;
  TH1F *LHBinnedHisto = new TH1F("LHBinnedHisto", "LHBinnedHisto", nbins, 0.0, 1.0);


  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    for(int iele=0; iele<nEle; iele++) {
      if (eleClassEle[iele]==40) LHBinnedHisto->Fill(1.0);
      else {
	LHBinnedHisto->Fill(eleLikelihoodEle[iele]);
      }
    }

  }


  // scan the likelihood to find the cut
  float nEntries = LHBinnedHisto->GetEntries();
  float efficiency = 0.0;
  int bin = nbins+1;

  while (efficiency < wantEfficiency) {

    float integral = LHBinnedHisto->Integral(bin,nbins+1);
    efficiency = integral / nEntries;
    std::cout << "integral = " << integral 
	      << "\tefficiency = " << efficiency
	      << "bin = " << bin << std::endl;
    
    bin--;

  }

  float equivalentLHCut = LHBinnedHisto->GetBinLowEdge(bin);

  std::cout << "Equivalent cut on LH is LH > " << equivalentLHCut << std::endl
	    << "efficiency is: " << efficiency << std::endl
	    << "while wanted is: " << wantEfficiency << std::endl;

}




void LikelihoodAnalysis::reproduceEgammaCutID() {

  CutBasedEleIDSelector EgammaLooseCutBasedID;
  EgammaLooseCutBasedID.Configure("config/looseEleId"); 

  int nEvtGoodElectron = 0;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    int iSelected = -1;
    for(int iele=0; iele<nEle; iele++) {

      TVector3 pTrkAtOuter(pxAtOuterEle[iele],pyAtOuterEle[iele],pzAtOuterEle[iele]);

      EgammaLooseCutBasedID.SetHOverE( eleHoEEle[iele] );
      EgammaLooseCutBasedID.SetS9S25( s9s25Ele[iele] );
      EgammaLooseCutBasedID.SetDEta( eleDeltaEtaAtVtxEle[iele] );
      EgammaLooseCutBasedID.SetDPhiIn( eleDeltaPhiAtVtxEle[iele] );
      EgammaLooseCutBasedID.SetDPhiOut( eleDeltaPhiAtCaloEle[iele] );
      EgammaLooseCutBasedID.SetInvEminusInvP( 1./eleCaloCorrEEle[iele]-1./eleTrackerPEle[iele] );
      EgammaLooseCutBasedID.SetBremFraction( fabs(eleTrackerPEle[iele]-pTrkAtOuter.Mag())/eleTrackerPEle[iele] );
      EgammaLooseCutBasedID.SetSigmaEtaEta( sqrt(covEtaEtaEle[iele]) );
      EgammaLooseCutBasedID.SetSigmaPhiPhi( sqrt(covPhiPhiEle[iele]) );
      EgammaLooseCutBasedID.SetEOverPout( eleCorrEoPoutEle[iele] );
      EgammaLooseCutBasedID.SetEOverPin( eleCorrEoPEle[iele] );
      EgammaLooseCutBasedID.SetElectronClass ( eleClassEle[iele] );
      EgammaLooseCutBasedID.SetEgammaCutBasedID ( eleIdCutBasedEle[iele] );
      EgammaLooseCutBasedID.SetLikelihood( eleLikelihoodEle[iele] );      

      bool isEleIDCutBased = EgammaLooseCutBasedID.output();

      if( isEleIDCutBased ) iSelected=iele;
    }

    nEvtGoodElectron++;

  }

  EgammaLooseCutBasedID.diplayEfficiencies();


}



void LikelihoodAnalysis::estimateIDEfficiency() {

  // hardcoded cuts
  float maxRelativeTrackerPtSum = 0.08;
  float minLikelihood = 0.01;

  int nbinsEta = 40;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *GenEta   = new TH1F( "GenEta",  "generated #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *IsoEta   = new TH1F( "IsoEta",  "isolated #eta",      nbinsEta, minEta, maxEta );
  TH1F *CutIdEta = new TH1F( "CutIdEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdEta  = new TH1F( "LHIdEta",  "LH ID #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 40;
  float minPt = 5.0;
  float maxPt = 100.;
  
  TH1F *GenPt   = new TH1F( "GenPt",  "generated p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *IsoPt   = new TH1F( "IsoPt",  "isolated p_{T} (GeV)",      nbinsPt, minPt, maxPt );
  TH1F *CutIdPt = new TH1F( "CutIdPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *LHIdPt  = new TH1F( "LHIdPt",  "LH ID p_{T} (GeV)",        nbinsPt, minPt, maxPt );

  // e+e- in the Z->e+e-
  int indexeplus = 7;
  int indexeminus = 8;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // to have only BARREL
    //     if (fabs(etaMc[indexeplus]) < 1.476) GenEta->Fill(etaMc[indexeplus]);
    //     if (fabs(etaMc[indexeminus]) < 1.476) GenEta->Fill(etaMc[indexeminus]);
    //     if (fabs(etaMc[indexeplus]) < 1.476) GenPt->Fill(pMc[indexeplus] * fabs(sin(thetaMc[indexeplus])));
    //     if (fabs(etaMc[indexeminus]) < 1.476) GenPt->Fill(pMc[indexeminus] * fabs(sin(thetaMc[indexeminus])));

    GenEta->Fill(etaMc[indexeplus]);
    GenEta->Fill(etaMc[indexeminus]);
    GenPt->Fill(pMc[indexeplus] * fabs(sin(thetaMc[indexeplus])));
    GenPt->Fill(pMc[indexeminus] * fabs(sin(thetaMc[indexeminus])));

    // to avoid double counting
    bool RecoFilled[2], IsoFilled[2], CutIdFilled[2], LHIdFilled[2];
    for (int mcele=0; mcele<2; mcele++){
      RecoFilled[mcele] = false;
      IsoFilled[mcele]   = false;
      CutIdFilled[mcele] = false;
      LHIdFilled[mcele] = false;
    }

    // loop over ALL reconstructed electrons
    for(int iele=0; iele<nEle; iele++) {

      // matching with mc
      int matchedMcEle = -1;
      float mcEta, mcPt, deltaR;

      TVector3 recoElectronPAtInner(pxAtInnerEle[iele],pyAtInnerEle[iele],pzAtInnerEle[iele]);

      TVector3 trueElePeplus, trueElePeminus;
      trueElePeplus.SetPtThetaPhi(pMc[indexeplus]*fabs(sin(thetaMc[indexeplus])),thetaMc[indexeplus],phiMc[indexeplus]);
      trueElePeminus.SetPtThetaPhi(pMc[indexeminus]*fabs(sin(thetaMc[indexeminus])),thetaMc[indexeminus],phiMc[indexeminus]);

      float deltaReplus = trueElePeplus.DeltaR(recoElectronPAtInner);
      float deltaReminus = trueElePeminus.DeltaR(recoElectronPAtInner);
      
      if (deltaReplus < deltaReminus) {
	deltaR = deltaReplus;
	mcEta = etaMc[indexeplus];
	mcPt = pMc[indexeplus]*fabs(sin(thetaMc[indexeplus])); 
	matchedMcEle = 0;
      }
      if (deltaReminus < deltaReplus) { 
	deltaR = deltaReminus; 
	mcEta = etaMc[indexeminus]; 
	mcPt = pMc[indexeminus]*fabs(sin(thetaMc[indexeminus])); 
	matchedMcEle = 1;
      }

      //      to have only BARREL
      // if ( fabs(etaEle[iele]) < 1.476 ) {

	if ( !RecoFilled[matchedMcEle] ) {
	    
	RecoEta->Fill(mcEta);
	RecoPt->Fill(mcPt);
	RecoFilled[matchedMcEle] = true;
	    
      }
	  
      if ( eleSumPt04Ele[iele] < maxRelativeTrackerPtSum ) {

	if ( !IsoFilled[matchedMcEle] ) {
	      
	  IsoEta->Fill(mcEta);
	  IsoPt->Fill(mcPt);
	  IsoFilled[matchedMcEle] = true;
	      
	}
	    
	if ( eleIdCutBasedEle[iele]  && !CutIdFilled[matchedMcEle] ) {
	      
	  CutIdEta->Fill(mcEta);
	  CutIdPt->Fill(mcPt);
	  CutIdFilled[matchedMcEle] = true;
	      
	}

	if ( eleLikelihoodEle[iele] > minLikelihood && !LHIdFilled[matchedMcEle] ) {
	      
	  LHIdEta->Fill(mcEta);
	  LHIdPt->Fill(mcPt);
	  LHIdFilled[matchedMcEle] = true;
	      
	}

      } // isolation
	  
      //      } // only BARREL
	
    } // loop ele

  } // loop events

  EfficiencyEvaluator ElectronEffEtaCutBased("ElectronEffEtaCutBased.root");
  ElectronEffEtaCutBased.AddNumerator(GenEta);
  ElectronEffEtaCutBased.AddNumerator(RecoEta);
  ElectronEffEtaCutBased.AddNumerator(IsoEta);
  ElectronEffEtaCutBased.AddNumerator(CutIdEta);
  ElectronEffEtaCutBased.SetDenominator(GenEta);
  ElectronEffEtaCutBased.ComputeEfficiencies();
  ElectronEffEtaCutBased.SetTitle("electron efficiency vs #eta");
  ElectronEffEtaCutBased.SetXaxisTitle("electron #eta");
  ElectronEffEtaCutBased.SetYaxisTitle("efficiency");
  ElectronEffEtaCutBased.SetYaxisMin(0.0);
  ElectronEffEtaCutBased.DrawAll();
  ElectronEffEtaCutBased.Write();

  EfficiencyEvaluator ElectronEffEtaLikelihood("ElectronEffEtaLikelihood.root");
  ElectronEffEtaLikelihood.AddNumerator(GenEta);
  ElectronEffEtaLikelihood.AddNumerator(RecoEta);
  ElectronEffEtaLikelihood.AddNumerator(IsoEta);
  ElectronEffEtaLikelihood.AddNumerator(LHIdEta);
  ElectronEffEtaLikelihood.SetDenominator(GenEta);
  ElectronEffEtaLikelihood.ComputeEfficiencies();
  ElectronEffEtaLikelihood.SetTitle("electron efficiency vs #eta");
  ElectronEffEtaLikelihood.SetXaxisTitle("electron #eta");
  ElectronEffEtaLikelihood.SetYaxisTitle("efficiency");
  ElectronEffEtaLikelihood.SetYaxisMin(0.0);
  ElectronEffEtaLikelihood.DrawAll();
  ElectronEffEtaLikelihood.Write();

  EfficiencyEvaluator ElectronEffPtCutBased("ElectronEffPtCutBased.root");
  ElectronEffPtCutBased.AddNumerator(GenPt);
  ElectronEffPtCutBased.AddNumerator(RecoPt);
  ElectronEffPtCutBased.AddNumerator(IsoPt);
  ElectronEffPtCutBased.AddNumerator(CutIdPt);
  ElectronEffPtCutBased.SetDenominator(GenPt);
  ElectronEffPtCutBased.ComputeEfficiencies();
  ElectronEffPtCutBased.SetTitle("electron efficiency vs p_{T}");
  ElectronEffPtCutBased.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPtCutBased.SetYaxisTitle("efficiency");
  ElectronEffPtCutBased.SetYaxisMin(0.0);
  ElectronEffPtCutBased.DrawAll();
  ElectronEffPtCutBased.Write();

  EfficiencyEvaluator ElectronEffPtLikelihood("ElectronEffPtLikelihood.root");
  ElectronEffPtLikelihood.AddNumerator(GenPt);
  ElectronEffPtLikelihood.AddNumerator(RecoPt);
  ElectronEffPtLikelihood.AddNumerator(IsoPt);
  ElectronEffPtLikelihood.AddNumerator(LHIdPt);
  ElectronEffPtLikelihood.SetDenominator(GenPt);
  ElectronEffPtLikelihood.ComputeEfficiencies();
  ElectronEffPtLikelihood.SetTitle("electron efficiency vs p_{T}");
  ElectronEffPtLikelihood.SetXaxisTitle("electron p_{T}");
  ElectronEffPtLikelihood.SetYaxisTitle("efficiency");
  ElectronEffPtLikelihood.SetYaxisMin(0.0);
  ElectronEffPtLikelihood.DrawAll();
  ElectronEffPtLikelihood.Write();


}



void LikelihoodAnalysis::estimateFakeRate(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.812;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *FakeableJetsEta   = new TH1F( "FakeableJetsEta",  "fakeable jets #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *CutIdEta = new TH1F( "CutIdEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH ID tight #eta",        nbinsEta, minEta, maxEta );
  TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH ID loose #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 20;
  float minPt = 5.0;
  float maxPt = 100.;
  
  TH1F *FakeableJetsPt   = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *CutIdPt = new TH1F( "CutIdPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH ID tight p_{T} (GeV)",        nbinsPt, minPt, maxPt );
  TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH ID loose p_{T} (GeV)",        nbinsPt, minPt, maxPt );

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    bool tauPresence=false;

    int mceleindex = -1;
    for(int imc=0; imc<50; imc++) {
      // not only ele from W->enu: there is enu emission (V_ud?) in madgraph
      if ( (fabs(idMc[imc])==11) ) {
        mceleindex=imc;
        break;
      }
      // since the stable particle list is truncated, if there is a tau 
      // not possible to say what happens...
      if ( (fabs(idMc[imc])==15) ) {
        tauPresence=true;
        break;
      }
    }

    if(tauPresence) continue;

    //debug
//     if(mceleindex==-1) {
//       for(int imc=0; imc<20; imc++) {
//         cout << "imc = " << imc << "\tidMc = " << idMc[imc] << "\tmothMc = " << mothMc[imc] << endl;
//       }
//     }
    //enddebug

    TVector3 mcEle(0,0,0);
    if(mceleindex>-1) mcEle = TVector3(pMc[mceleindex]*cos(phiMc[mceleindex])*sin(thetaMc[mceleindex]), 
                                       pMc[mceleindex]*sin(phiMc[mceleindex])*sin(thetaMc[mceleindex]),
                                       pMc[mceleindex]*cos(thetaMc[mceleindex]));
    
//     int nJetsReco=0;
//     int nEleReco=0;

    for ( int jet=0; jet<nJet; jet++ ) {

      TVector3 p3Jet(pxJet[jet],pyJet[jet],pzJet[jet]);

       if ( fabs(etaJet[jet]) < 2.5 && etJet[jet] > 10.0 ) {
        
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Jet.DeltaR(mcEle);

        // remove from denominator the electron reconstructed as jet
        if( deltaR>0.3 ) { 
          FakeableJetsEta->Fill( etaJet[jet] );
          FakeableJetsPt->Fill( etJet[jet] );
          //          nJetsReco++;
        }

      }

    }

    for ( int ele=0; ele<nEle; ele++ ) {
      
      if ( fabs(etaEle[ele]) < 2.5 && etEle[ele] > 10.0 ) {

        TVector3 p3Ele(pxEle[ele], pyEle[ele], pzEle[ele]);
      
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Ele.DeltaR(mcEle);
        
        if(deltaR<0.3) continue;

        float dREleJet_min = 1000;
        int closestJet=-1;

        for ( int jet=0; jet<nJet; jet++ ) {

          if ( fabs(etaJet[jet]) < 2.5 && etJet[jet] > 10.0 ) {          

            TVector3 p3Jet(pxJet[jet],pyJet[jet],pzJet[jet]);
            
            float dREleJet = p3Jet.DeltaR(p3Ele);
            if(dREleJet<dREleJet_min) {
              closestJet=jet;
              dREleJet_min=dREleJet;
            }

          }

        }

        //        nEleReco++;

        if(closestJet > -1) {

          float etFake=etJet[closestJet];
          float etaFake=etaJet[closestJet];

          RecoEta->Fill(etaFake);
          RecoPt->Fill(etFake);

          if ( eleIdCutBasedEle[ele] ) {

            CutIdEta->Fill(etaFake);
            CutIdPt->Fill(etFake);

          }

          if ( eleLikelihoodEle[ele] > minLikelihoodLoose ) {

            LHIdLooseEta->Fill(etaFake);
            LHIdLoosePt->Fill(etFake);

          }

          if ( eleLikelihoodEle[ele] > minLikelihoodTight ) {

            LHIdTightEta->Fill(etaFake);
            LHIdTightPt->Fill(etFake);

          }

        }

      } // electron acceptance & pt cut
      
    } // loop ele

    //    cout << "nEleReco = " << nEleReco << " nJetsReco = " << nJetsReco << endl;
    
  } // loop events
  

  char filename[200];
  sprintf(filename,"%s-EleMisidEta.root",outname);
  EfficiencyEvaluator ElectronFakeRateEta(filename);
  ElectronFakeRateEta.AddNumerator(FakeableJetsEta);
  ElectronFakeRateEta.AddNumerator(RecoEta);
  ElectronFakeRateEta.AddNumerator(CutIdEta);
  ElectronFakeRateEta.AddNumerator(LHIdLooseEta);
  ElectronFakeRateEta.AddNumerator(LHIdTightEta);
  ElectronFakeRateEta.SetDenominator(FakeableJetsEta);
  ElectronFakeRateEta.ComputeEfficiencies();
  ElectronFakeRateEta.SetTitle("jet fake probability vs #eta");
  ElectronFakeRateEta.SetXaxisTitle("#eta of closest jet");
  ElectronFakeRateEta.SetYaxisTitle("jet #rightarrow fake e probability");
  ElectronFakeRateEta.SetYaxisMin(0.0);
  ElectronFakeRateEta.Write();

  sprintf(filename,"%s-EleMisidPt.root",outname);
  EfficiencyEvaluator ElectronFakeRatePt(filename);
  ElectronFakeRatePt.AddNumerator(FakeableJetsPt);
  ElectronFakeRatePt.AddNumerator(RecoPt);
  ElectronFakeRatePt.AddNumerator(CutIdPt);
  ElectronFakeRatePt.AddNumerator(LHIdLoosePt);
  ElectronFakeRatePt.AddNumerator(LHIdTightPt);
  ElectronFakeRatePt.SetDenominator(FakeableJetsPt);
  ElectronFakeRatePt.ComputeEfficiencies();
  ElectronFakeRatePt.SetTitle("jet fake probability vs p_{T}");
  ElectronFakeRatePt.SetXaxisTitle("p_{T} of closest jet [GeV]");
  ElectronFakeRatePt.SetYaxisTitle("jet #rightarrow fake e probability");
  ElectronFakeRatePt.SetYaxisMin(0.0);
  ElectronFakeRatePt.Write();

}



bool LikelihoodAnalysis::getCustomEleID(int eleIndex, const char *fileCuts, const char *fileSwitches) {

}
