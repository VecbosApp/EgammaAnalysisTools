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
      if (classificationEle[iele]==40) LHBinnedHisto->Fill(1.0);
      else {
	LHBinnedHisto->Fill(eleIdLikelihoodEle[iele]);
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

      EgammaLooseCutBasedID.SetHOverE( hOverEEle[iele] );
      EgammaLooseCutBasedID.SetS9S25( s9s25Ele[iele] );
      EgammaLooseCutBasedID.SetDEta( deltaEtaAtVtxEle[iele] );
      EgammaLooseCutBasedID.SetDPhiIn( deltaPhiAtVtxEle[iele] );
      EgammaLooseCutBasedID.SetDPhiOut( deltaPhiAtCaloEle[iele] );
      EgammaLooseCutBasedID.SetInvEminusInvP( 1./ecalEle[iele]-1./momentumEle[iele] );
      EgammaLooseCutBasedID.SetBremFraction( fabs(momentumEle[iele]-pTrkAtOuter.Mag())/momentumEle[iele] );
      EgammaLooseCutBasedID.SetSigmaEtaEta( sqrt(covEtaEtaEle[iele]) );
      EgammaLooseCutBasedID.SetSigmaPhiPhi( sqrt(covPhiPhiEle[iele]) );
      EgammaLooseCutBasedID.SetEOverPout( eSeedOverPoutEle[iele] );
      EgammaLooseCutBasedID.SetEOverPin( eSuperClusterOverPEle[iele] );
      EgammaLooseCutBasedID.SetElectronClass ( classificationEle[iele] );
      EgammaLooseCutBasedID.SetEgammaCutBasedID ( eleIdCutsEle[iele] );
      EgammaLooseCutBasedID.SetLikelihood( eleIdLikelihoodEle[iele] );      

      bool isEleIDCutBased = EgammaLooseCutBasedID.output();

      if( isEleIDCutBased ) iSelected=iele;
    }

    nEvtGoodElectron++;

  }

  EgammaLooseCutBasedID.diplayEfficiencies();


}



void LikelihoodAnalysis::estimateIDEfficiency(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.812;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *GenEta   = new TH1F( "GenEta",  "generated #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *NarrowEta = new TH1F( "NarrowEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *BigBremEta = new TH1F( "BigBremEta", "bigbrem@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *IsoEta   = new TH1F( "IsoEta",  "isolated #eta",      nbinsEta, minEta, maxEta );
  TH1F *CutIdEta = new TH1F( "CutIdEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH Loose ID #eta",        nbinsEta, minEta, maxEta );
  TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH Tight ID #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 20;
  float minPt = 10.0;
  float maxPt = 100.;
  
  TH1F *GenPt   = new TH1F( "GenPt",  "generated p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *NarrowPt = new TH1F( "NarrowPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *BigBremPt = new TH1F( "BigBremPt", "bigbrem@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *IsoPt   = new TH1F( "IsoPt",  "isolated p_{T} (GeV)",      nbinsPt, minPt, maxPt );
  TH1F *CutIdPt = new TH1F( "CutIdPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH Loose ID p_{T} (GeV)",        nbinsPt, minPt, maxPt );
  TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH Tight ID p_{T} (GeV)",        nbinsPt, minPt, maxPt );

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    int mceleindex = -1;
    for(int imc=0; imc<50; imc++) {
      if ( fabs(idMc[imc])==11 && fabs(idMc[mothMc[imc]])==24 ) mceleindex = imc;
    }
    
    if(mceleindex==-1) continue;

    TVector3 mcEle(pMc[mceleindex]*cos(phiMc[mceleindex])*sin(thetaMc[mceleindex]), 
                   pMc[mceleindex]*sin(phiMc[mceleindex])*sin(thetaMc[mceleindex]),
                   pMc[mceleindex]*cos(thetaMc[mceleindex]));

    float mcEta=etaMc[mceleindex];
    float mcPt=pMc[mceleindex] * fabs(sin(thetaMc[mceleindex]));

    // exclude forward electrons
    if(mcPt < minPt || fabs(mcEta) > maxEta) continue;

    GenEta->Fill(mcEta);
    GenPt->Fill(mcPt);

    // loop over ALL reconstructed electrons and find the closest one to the generated one
    float deltaR_min=0.3;
    int matchedRecoEle=-1;
    for(int iele=0; iele<nEle; iele++) {

      TVector3 eleP3(pxEle[iele],pyEle[iele],pzEle[iele]);
      float deltaR = eleP3.DeltaR(mcEle);
      if(deltaR < deltaR_min) {
        matchedRecoEle=iele;
        deltaR_min=deltaR;
      }

    }

    if(matchedRecoEle > -1) {

      RecoEta->Fill(mcEta);
      RecoPt->Fill(mcPt);

      int fullclass = classificationEle[matchedRecoEle];

      if ( fullclass == 0 || fullclass == 100 ) {
        GoldenEta->Fill(mcEta);
        GoldenPt->Fill(mcPt);
      } else if ( fullclass == 10 || fullclass == 110 ) {
        BigBremEta->Fill(mcEta);
        BigBremPt->Fill(mcPt);
      } else if ( fullclass == 20 || fullclass == 120 ) {
        NarrowEta->Fill(mcEta);
        NarrowPt->Fill(mcPt);
      } else if ( (fullclass >= 30 && fullclass <= 40) ||
                  (fullclass >= 130 && fullclass <= 140) ) {
        ShoweringEta->Fill(mcEta);
        ShoweringPt->Fill(mcPt);
      }
      
      if ( eleIdCutsEle[matchedRecoEle] ) {
        
        CutIdEta->Fill(mcEta);
        CutIdPt->Fill(mcPt);
        
      }
      
      if ( eleIdLikelihoodEle[matchedRecoEle] > minLikelihoodLoose ) {
        
        LHIdLooseEta->Fill(mcEta);
        LHIdLoosePt->Fill(mcPt);
        
      }
      
      if ( eleIdLikelihoodEle[matchedRecoEle] > minLikelihoodTight ) {
        
        LHIdTightEta->Fill(mcEta);
        LHIdTightPt->Fill(mcPt);
        
      }
      
    }

  } // loop events

  char filename[200];
  sprintf(filename,"%s-EleEfficiencyEta.root",outname);
  EfficiencyEvaluator ElectronEffEta(filename);
  ElectronEffEta.AddNumerator(GenEta);
  ElectronEffEta.AddNumerator(RecoEta);
  ElectronEffEta.AddNumerator(GoldenEta);
  ElectronEffEta.AddNumerator(BigBremEta);
  ElectronEffEta.AddNumerator(NarrowEta);
  ElectronEffEta.AddNumerator(ShoweringEta);
  ElectronEffEta.AddNumerator(CutIdEta);
  ElectronEffEta.AddNumerator(LHIdLooseEta);
  ElectronEffEta.AddNumerator(LHIdTightEta);
  ElectronEffEta.SetDenominator(GenEta);
  ElectronEffEta.ComputeEfficiencies();
  ElectronEffEta.SetTitle("electron efficiency vs #eta");
  ElectronEffEta.SetXaxisTitle("electron #eta");
  ElectronEffEta.SetYaxisTitle("efficiency");
  ElectronEffEta.SetYaxisMin(0.0);
  ElectronEffEta.Write();

  sprintf(filename,"%s-EleEfficiencyPt.root",outname);
  EfficiencyEvaluator ElectronEffPt(filename);
  ElectronEffPt.AddNumerator(GenPt);
  ElectronEffPt.AddNumerator(RecoPt);
  ElectronEffPt.AddNumerator(GoldenPt);
  ElectronEffPt.AddNumerator(BigBremPt);
  ElectronEffPt.AddNumerator(NarrowPt);
  ElectronEffPt.AddNumerator(ShoweringPt);
  ElectronEffPt.AddNumerator(CutIdPt);
  ElectronEffPt.AddNumerator(LHIdLoosePt);
  ElectronEffPt.AddNumerator(LHIdTightPt);
  ElectronEffPt.SetDenominator(GenPt);
  ElectronEffPt.ComputeEfficiencies();
  ElectronEffPt.SetTitle("electron efficiency vs p_{T}");
  ElectronEffPt.SetXaxisTitle("electron p_{T} (GeV)");
  ElectronEffPt.SetYaxisTitle("efficiency");
  ElectronEffPt.SetYaxisMin(0.0);
  ElectronEffPt.Write();

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
  TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *NarrowEta = new TH1F( "NarrowEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *BigBremEta = new TH1F( "BigBremEta", "bigbrem@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *CutIdEta = new TH1F( "CutIdEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH ID tight #eta",        nbinsEta, minEta, maxEta );
  TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH ID loose #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 20;
  float minPt = 10.0;
  float maxPt = 100.;
  
  TH1F *FakeableJetsPt   = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *NarrowPt = new TH1F( "NarrowPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *BigBremPt = new TH1F( "BigBremPt", "bigbrem@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
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

    for ( int jet=0; jet<nSisConeJet; jet++ ) {

      TVector3 p3Jet(pxSisConeJet[jet],pySisConeJet[jet],pzSisConeJet[jet]);

       if ( fabs(etaSisConeJet[jet]) < 2.5 && etSisConeJet[jet] > 10.0 ) {
        
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Jet.DeltaR(mcEle);

        // remove from denominator the electron reconstructed as jet
        if( deltaR>0.3 ) { 
          FakeableJetsEta->Fill( etaSisConeJet[jet] );
          FakeableJetsPt->Fill( etSisConeJet[jet] );
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

        for ( int jet=0; jet<nSisConeJet; jet++ ) {

          if ( fabs(etaSisConeJet[jet]) < 2.5 && etSisConeJet[jet] > 10.0 ) {          

            TVector3 p3Jet(pxSisConeJet[jet],pySisConeJet[jet],pzSisConeJet[jet]);
            
            float dREleJet = p3Jet.DeltaR(p3Ele);
            if(dREleJet<dREleJet_min) {
              closestJet=jet;
              dREleJet_min=dREleJet;
            }

          }

        }

        //        nEleReco++;

        if(closestJet > -1) {

          float etFake=etSisConeJet[closestJet];
          float etaFake=etaSisConeJet[closestJet];

          RecoEta->Fill(etaFake);
          RecoPt->Fill(etFake);

          int fullclass = classificationEle[ele];
          
          if ( fullclass == 0 || fullclass == 100 ) {
            GoldenEta->Fill(etaFake);
            GoldenPt->Fill(etFake);
          } else if ( fullclass == 10 || fullclass == 110 ) {
            BigBremEta->Fill(etaFake);
            BigBremPt->Fill(etFake);
          } else if ( fullclass == 20 || fullclass == 120 ) {
            NarrowEta->Fill(etaFake);
            NarrowPt->Fill(etFake);
          } else if ( (fullclass >= 30 && fullclass <= 40) ||
                      (fullclass >= 130 && fullclass <= 140) ) {
            ShoweringEta->Fill(etaFake);
            ShoweringPt->Fill(etFake);
          }

          if ( eleIdCutsEle[ele] ) {

            CutIdEta->Fill(etaFake);
            CutIdPt->Fill(etFake);

          }

          if ( eleIdLikelihoodEle[ele] > minLikelihoodLoose ) {

            LHIdLooseEta->Fill(etaFake);
            LHIdLoosePt->Fill(etFake);

          }

          if ( eleIdLikelihoodEle[ele] > minLikelihoodTight ) {

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
  ElectronFakeRateEta.AddNumerator(GoldenEta);
  ElectronFakeRateEta.AddNumerator(BigBremEta);
  ElectronFakeRateEta.AddNumerator(NarrowEta);
  ElectronFakeRateEta.AddNumerator(ShoweringEta);
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
  ElectronFakeRatePt.AddNumerator(GoldenPt);
  ElectronFakeRatePt.AddNumerator(BigBremPt);
  ElectronFakeRatePt.AddNumerator(NarrowPt);
  ElectronFakeRatePt.AddNumerator(ShoweringPt);
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
