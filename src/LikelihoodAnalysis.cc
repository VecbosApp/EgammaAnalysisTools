#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"

using namespace bits;
using namespace std;

LikelihoodAnalysis::LikelihoodAnalysis(TTree *tree)
  : Egamma(tree) {

  _isData = false;

  // single electron efficiency
  cout << "=== CONFIGURING WP80 SYMMETRIC ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/WP80");
  //  EgammaCutBasedID.ConfigureEcalCleaner("config/");

  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *EBlt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EBgt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = true;
  defaultSwitches.m_useSigmaPhiPhi = true;
  LH = new ElectronLikelihood(&(*EBlt15dir), &(*EElt15dir), &(*EBgt15dir), &(*EEgt15dir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);
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
      TVector3 pEle(pxEle[iele],pyEle[iele],pzEle[iele]);
      if(pEle.Pt()>20.0) LHBinnedHisto->Fill(likelihoodRatio(iele,*LH));
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

  TFile *outfile = TFile::Open("likHisto.root","recreate");
  LHBinnedHisto->Write();
  outfile->Close();

}




void LikelihoodAnalysis::reproduceEgammaCutID() {

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

      bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
      isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
      isEleID(&EgammaCutBasedID,iele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

      if( isEleIDCutBased ) iSelected=iele;
    }

    nEvtGoodElectron++;

  }

  EgammaCutBasedID.displayEfficiencies();


}



void LikelihoodAnalysis::estimateIDEfficiency(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *GenEta   = new TH1F( "GenEta",  "generated #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *IsoEta   = new TH1F( "IsoEta",  "isolated #eta",      nbinsEta, minEta, maxEta );
  TH1F *CutIdWP80Eta = new TH1F( "CutIdWP80Eta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *CutIdTightCICEta = new TH1F( "CutIdTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *CutIdSuperTightCICEta = new TH1F( "CutIdSuperTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH Loose ID #eta",        nbinsEta, minEta, maxEta );
  TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH Tight ID #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 20;
  float minPt = 20.0;
  float maxPt = 100.;
  
  TH1F *GenPt   = new TH1F( "GenPt",  "generated p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *IsoPt   = new TH1F( "IsoPt",  "isolated p_{T} (GeV)",      nbinsPt, minPt, maxPt );
  TH1F *CutIdWP80Pt = new TH1F( "CutIdWP80Pt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
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

      // Golden means 0 brem clusters, showering otherwise
      if ( nbremsEle[matchedRecoEle] == 0 ) {
        GoldenEta->Fill(mcEta);
        GoldenPt->Fill(mcPt);
      } else {
        ShoweringEta->Fill(mcEta);
        ShoweringPt->Fill(mcPt);
      }

      bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
      isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
      isEleID(&EgammaCutBasedID,matchedRecoEle,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);

      if ( isEleIDCutBased ) {

        CutIdWP80Eta->Fill(mcEta);
        CutIdWP80Pt->Fill(mcPt);

      }
      
      Utils anaUtils;

      if ( anaUtils.electronIdVal(eleIdCutsEle[matchedRecoEle], eleIdTightCIC) ) {

        CutIdTightCICEta->Fill(mcEta);
        CutIdTightCICPt->Fill(mcPt);
        
      }

      if ( anaUtils.electronIdVal(eleIdCutsEle[matchedRecoEle], eleIdSuperTightCIC) ) {

        CutIdSuperTightCICEta->Fill(mcEta);
        CutIdSuperTightCICPt->Fill(mcPt);
        
      }
      
      float lik = likelihoodRatio(matchedRecoEle,*LH);
      bool eleInGap = anaUtils.fiducialFlagECAL(fiducialFlagsEle[matchedRecoEle], isEBEEGap);

      if ( lik > minLikelihoodLoose && !eleInGap) {
        
        LHIdLooseEta->Fill(mcEta);
        LHIdLoosePt->Fill(mcPt);
        
      }
      
      if ( lik > minLikelihoodTight && !eleInGap) {
        
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
  ElectronEffEta.AddNumerator(ShoweringEta);
  ElectronEffEta.AddNumerator(CutIdWP80Eta);
  ElectronEffEta.AddNumerator(CutIdTightCICEta);
  ElectronEffEta.AddNumerator(CutIdSuperTightCICEta);
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
  ElectronEffPt.AddNumerator(ShoweringPt);
  ElectronEffPt.AddNumerator(CutIdWP80Pt);
  ElectronEffPt.AddNumerator(CutIdTightCICPt);
  ElectronEffPt.AddNumerator(CutIdSuperTightCICPt);
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
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *FakeableJetsEta   = new TH1F( "FakeableJetsEta",  "fakeable jets #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *CutIdWP80Eta = new TH1F( "CutIdWP80Eta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *CutIdTightCICEta = new TH1F( "CutIdTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *CutIdSuperTightCICEta = new TH1F( "CutIdSuperTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH ID tight #eta",        nbinsEta, minEta, maxEta );
  TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH ID loose #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 20;
  float minPt = 20.0;
  float maxPt = 100.;
  
  TH1F *FakeableJetsPt   = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *CutIdWP80Pt = new TH1F( "CutIdWP80Pt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
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

    for ( int jet=0; jet<nAK5PFJet; jet++ ) {

      TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);

      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 20.0 ) {
        
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Jet.DeltaR(mcEle);

        // remove from denominator the electron reconstructed as jet
        if( deltaR>0.5 ) { 
          FakeableJetsEta->Fill( p3Jet.Eta() );
          FakeableJetsPt->Fill( p3Jet.Pt() );
          //          nJetsReco++;
        }

      }

    }

    for ( int ele=0; ele<nEle; ele++ ) {

      TVector3 p3Ele(pxEle[ele], pyEle[ele], pzEle[ele]);
      
      if ( fabs(etaEle[ele]) < 2.5 && p3Ele.Pt() > 20.0 ) {

      
        float deltaR = 1000;
        if(mceleindex>-1) deltaR = p3Ele.DeltaR(mcEle);
        
        if(deltaR<0.5) continue;

        float dREleJet_min = 1000;
        int closestJet=-1;

        for ( int jet=0; jet<nAK5PFJet; jet++ ) {

          TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);

          if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 20.0 ) {          
            
            float dREleJet = p3Jet.DeltaR(p3Ele);
            if(dREleJet<dREleJet_min) {
              closestJet=jet;
              dREleJet_min=dREleJet;
            }

          }

        }

        //        nEleReco++;

        if(closestJet > -1) {

          TVector3 p3ClosestJet(pxAK5PFJet[closestJet],pyAK5PFJet[closestJet],pzAK5PFJet[closestJet]);

          float etFake=p3ClosestJet.Pt();
          float etaFake=p3ClosestJet.Eta();

          RecoEta->Fill(etaFake);
          RecoPt->Fill(etFake);

          if ( nbremsEle[ele] == 0 ) {
            GoldenEta->Fill(etaFake);
            GoldenPt->Fill(etFake);
          } else {
            ShoweringEta->Fill(etaFake);
            ShoweringPt->Fill(etFake);
          }

          bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
          isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
          isEleID(&EgammaCutBasedID,ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
          
          if ( isEleIDCutBased ) {

            CutIdWP80Eta->Fill(etaFake);
            CutIdWP80Pt->Fill(etFake);

          }
          
          Utils anaUtils;
          
          if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdTightCIC) ) {
            
            CutIdTightCICEta->Fill(etaFake);
            CutIdTightCICPt->Fill(etFake);
            
          }

          if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdSuperTightCIC) ) {
            
            CutIdSuperTightCICEta->Fill(etaFake);
            CutIdSuperTightCICPt->Fill(etFake);
            
          }

          float lik = likelihoodRatio(ele,*LH);
          bool eleInGap = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEBEEGap);

          if ( lik > minLikelihoodLoose && !eleInGap) {

            LHIdLooseEta->Fill(etaFake);
            LHIdLoosePt->Fill(etFake);

          }

          if ( lik > minLikelihoodTight && !eleInGap) {

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
  ElectronFakeRateEta.AddNumerator(ShoweringEta);
  ElectronFakeRateEta.AddNumerator(RecoEta);
  ElectronFakeRateEta.AddNumerator(CutIdWP80Eta);
  ElectronFakeRateEta.AddNumerator(CutIdTightCICEta);
  ElectronFakeRateEta.AddNumerator(CutIdSuperTightCICEta);
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
  ElectronFakeRatePt.AddNumerator(ShoweringPt);
  ElectronFakeRatePt.AddNumerator(RecoPt);
  ElectronFakeRatePt.AddNumerator(CutIdWP80Pt);
  ElectronFakeRatePt.AddNumerator(CutIdTightCICPt);
  ElectronFakeRatePt.AddNumerator(CutIdSuperTightCICPt);
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

void LikelihoodAnalysis::estimateFakeRateQCD(const char *outname) {

  // hardcoded cuts
  float minLikelihoodTight = 0.75;
  float minLikelihoodLoose = 0.056;

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
    
  TH1F *FakeableJetsEta   = new TH1F( "FakeableJetsEta",  "fakeable jets #eta",     nbinsEta, minEta, maxEta );
  TH1F *RecoEta  = new TH1F( "RecoEta", "reconstructed #eta", nbinsEta, minEta, maxEta );
  TH1F *GoldenEta = new TH1F( "GoldenEta", "golden@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *ShoweringEta = new TH1F( "ShoweringEta", "narrow@reconstruction #eta", nbinsEta, minEta, maxEta );
  TH1F *CutIdWP80Eta = new TH1F( "CutIdWP80Eta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *CutIdTightCICEta = new TH1F( "CutIdTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *CutIdSuperTightCICEta = new TH1F( "CutIdSuperTightCICEta", "cut ID #eta",       nbinsEta, minEta, maxEta );
  TH1F *LHIdTightEta  = new TH1F( "LHIdTightEta",  "LH ID tight #eta",        nbinsEta, minEta, maxEta );
  TH1F *LHIdLooseEta  = new TH1F( "LHIdLooseEta",  "LH ID loose #eta",        nbinsEta, minEta, maxEta );

  int nbinsPt = 20;
  float minPt = 20.0;
  float maxPt = 100.;
  
  TH1F *FakeableJetsPt   = new TH1F( "FakeableJetsPt",  "fakeable jets p_{T} (GeV)",     nbinsPt, minPt, maxPt );
  TH1F *RecoPt  = new TH1F( "RecoPt", "reconstructed p_{T} (GeV)", nbinsPt, minPt, maxPt );
  TH1F *GoldenPt = new TH1F( "GoldenPt", "golden@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *ShoweringPt = new TH1F( "ShoweringPt", "narrow@reconstruction vs p_{T}", nbinsPt, minPt, maxPt );
  TH1F *CutIdWP80Pt = new TH1F( "CutIdWP80Pt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *CutIdTightCICPt = new TH1F( "CutIdTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *CutIdSuperTightCICPt = new TH1F( "CutIdSuperTightCICPt", "cut ID p_{T} (GeV)",       nbinsPt, minPt, maxPt );
  TH1F *LHIdTightPt  = new TH1F( "LHIdTightPt",  "LH ID tight p_{T} (GeV)",        nbinsPt, minPt, maxPt );
  TH1F *LHIdLoosePt  = new TH1F( "LHIdLoosePt",  "LH ID loose p_{T} (GeV)",        nbinsPt, minPt, maxPt );

  // json 
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask
    bool newTriggerMask = false;
    if(_isData) newTriggerMask = true;
    reloadTriggerMask(newTriggerMask);
    
    // Good Run selection 
    if (_isData && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_isData && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    // QCD trigger 
    requiredTriggers.push_back("HLT_Jet15U");
    requiredTriggers.push_back("HLT_Jet30U");
    requiredTriggers.push_back("HLT_Jet50U");
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( !passedHLT ) continue;   

    // look for the leading jet (not considered as fakeable object to remove trigger bias)
    float maxEt = -1;
    int leadingJet = -1;
    for ( int jet=0; jet<nAK5PFJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 20.0 && p3Jet.Pt() > maxEt) {
        maxEt = p3Jet.Pt();
        leadingJet = jet;
      }
    }
    
    // consider the other as fakes (denominator)
    for ( int jet=0; jet<nAK5PFJet && jet!=leadingJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 20.0 ) {
          FakeableJetsEta->Fill( p3Jet.Eta() );
          FakeableJetsPt->Fill( p3Jet.Pt() );
      }
    }

    // consider the electrons near a jet (numerator)
    for ( int ele=0; ele<nEle; ele++ ) {
      TVector3 p3Ele(pxEle[ele], pyEle[ele], pzEle[ele]);
      if ( fabs(etaEle[ele]) < 2.5 && p3Ele.Pt() > 20.0 ) {

        float dREleJet_min = 1000;
        int closestJet=-1;

        for ( int jet=0; jet<nAK5PFJet && jet!=leadingJet; jet++ ) {
          TVector3 p3Jet(pxAK5PFJet[jet],pyAK5PFJet[jet],pzAK5PFJet[jet]);
          if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 20.0 ) {          
            float dREleJet = p3Jet.DeltaR(p3Ele);
            if(dREleJet<dREleJet_min) {
              closestJet=jet;
              dREleJet_min=dREleJet;
            }
          }
        }

        //        nEleReco++;

        if(closestJet > -1) {

          TVector3 p3ClosestJet(pxAK5PFJet[closestJet],pyAK5PFJet[closestJet],pzAK5PFJet[closestJet]);

          float etFake=p3ClosestJet.Pt();
          float etaFake=p3ClosestJet.Eta();

          RecoEta->Fill(etaFake);
          RecoPt->Fill(etFake);

          if ( nbremsEle[ele] == 0 ) {
            GoldenEta->Fill(etaFake);
            GoldenPt->Fill(etFake);
          } else {
            ShoweringEta->Fill(etaFake);
            ShoweringPt->Fill(etFake);
          }

          bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
          isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
          isEleID(&EgammaCutBasedID,ele,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
          
          if ( isEleIDCutBased ) {

            CutIdWP80Eta->Fill(etaFake);
            CutIdWP80Pt->Fill(etFake);

          }
          
          Utils anaUtils;
          
          if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdTightCIC) ) {
            
            CutIdTightCICEta->Fill(etaFake);
            CutIdTightCICPt->Fill(etFake);
            
          }

          if ( anaUtils.electronIdVal(eleIdCutsEle[ele], eleIdSuperTightCIC) ) {
            
            CutIdSuperTightCICEta->Fill(etaFake);
            CutIdSuperTightCICPt->Fill(etFake);
            
          }

          float lik = likelihoodRatio(ele,*LH);
          bool eleInGap = anaUtils.fiducialFlagECAL(fiducialFlagsEle[ele], isEBEEGap);

          if ( lik > minLikelihoodLoose && !eleInGap) {

            LHIdLooseEta->Fill(etaFake);
            LHIdLoosePt->Fill(etFake);

          }

          if ( lik > minLikelihoodTight && !eleInGap) {

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
  ElectronFakeRateEta.AddNumerator(ShoweringEta);
  ElectronFakeRateEta.AddNumerator(RecoEta);
  ElectronFakeRateEta.AddNumerator(CutIdWP80Eta);
  ElectronFakeRateEta.AddNumerator(CutIdTightCICEta);
  ElectronFakeRateEta.AddNumerator(CutIdSuperTightCICEta);
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
  ElectronFakeRatePt.AddNumerator(ShoweringPt);
  ElectronFakeRatePt.AddNumerator(RecoPt);
  ElectronFakeRatePt.AddNumerator(CutIdWP80Pt);
  ElectronFakeRatePt.AddNumerator(CutIdTightCICPt);
  ElectronFakeRatePt.AddNumerator(CutIdSuperTightCICPt);
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


void LikelihoodAnalysis::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);

  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];

  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxSC[sc];
      e4SwissCross = e4SwissCrossSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagSC[sc];
      seedTime = timeSC[sc];
      seedChi2 = chi2SC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->applyElectronIDOnPFlowElectrons(true);
  selector->SetHOverE( HoE );
  selector->SetS9S25( s9s25 );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetDPhiOut( dphiout );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetSigmaPhiPhi( spp );
  selector->SetEOverPout( eopout );
  selector->SetEOverPin( eop );
  selector->SetElectronClass ( classificationEle[eleIndex] );
  selector->SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  selector->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  selector->SetEcalIsolation( dr03EcalRecHitSumEtEle[eleIndex] );
  selector->SetTrkIsolation( dr03TkSumPtEle[eleIndex] );
  selector->SetHcalIsolation( dr03HcalTowerSumEtEle[eleIndex] );
  selector->SetCombinedIsolation( (dr03TkSumPtEle[eleIndex] + 
                                   TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + 
                                   dr03HcalTowerSumEtEle[eleIndex]) / pEle.Pt() );
  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
//   selector->m_cleaner->SetE1(e1);
//   selector->m_cleaner->SetE4SwissCross(e4SwissCross);
//   selector->m_cleaner->SetFiducialFlag(fidFlagSC);
//   selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
//   selector->m_cleaner->SetSeedTime(seedTime);
//   selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();

}

/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float LikelihoodAnalysis::SigmaiEiE(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

/// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
float LikelihoodAnalysis::SigmaiPiP(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

float LikelihoodAnalysis::likelihoodRatio(int eleIndex, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  
  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  measurements.pt = pEle.Pt();
  measurements.subdet = (fabs(etaEle[eleIndex])<1.479) ? 0 : 1;
  measurements.deltaPhi = deltaPhiAtVtxEle[eleIndex];
  measurements.deltaEta = deltaEtaAtVtxEle[eleIndex];
  measurements.eSeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  measurements.eSuperClusterOverP = eSuperClusterOverPEle[eleIndex];
  measurements.hadronicOverEm = hOverEEle[eleIndex];
  measurements.sigmaIEtaIEta = SigmaiEiE(eleIndex);
  measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
  measurements.fBrem = fbremEle[eleIndex];
  measurements.nBremClusters = nbremsEle[eleIndex];
  return lh.result(measurements);
}
