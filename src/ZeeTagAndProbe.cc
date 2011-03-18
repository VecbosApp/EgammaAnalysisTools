#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "EgammaAnalysisTools/include/ZeeTagAndProbe.hh"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"
#include <CommonTools/include/TriggerMask.hh>

using namespace bits;
using namespace std;

ZeeTagAndProbe::ZeeTagAndProbe(TTree *tree)
  : Egamma(tree) {

  jsonFile = "";
  lastFile = "";
  
  std::string fileCuts("config/ZeeTagAndProbe/cuts.txt");
  std::string fileSwitches("config/ZeeTagAndProbe/switches.txt");
  
  m_selection = new Selection(fileCuts,fileSwitches);
  configSelection(m_selection);

  // data or MC
  isData_ = m_selection->getSwitch("isData");

  // reading GoodRUN LS 
  std::cout << "[GoodRunLS]::goodRunLS is " << m_selection->getSwitch("goodRunLS") << " isData is " <<  m_selection->getSwitch("isData") << std::endl;

  // to read good run list
  if (m_selection->getSwitch("goodRunLS") && m_selection->getSwitch("isData")) {
    std::string goodRunGiasoneFile = "config/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  // single electron efficiency
  EgammaCutBasedIDWPs.push_back("WP95");
  EgammaCutBasedIDWPs.push_back("WP90");
  EgammaCutBasedIDWPs.push_back("WP80");
  EgammaCutBasedIDWPs.push_back("WP70");

  EgammaCiCBasedIDWPs.push_back("CiCVeryLoose");
  EgammaCiCBasedIDWPs.push_back("CiCLoose");
  EgammaCiCBasedIDWPs.push_back("CiCMedium");
  EgammaCiCBasedIDWPs.push_back("CiCTight");
  EgammaCiCBasedIDWPs.push_back("CiCSuperTight");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight2");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight3");
  EgammaCiCBasedIDWPs.push_back("CiCHyperTight4");

  EgammaLHBasedIDWPs.push_back("LHVeryLoose");
  EgammaLHBasedIDWPs.push_back("LHLoose");
  EgammaLHBasedIDWPs.push_back("LHMedium");
  EgammaLHBasedIDWPs.push_back("LHTight");
  EgammaLHBasedIDWPs.push_back("LHHyperTight");

  // single electron efficiency
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;
      
      char configDir[50];
      sprintf(configDir,"config/%s",EgammaCutBasedIDWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaCutBasedIDWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
      //  aSelector.ConfigureEcalCleaner("config/");
      EgammaCutBasedID.push_back(aSelector);
    }  

  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;

      char configDir[50];
      sprintf(configDir,"config/%s",EgammaLHBasedIDWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaLHBasedIDWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
  //  aSelector.ConfigureEcalCleaner("config/");
      EgammaLHBasedID.push_back(aSelector);
    }  

  for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
    {
      CiCBasedEleSelector aSelector;

      char configDir[50];
      std::cout << "===== Configuring " <<  EgammaCiCBasedIDWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.Configure(EgammaCiCBasedIDWPs[i],1,1,2);
  //  aSelector.ConfigureEcalCleaner("config/");
      EgammaCiCBasedID.push_back(aSelector);
    }  


  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *EBlt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EBgt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;

  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useHoverE = false;

  LH = new ElectronLikelihood(&(*EBlt15dir), &(*EElt15dir), &(*EBgt15dir), &(*EEgt15dir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);


}

ZeeTagAndProbe::~ZeeTagAndProbe() { 
  
  delete m_selection;
}


// PDF for probe electrons within the acceptance and loose isolated in the tracker
// the tag is the one with the best match to the Z
// the tag must be within the acceptance, tracker isolated and loose identified
void ZeeTagAndProbe::Loop(const char *treefilesuffix) {
  
  if(fChain == 0) return;
  
  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addElectronIdBits();
  reducedTree.addIsolations();
  reducedTree.addMore();
  
  // counters
  int allevents   = 0;
  int trigger     = 0;
  int twoele      = 0;
  int eleTot      = 0;
  int eleEta      = 0;
  int elePt       = 0;
  int invmass     = 0;
  int tagId       = 0;
  int tagIsol     = 0;
  int probeIsol   = 0;
  int tagProbeTot = 0;

  // json 
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask
    reloadTriggerMask(true);
    
    // Good Run selection 
    if (isData_ && m_selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (isData_ && m_selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    allevents++;
    
    // trigger
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerSignal") && !passedHLT ) continue;   
    trigger++;

    // best tag-probe pair = mee closest to Z mass
    float minpull = 100000.;
    float mass    = 1000.;
    float okmass  = 1000.;

    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);

      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);
	
	eleTot++;

        if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[iele1]) || !m_selection->passCut("etaEleAcc",etaEle[iele2]) ) ) continue;
	eleEta++;

        if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc",electron1.Pt()) || !m_selection->passCut("ptEleAcc",electron2.Pt()) ) ) continue;
	elePt++;
	
        mass = (electron1+electron2).M();
        float pull=fabs(mass-91.1876);
	if(pull < minpull) {
	  okmass  = mass;
          minpull = pull;
          electrons[0] = iele1;
          electrons[1] = iele2;
        }
      }
    }

    if (okmass<999) twoele++;

    // start the tag & probe
    if( m_selection->passCut("meeWindow",okmass) ) {

      if ( okmass>120 || okmass<60 ) cout << "BACO!" << endl;
      invmass++;
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele==1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }

        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],    pyEle[tag],  pzEle[tag],  energyEle[tag]);

	// various about probe
        int charge = chargeEle[probe];
        float pt   = probeP4.Pt();
        float eta  = etaEle[probe];

        // apply the electron ID tight (WP80) on the tag electron
	bool tagIdentified, tagIsolated, tagConvRej;
        tagIdentified = tagIsolated = tagConvRej = false;
        isEleID(&EgammaCutBasedID[2],tag,&tagIdentified,&tagIsolated,&tagConvRej);
        if (tagIdentified && tagIsolated && tagConvRej) {

          // look at the probe
          int CutBasedId[4], CutBasedIdOnlyID[4], CutBasedIdOnlyIso[4], CutBasedIdOnlyConv[4];
          int CiCBasedId[9], CiCBasedIdOnlyID[9], CiCBasedIdOnlyIso[9], CiCBasedIdOnlyConv[9];
          int LHBasedId[5], LHBasedIdOnlyID[5], LHBasedIdOnlyIso[5], LHBasedIdOnlyConv[5];
          
          for (int icut=0;icut<EgammaCutBasedIDWPs.size();++icut)
          {
            bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
            isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
            isEleID(&EgammaCutBasedID[icut],probe,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
            CutBasedIdOnlyID[icut] = (isEleIDCutBased) ? 1 : 0;
            CutBasedIdOnlyIso[icut] = (isIsolCutBased) ? 1 : 0;
            CutBasedIdOnlyConv[icut] = (isConvRejCutBased) ? 1 : 0;
            CutBasedId[icut] = (isEleIDCutBased && isIsolCutBased && isConvRejCutBased) ? 1 : 0;
          }

          for (int icut=0;icut<EgammaLHBasedIDWPs.size();++icut)
            {
              bool isEleIDLHBased, isIsolLHBased, isConvRejLHBased;
              isEleIDLHBased = isIsolLHBased = isConvRejLHBased = false;
              isEleID(&EgammaLHBasedID[icut],probe,&isEleIDLHBased,&isIsolLHBased,&isConvRejLHBased);
              LHBasedIdOnlyID[icut] = (isEleIDLHBased) ? 1 : 0;
              LHBasedIdOnlyIso[icut] = (isIsolLHBased) ? 1 : 0;
              LHBasedIdOnlyConv[icut] = (isConvRejLHBased) ? 1 : 0;
              LHBasedId[icut] = (isEleIDLHBased && isIsolLHBased && isConvRejLHBased) ? 1 : 0;
            }
              
          for (int icut=0;icut<EgammaCiCBasedIDWPs.size();++icut)
            {
              bool isEleIDCiCBased, isIsolCiCBased, isConvRejCiCBased;
              isEleIDCiCBased = isIsolCiCBased = isConvRejCiCBased = false;
              isEleID(&EgammaCiCBasedID[icut],probe,&isEleIDCiCBased,&isIsolCiCBased,&isConvRejCiCBased);
              CiCBasedIdOnlyID[icut] = (isEleIDCiCBased) ? 1 : 0;
              CiCBasedIdOnlyIso[icut] = (isIsolCiCBased) ? 1 : 0;
              CiCBasedIdOnlyConv[icut] = (isConvRejCiCBased) ? 1 : 0;
              CiCBasedId[icut] = (isEleIDCiCBased && isIsolCiCBased && isConvRejCiCBased) ? 1 : 0;
            }

          // some eleID variables
          float HoE, s9s25, deta, dphi, fbrem, see, spp, eopout, eop, nbrems, recoFlag;
          bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[probe], isEcalDriven);
          HoE = hOverEEle[probe];
          deta = deltaEtaAtVtxEle[probe];
          dphi = deltaPhiAtVtxEle[probe];
          fbrem = fbremEle[probe];
          nbrems = nbremsEle[probe];
          eopout = eSeedOverPoutEle[probe];
          eop = eSuperClusterOverPEle[probe];
          if(ecaldriven) {
            int sc = superClusterIndexEle[probe];
            s9s25 = e3x3SC[sc]/e5x5SC[sc];
            see = sqrt(covIEtaIEtaSC[sc]);
            spp = sqrt(covIPhiIPhiSC[sc]);
            recoFlag = recoFlagSC[sc];
          } else {
            int sc = PFsuperClusterIndexEle[probe];
            if(sc>-1) {
              s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
              see = sqrt(covIEtaIEtaPFSC[sc]);
              spp = sqrt(covIPhiIPhiPFSC[sc]);
              recoFlag = recoFlagPFSC[sc];
            } else {
              s9s25 = 999.;
              see = 999.;
              spp = 999.;
            }
          }

          float lh=likelihoodRatio(probe,*LH);

          // fill the reduced tree
	  reducedTree.fillVariables(eopout,eop,HoE,deta,dphi,s9s25,see,spp,fbrem,nbrems,pt,eta,charge);
          reducedTree.fillAttributesSignal(okmass);
          reducedTree.fillIsolations(dr03TkSumPtEle[probe] - rhoFastjet*TMath::Pi()*0.3*0.3,
                                     dr03EcalRecHitSumEtEle[probe] - rhoFastjet*TMath::Pi()*0.3*0.3,
                                     dr03HcalTowerSumEtFullConeEle[probe] - rhoFastjet*TMath::Pi()*0.3*0.3);
          reducedTree.fillMore(nPV,rhoFastjet);
          reducedTree.fillCutBasedIDBits(CutBasedId,CutBasedIdOnlyID,CutBasedIdOnlyIso,CutBasedIdOnlyConv);
          reducedTree.fillLHBasedIDBits(LHBasedId,LHBasedIdOnlyID,LHBasedIdOnlyIso,LHBasedIdOnlyConv);
          reducedTree.fillCiCBasedIDBits(CiCBasedId,CiCBasedIdOnlyID,CiCBasedIdOnlyIso,CiCBasedIdOnlyConv);
          reducedTree.store();

        } // tag identified
	
      } // loop over the 2 Z electrons
      
    } // end tag and probe
    
  } // loop over events
  
  cout << "statistics from Tag and Probe: " << endl;
  cout << "allevents   = " << allevents << endl;
  cout << "trigger     = " << trigger << endl;
  cout << "twoele      = " << twoele  << endl;
  cout << "invmass     = " << invmass << endl;
  cout << "tagProbeTot = " << tagProbeTot << endl;
  cout << "tagId       = " << tagId       << endl;
  cout << "tagIsol     = " << tagIsol     << endl;
  cout << "probeIsol   = " << probeIsol   << endl;  
  cout << "statistics from Tag and Probe - electrons: " << endl;
  cout << "eleTot      = " << eleTot  << endl;
  cout << "eleEta      = " << eleEta  << endl;
  cout << "elePt       = " << elePt   << endl;

  reducedTree.save();
}


void ZeeTagAndProbe::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

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
      e1 = eMaxPFSC[sc];
      e4SwissCross = e4SwissCrossPFSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagPFSC[sc];
      seedTime = timePFSC[sc];
      seedChi2 = chi2PFSC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  float lh=likelihoodRatio(eleIndex,*LH);
  bool isEleEB= anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex], isEB);
  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->applyElectronIDOnPFlowElectrons(true);
  selector->SetHOverE( HoE );
  selector->SetS9S25( s9s25 );
  selector->SetNBrem( nbremsEle[eleIndex] );
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
  selector->SetLikelihood( lh );
  selector->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetTrkIsolation( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  float combinedIso = 0.0;
  if (isEleEB) combinedIso = dr03TkSumPtEle[eleIndex] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else combinedIso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  selector->SetCombinedIsolation( (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / pEle.Pt() ); 

  selector->SetMissingHits( trackLostHitsGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
  //selector->m_cleaner->SetE1(e1);
  //selector->m_cleaner->SetE4SwissCross(e4SwissCross);
  //selector->m_cleaner->SetFiducialFlag(fidFlagSC);
  //selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
  //selector->m_cleaner->SetSeedTime(seedTime);
  //selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();

}

void ZeeTagAndProbe::isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);

  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eseedopin, eop, eseed,  sce, scet,sceta, rawe;
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
    eseed = seedEnergySC[sc];
    sce = energySC[sc];
    rawe = rawEnergySC[sc];
    scet = energySC[sc]/TMath::CosH(etaSC[sc]);
    sceta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxPFSC[sc];
      e4SwissCross = e4SwissCrossPFSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagPFSC[sc];
      seedTime = timePFSC[sc];
      seedChi2 = chi2PFSC[sc];
      eseed = seedEnergyPFSC[sc];
      sce = energyPFSC[sc];
      rawe = rawEnergyPFSC[sc];
      scet = energyPFSC[sc]/TMath::CosH(etaPFSC[sc]);
      sceta = etaPFSC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  
  eseedopin = eSeedOverPoutEle[eleIndex]*(1-fbremEle[eleIndex]);
  //  std::cout << eop << "," <<  eSeedOverPoutEle[eleIndex] << ","<< eseedopin*sce/rawe << "," << fbremEle[eleIndex] << "," << nbremsEle[eleIndex] << std::endl;;
  
  //  float lh=likelihoodRatio(eleIndex,*LH);
  selector->reset();
  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->SetSCEt(scet);
  selector->SetSCEta(sceta);
  selector->SetHOverE( HoE );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetEOverPin( eop );
  selector->SetESeedOverPin( eseedopin );
  selector->SetEcalIsolation( dr04EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4 );
  selector->SetTrkIsolation( dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3 );
  selector->SetHcalIsolation( dr04HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4);
  selector->SetMissingHits( trackLostHitsGsfTrack[gsf] );
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
  *eleIdOutput = selector->outputEleId();
  *isolOutput = selector->outputIso();
  *convRejOutput = selector->outputConv();

}


void ZeeTagAndProbe::configSelection(Selection* selection) {

  m_selection->addSwitch("isData");
  m_selection->addSwitch("goodRunLS");
  m_selection->addCut("meeWindow");
  m_selection->addCut("etaEleAcc");
  m_selection->addCut("ptEleAcc");
  m_selection->summary();

}

