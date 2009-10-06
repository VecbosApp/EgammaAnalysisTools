#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "EgammaAnalysisTools/include/LHPdfsProducer.hh"

using namespace bits;

LHPdfsProducer::LHPdfsProducer(TTree *tree)
  : EgammaBase(tree) {
  
  std::string fileCuts("config/LHPdfsProducer/cuts.txt");
  std::string fileSwitches("config/LHPdfsProducer/switches.txt");
  
  m_selection = new Selection(fileCuts,fileSwitches);
  m_selection->addSwitch("requireTriggerSignal");
  m_selection->addSwitch("requireTriggerQCDBack");
  m_selection->addSwitch("applyIsolationOnProbe");
  m_selection->addCut("meeWindow");
  m_selection->addCut("etaEleAcc");
  m_selection->addCut("ptEleAcc");
  m_selection->addCut("etaJetAcc");
  m_selection->addCut("ptJetAcc");
  m_selection->addCut("relSumPtTracks");
  m_selection->addCut("jetDeltaPhi");
  m_selection->addCut("jetInvMass");
  m_selection->summary();
  
  // single electron efficiency            
  EgammaCutBasedID.Configure("config/looseEleId/");
}

LHPdfsProducer::~LHPdfsProducer() { }


// PDF for probe electrons within the acceptance and loose isolated in the tracker
// the tag is the one with the best match to the Z
// the tag must be within the acceptance, tracker isolated and loose identified
void LHPdfsProducer::LoopZTagAndProbe(const char *treefilesuffix) {
  
  if(fChain == 0) return;
  
  bookHistos();

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addCategories();
  reducedTree.addMore();        // to find the best cut
  
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

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    Utils anaUtils;
    
    allevents++;

    // trigger: for october exercise skimmed samples, with HLT15 passed: can be switched off.
    if ( m_selection->getSwitch("requireTriggerSignal") ) { 
      bool passedHLT = anaUtils.getTriggersOR(m_requiredSignalTriggers, firedTrg);
      if ( !passedHLT ) continue;   
    }
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
        m_Zmass->Fill(mass);
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

      if ( okmass>110 || okmass<60 ) cout << "BACO!" << endl;
      invmass++;
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele=1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }
	
        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],    pyEle[tag],  pzEle[tag],  energyEle[tag]);

	// various about probe
        int charge = chargeEle[probe];
        float pt   = probeP4.Pt();
        float eta  = etaEle[probe];

        /// define the bins in which can be splitted the PDFs
        int iecal  = (fabs(etaEle[probe])<1.479) ? 0 : 1;
        int iptbin = (probeP4.Pt()<15.0) ? 0 : 1;
        int fullclassRaw = classificationEle[probe];
        int iclass     = -1;
        int ifullclass = -1;
        if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
        else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
        else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
        else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
        if (iclass>-1) tagProbeTot++;

        // apply the electron ID loose on the tag electron
	// float tagIdentified = anaUtils.electronIdVal(eleIdCutsEle[tag],eleIdRobustLoose);  
	float tagIdentified = isEleID(tag);
        if (tagIdentified) tagId++;

        // apply tracker isolation on the tag electron 
	float tagPt           = tagP4.Pt();
	float relativeIsolTag = dr04TkSumPtEle[tag]/tagPt;
        bool tagIsolated      = ( m_selection->passCut("relSumPtTracks",relativeIsolTag) );
	if (tagIsolated) tagIsol++;

        // apply tracker isolation on the probe electron
        bool probeIsolated = true;
	float probePt      = probeP4.Pt();
	float relativeIsolProbe = dr04TkSumPtEle[probe]/probePt;	
        if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
          probeIsolated = ( m_selection->passCut("relSumPtTracks",relativeIsolProbe) );
        }
	if(probeIsolated) probeIsol++;


	// some eleID variables
        float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[probe]));
        float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiEle[probe]));
        float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiEle[probe]));
        float s1s9          = s1s9Ele[probe];
        float s9s25         = s9s25Ele[probe];
        float lat           = latEle[probe];
        float etaLat        = etaLatEle[probe];
        float phiLat        = phiLatEle[probe];
        float a20           = a20Ele[probe];
        float a42           = a42Ele[probe];
	double dPhiCalo     = deltaPhiAtCaloEle[probe];
	double dPhiVtx      = deltaPhiAtVtxEle[probe];
	double dEtaVtx      = deltaEtaAtVtxEle[probe];
	double EoPout       = eSeedOverPoutEle[probe];
	double EoP          = eSuperClusterOverPEle[probe];
	double HoE          = hOverEEle[probe];

	
        /// fill the electron ID pdfs only if:
        /// the tag is loose isolated and identified (ALWAYS)
        /// the probe is loose isolated              (ONLY IF REQUIRED)
	if( tagIsolated && tagIdentified && probeIsolated && iclass>-1) {   
          
          dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
          dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
          dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
          EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
          HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
          sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
          sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
          sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
          s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
          s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );
          LATUnsplitEle           [iecal][iptbin] -> Fill ( lat );
          etaLATUnsplitEle        [iecal][iptbin] -> Fill ( etaLat );
          phiLATUnsplitEle        [iecal][iptbin] -> Fill ( phiLat );
          a20UnsplitEle           [iecal][iptbin] -> Fill ( a20 );
          a42UnsplitEle           [iecal][iptbin] -> Fill ( a42 );

          dPhiCaloClassEle      [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
          dPhiVtxClassEle       [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
          dEtaClassEle          [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
          EoPoutClassEle        [iecal][iptbin][iclass] -> Fill ( EoPout );
          HoEClassEle           [iecal][iptbin][iclass] -> Fill ( HoE );
          sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
          sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
          sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
          s1s9ClassEle          [iecal][iptbin][iclass] -> Fill ( s1s9 );
          s9s25ClassEle         [iecal][iptbin][iclass] -> Fill ( s9s25 );
          LATClassEle           [iecal][iptbin][iclass] -> Fill ( lat );
          etaLATClassEle        [iecal][iptbin][iclass] -> Fill ( etaLat );
          phiLATClassEle        [iecal][iptbin][iclass] -> Fill ( phiLat );
          a20ClassEle           [iecal][iptbin][iclass] -> Fill ( a20 );
          a42ClassEle           [iecal][iptbin][iclass] -> Fill ( a42 );

          dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
          dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
          dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
          EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
          HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
          sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
          sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
          sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
          s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
          s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
          LATFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( lat );
          etaLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( etaLat );
          phiLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( phiLat );
          a20FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a20 );
          a42FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a42 );
	  
          // fill the reduced tree
	  reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,s1s9,sigmaIEtaIEta);
          reducedTree.fillAttributesSignal(charge,eta,pt,okmass);
          reducedTree.fillCategories(iecal,iptbin,iclass);
          reducedTree.fillMore(relativeIsolTag, relativeIsolProbe);
          reducedTree.store();

        } // fill histograms
	
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


// PDF for the signal from a pure Z MC sample. Same requests as above
void LHPdfsProducer::LoopZ(const char *treefilesuffix) {

  if(fChain == 0) return;

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addCategories();
  
  // counters
  int allevents    = 0;
  int taus         = 0;
  int mc           = 0;
  int trigger      = 0;
  int twoele       = 0;
  int foundReco1   = 0;
  int foundReco2   = 0;
  int eleEta1      = 0;
  int elePt1       = 0;
  int tagProbeTot1 = 0;
  int probeIsol1   = 0;
  int eleEta2      = 0;
  int elePt2       = 0;
  int tagProbeTot2 = 0;
  int probeIsol2   = 0;

  // loop over events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    allevents++;

    // the stable particle list is truncated, if there is a tau not possible to say what happens...                                   
    bool tauPresence=false;
    for(int iMc=0; iMc<50; iMc++) {
      if ( (fabs(idMc[iMc])==15) ) { tauPresence=true; break; }
    }
    taus++;

    // to find the real electrons from Z
    int mcInd1 = -1;
    int mcInd2 = -1;
    for(int iMc=0; iMc<nMc; iMc++) {
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1==-1 )              { mcInd1=iMc; continue; }
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1!=-1 && mcInd2==-1) { mcInd2=iMc; break; }
    }
    TVector3 _mcEle1(0,0,0);
    TVector3 _mcEle2(0,0,0);
    if(mcInd1>-1) _mcEle1 = TVector3(pMc[mcInd1]*cos(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*sin(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*cos(thetaMc[mcInd1]));
    if(mcInd2>-1) _mcEle2 = TVector3(pMc[mcInd2]*cos(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*sin(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*cos(thetaMc[mcInd2]));
    if (mcInd1<0 || mcInd2<0) continue; 
    mc++;

    // trigger: for october exercise skimmed samples, with HLT15 passed: can be switched off.
    Utils anaUtils;
    if ( m_selection->getSwitch("requireTriggerSignal") ) { 
      bool passedHLT = anaUtils.getTriggersOR(m_requiredSignalTriggers, firedTrg);
      if ( !passedHLT ) continue;   
    }
    trigger++;

    // electrons matching MC truth
    float deltaRmin_mc1   = 999.;
    float deltaRmin_mc2   = 999.;
    int theClosestEle_mc1 = -1;
    int theClosestEle_mc2 = -1;
    for(int theEle=0; theEle<nEle; theEle++){
      
      TVector3 _p3Ele(pxEle[theEle],pyEle[theEle],pzEle[theEle]);
      float deltaR_mc1 = _p3Ele.DeltaR(_mcEle1);
      float deltaR_mc2 = _p3Ele.DeltaR(_mcEle2);
      
      if (deltaR_mc1<deltaRmin_mc1 && deltaR_mc1<0.3){
	deltaRmin_mc1     = deltaR_mc1;
	theClosestEle_mc1 = theEle;
      }

      if (deltaR_mc2<deltaRmin_mc2 && deltaR_mc2<0.3){
	deltaRmin_mc2     = deltaR_mc2;
	theClosestEle_mc2 = theEle;
      }
    }
    

    // fill PDFs using the electrons in the acceptance and matching the MC truth
    if (theClosestEle_mc1>-1) {
      
      foundReco1++;

      TLorentzVector electron(pxEle[theClosestEle_mc1],pyEle[theClosestEle_mc1],pzEle[theClosestEle_mc1],energyEle[theClosestEle_mc1]);
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[theClosestEle_mc1])) ) continue;
      eleEta1++;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc", electron.Pt())) ) continue;
      elePt1++;

      // various
      int charge = chargeEle[theClosestEle_mc1];
      float pt   = electron.Pt();
      float eta  = etaEle[theClosestEle_mc1];

      // define the bins in which can be splitted the PDFs
      int iecal  = (fabs( etaEle[theClosestEle_mc1])<1.479) ? 0 : 1;
      int iptbin = (pt<15.0) ? 0 : 1;
      int fullclassRaw = classificationEle[theClosestEle_mc1];
      int iclass = -1;
      int ifullclass = -1;
      if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
      else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
      else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
      else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
      if (iclass>-1) tagProbeTot1++;

      // apply loose tracker isolation on the first electron
      bool isolated1 = true;
      float relativeIsol1 = dr04TkSumPtEle[theClosestEle_mc1]/pt;	
      if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
	isolated1 = ( m_selection->passCut("relSumPtTracks",relativeIsol1) );
      }
      probeIsol1++;

      // some eleID variables
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[theClosestEle_mc1]));
      float s1s9          = s1s9Ele[theClosestEle_mc1];
      float s9s25         = s9s25Ele[theClosestEle_mc1];
      double dPhiVtx      = deltaPhiAtVtxEle[theClosestEle_mc1];
      double dEtaVtx      = deltaEtaAtVtxEle[theClosestEle_mc1];
      double EoPout       = eSeedOverPoutEle[theClosestEle_mc1];
      double EoP          = eSuperClusterOverPEle[theClosestEle_mc1];
      double HoE          = hOverEEle[theClosestEle_mc1];

      // fill the reduced tree     
      if( isolated1 && iclass>-1 ) {
	reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,s1s9,sigmaIEtaIEta);
	reducedTree.fillAttributesSignal(charge,eta,pt,9999.);
	reducedTree.fillCategories(iecal,iptbin,iclass);
	reducedTree.store();
      } // fill tree
    } // ok 1st electron


    // fill PDFs using the electrons in the acceptance and matching the MC truth
    if (theClosestEle_mc2>-1) {
      
      foundReco2++;

      TLorentzVector electron(pxEle[theClosestEle_mc2],pyEle[theClosestEle_mc2],pzEle[theClosestEle_mc2],energyEle[theClosestEle_mc2]);
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[theClosestEle_mc2])) ) continue;
      eleEta2++;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc", electron.Pt())) ) continue;
      elePt2++;

      // various
      int charge = chargeEle[theClosestEle_mc2];
      float pt   = electron.Pt();
      float eta  = etaEle[theClosestEle_mc2];

      // define the bins in which can be splitted the PDFs
      int iecal  = (fabs( etaEle[theClosestEle_mc2])<1.479) ? 0 : 1;
      int iptbin = (pt<15.0) ? 0 : 1;
      int fullclassRaw = classificationEle[theClosestEle_mc2];
      int iclass = -1;
      int ifullclass = -1;
      if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
      else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
      else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
      else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
      if (iclass>-1) tagProbeTot2++;

      // apply loose tracker isolation on the first electron
      bool isolated2 = true;
      float relativeIsol2 = dr04TkSumPtEle[theClosestEle_mc2]/pt;	
      if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
	isolated2 = ( m_selection->passCut("relSumPtTracks",relativeIsol2) );
      }
      probeIsol2++;

      // some eleID variables
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[theClosestEle_mc2]));
      float s1s9          = s1s9Ele[theClosestEle_mc2];
      float s9s25         = s9s25Ele[theClosestEle_mc2];
      double dPhiVtx      = deltaPhiAtVtxEle[theClosestEle_mc2];
      double dEtaVtx      = deltaEtaAtVtxEle[theClosestEle_mc2];
      double EoPout       = eSeedOverPoutEle[theClosestEle_mc2];
      double EoP          = eSuperClusterOverPEle[theClosestEle_mc2];
      double HoE          = hOverEEle[theClosestEle_mc2];

      // fill the reduced tree     
      if( isolated2 && iclass>-1 ) {
	reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,s1s9,sigmaIEtaIEta);
	reducedTree.fillAttributesSignal(charge,eta,pt,9999.);
	reducedTree.fillCategories(iecal,iptbin,iclass);
	reducedTree.store();
      } // fill tree
    } // ok 2nd electron
    


  } // loop over events
  

  cout << "statistics from MC: " << endl;
  cout << "allevents    = " << allevents    << endl;
  cout << "taus         = " << taus         << endl;  
  cout << "mc           = " << mc           << endl;  
  cout << "trigger      = " << trigger      << endl;
  cout << "foundReco1   = " << foundReco1   << endl;
  cout << "tagProbeTot1 = " << tagProbeTot1 << endl;
  cout << "probeIsol1   = " << probeIsol1   << endl;  
  cout << "foundReco2   = " << foundReco2   << endl;
  cout << "tagProbeTot2 = " << tagProbeTot2 << endl;
  cout << "probeIsol2   = " << probeIsol2   << endl;  

  cout << "statistics from MC - electrons: " << endl;
  cout << "eleEta1      = " << eleEta1      << endl;
  cout << "elePt1       = " << elePt1       << endl;
  cout << "eleEta2      = " << eleEta2      << endl;
  cout << "elePt2       = " << elePt2       << endl;

  reducedTree.save();
}

void LHPdfsProducer::LoopQCD() {
  
  if(fChain == 0) return;
  
  bookHistos();
  
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;


    
    // trigger
    Utils anaUtils;
    // bool passedHLT = anaUtils.getTriggersOR(m_requiredBackgroundTriggers, firedTrg);
    // if(!passedHLT) continue;

    // fill the PDFs for QCD with all the (isolated) reco'ed electrons
    for(int iele=0;iele<nEle;iele++) {
      
      if( m_selection->getSwitch("etaEleAcc") && 
          ! m_selection->passCut("etaEleAcc",etaEle[iele]) ) continue;
      
      if ( m_selection->getSwitch("applyIsolationOnProbe") &&
           ! m_selection->passCut("relSumPtTracks",dr04TkSumPtEle[iele]) ) continue;
      
      TLorentzVector eleP4(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);

      /// define the bins in which can be splitted the PDFs
      int iecal  = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = classificationEle[iele];
        
      int iclass = -1;
      int ifullclass = -1;
      if ( fullclassRaw == 0 || fullclassRaw == 100 ) { // golden
        iclass = 0;
        ifullclass = 0;
      }
      else if ( fullclassRaw == 10 || fullclassRaw == 110 ) { // bigbrem
        iclass = 0;
        ifullclass = 1;
      }
      else if ( fullclassRaw == 20 || fullclassRaw == 120 ) { // narrow
        iclass = 0;
        ifullclass = 2;
      }
      else if ( (fullclassRaw >= 30 && fullclassRaw <= 40) ||
                (fullclassRaw >= 130 && fullclassRaw <= 140) ) { // showering + cracks
        iclass = 1;
        ifullclass = 3;
      }
        
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[iele]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiEle[iele]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiEle[iele]));
      float s1s9   = s1s9Ele[iele];
      float s9s25  = s9s25Ele[iele];
      float lat    = latEle[iele];
      float etaLat = etaLatEle[iele];
      float phiLat = phiLatEle[iele];
      float a20    = a20Ele[iele];
      float a42    = a42Ele[iele];
      double dPhiCalo = deltaPhiAtCaloEle[iele];
      double dPhiVtx  = deltaPhiAtVtxEle[iele];
      double dEtaVtx  = deltaEtaAtVtxEle[iele];
      double EoPout   = eSeedOverPoutEle[iele];
      double EoP      = eSuperClusterOverPEle[iele];
      double HoE = hOverEEle[iele];

      dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
      EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );
      LATUnsplitEle           [iecal][iptbin] -> Fill ( lat );
      etaLATUnsplitEle        [iecal][iptbin] -> Fill ( etaLat );
      phiLATUnsplitEle        [iecal][iptbin] -> Fill ( phiLat );
      a20UnsplitEle           [iecal][iptbin] -> Fill ( a20 );
      a42UnsplitEle           [iecal][iptbin] -> Fill ( a42 );

      dPhiCaloClassEle        [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle         [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle            [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
      EoPoutClassEle          [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle             [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle            [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle           [iecal][iptbin][iclass] -> Fill ( s9s25 );
      LATClassEle             [iecal][iptbin][iclass] -> Fill ( lat );
      etaLATClassEle          [iecal][iptbin][iclass] -> Fill ( etaLat );
      phiLATClassEle          [iecal][iptbin][iclass] -> Fill ( phiLat );
      a20ClassEle             [iecal][iptbin][iclass] -> Fill ( a20 );
      a42ClassEle             [iecal][iptbin][iclass] -> Fill ( a42 );

      dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
      EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
      LATFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( lat );
      etaLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( etaLat );
      phiLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( phiLat );
      a20FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a20 );
      a42FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a42 );

    } // loop on electrons

  } // loop on events
  
}


void LHPdfsProducer::LoopQCDTagAndProbe(const char *treefilesuffix) {

  if(fChain == 0) return;
  
  bookHistos();
  
  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesBackground();
  reducedTree.addCategories();

  // counters
  int allevents   = 0;
  int trigger     = 0;
  int oneele      = 0;
  int onejet      = 0;
  int tagandprobe = 0;
  int nocrack     = 0;
  int eleTot      = 0;
  int deltaphi    = 0;
  int invmass     = 0;

  // loop over events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    allevents++;
    
    // QCD trigger
    Utils anaUtilsQCD;
    bool passedHLTQCD = anaUtilsQCD.getTriggersOR(m_requiredBackgroundTriggers, firedTrg);
    if ( m_selection->getSwitch("requireTriggerQCDBack") && !passedHLTQCD ) continue;   
    trigger++;

    // electrons and jets within acceptance
    vector<int> probeCandidates;
    vector<int> tagCandidates;
    
    // selecting electrons in the acceptance and possibly loose isolated to work as probe
    for(int iele=0; iele<nEle; iele++) {
      TLorentzVector electron(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);
      float relativeIsol = dr04TkSumPtEle[iele]/electron.Pt();
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[iele]) ) ) continue;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc",electron.Pt()) ) ) continue;      
      if( m_selection->getSwitch("applyIsolationOnProbe") && !m_selection->passCut("relSumPtTracks",relativeIsol) ) continue;
      probeCandidates.push_back(iele);
    }
    if (probeCandidates.size()>0) oneele++;

    // selecting jets in the acceptance to work as tag
    for (int ijet=0; ijet<nSisConeCorrJet; ijet++){    
      if ( m_selection->getSwitch("etaJetAcc") && !m_selection->passCut("etaJetAcc",etaSisConeCorrJet[ijet]) ) continue;
      if ( m_selection->getSwitch("etJetAcc")  && !m_selection->passCut("etJetAcc",etSisConeCorrJet[ijet]) )   continue;
      tagCandidates.push_back(ijet);
    }
    if (tagCandidates.size()>0) onejet++;


    // to select the highest Et jet in the event within the acceptance in eta (tag) matching a probe electron
    int theTag   = -999;
    int theProbe = -999;
    float tagEt  = -999.;
    
    for (int iJet=0; iJet<tagCandidates.size(); iJet++){ 

      TLorentzVector p4Jet;
      TVector3 p3Jet(pxSisConeCorrJet[tagCandidates[iJet]],pySisConeCorrJet[tagCandidates[iJet]],pzSisConeCorrJet[tagCandidates[iJet]]);
      p4Jet.SetXYZT (pxSisConeCorrJet[tagCandidates[iJet]],pySisConeCorrJet[tagCandidates[iJet]],pzSisConeCorrJet[tagCandidates[iJet]], energySisConeCorrJet[tagCandidates[iJet]]);

      // we use the highest ET jet (with a potential probe) as a tag 
      if (etSisConeCorrJet[tagCandidates[iJet]]<tagEt) continue;

      // we choose the probe with max Dphi wrt the tag
      float dPhiMax = -999.;

      for (int iEle=0; iEle<probeCandidates.size(); iEle++){       

	eleTot++;

	TLorentzVector p4Ele;	
	TVector3 p3Ele(pxEle[probeCandidates[iEle]],pyEle[probeCandidates[iEle]],pzEle[probeCandidates[iEle]]);
	p4Ele.SetXYZT (pxEle[probeCandidates[iEle]],pyEle[probeCandidates[iEle]],pzEle[probeCandidates[iEle]],energyEle[probeCandidates[iEle]]);

	// the electron and the tag must not be the same thing
	float deltaR = p3Ele.DeltaR(p3Jet);
	if (deltaR<0.3) continue;

	// minimal requirements on invariant mass and separation
	float deltaPhi = fabs(p3Jet.DeltaPhi(p3Ele));
	float invMass  = (p4Jet+p4Ele).M();
	if ( m_selection->getSwitch("jetDeltaPhi") && !m_selection->passCut("jetDeltaPhi", deltaPhi) ) continue;
	deltaphi++;
	if ( m_selection->getSwitch("jetInvMass")  && m_selection->passCut("jetInvMass", invMass) ) continue;
	invmass++;

	// we have a potential probe and this is highest ET jet up to now -> it's our tag
	tagEt  = etSisConeCorrJet[tagCandidates[iJet]];
	theTag = tagCandidates[iJet];
	
	// in case of several probes for this tag
	if (deltaPhi>dPhiMax) { 
	  dPhiMax  = deltaPhi;
	  theProbe = probeCandidates[iEle];
	}
      }
    }
  
    // we need a tag and a probe
    if (theTag<-800 || theProbe<-800) continue;
    tagandprobe++;
        
    // variables for the tree
    TLorentzVector p4Tag, p4Probe;	
    TVector3 p3Probe(pxEle[theProbe],pyEle[theProbe],pzEle[theProbe]);
    p4Probe.SetXYZT (pxEle[theProbe],pyEle[theProbe],pzEle[theProbe],energyEle[theProbe]);
    TVector3 p3Tag(pxSisConeCorrJet[theTag],pySisConeCorrJet[theTag],pzSisConeCorrJet[theTag]);
    p4Tag.SetXYZT (pxSisConeCorrJet[theTag],pySisConeCorrJet[theTag],pzSisConeCorrJet[theTag], energySisConeCorrJet[theTag]);
    float theDeltaPhi = fabs(p3Tag.DeltaPhi(p3Probe));
    float theInvMass  = (p4Tag+p4Probe).M();
    float theMet      = etMet[0];

    if (theDeltaPhi<0.5) { 
      cout << "probe: " << etaEle[theProbe]          << " " << phiEle[theProbe]          << " " << etEle[theProbe] << endl;
      cout << "tag "    << etaSisConeCorrJet[theTag] << " " << phiSisConeCorrJet[theTag] << " " << etSisConeCorrJet[theTag] << endl;
      cout << "dPhi = " << theDeltaPhi << endl;
      cout << jentry << endl;
    }

    // others
    int charge = chargeEle[theProbe];
    float pt   = p4Probe.Pt();
    float eta  = etaEle[theProbe];

    // fill the PDFs for QCD with the probe
    TLorentzVector eleP4(pxEle[theProbe],pyEle[theProbe],pzEle[theProbe],energyEle[theProbe]);

    /// define the bins in which can be splitted the PDFs
    int iecal  = (fabs( etaEle[theProbe])<1.479) ? 0 : 1;
    int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
    int fullclassRaw = classificationEle[theProbe];
    int iclass = -1;
    int ifullclass = -1;    
    if      ( fullclassRaw == GOLDEN ) { iclass = 0; ifullclass = 0; }
    else if ( fullclassRaw == BIGBREM ) { iclass = 0; ifullclass = 1; }
    else if ( fullclassRaw == NARROW ) { iclass = 0; ifullclass = 2; }
    else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
    if (iclass>-1) nocrack++;
          
    if (iclass>-1) {
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[theProbe]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiEle[theProbe]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiEle[theProbe]));
      float s1s9      = s1s9Ele[theProbe];
      float s9s25     = s9s25Ele[theProbe];
      float lat       = latEle[theProbe];
      float etaLat    = etaLatEle[theProbe];
      float phiLat    = phiLatEle[theProbe];
      float a20       = a20Ele[theProbe];
      float a42       = a42Ele[theProbe];
      double dPhiCalo = deltaPhiAtCaloEle[theProbe];
      double dPhiVtx  = deltaPhiAtVtxEle[theProbe];
      double dEtaVtx  = deltaEtaAtVtxEle[theProbe];
      double EoPout   = eSeedOverPoutEle[theProbe];
      double EoP      = eSuperClusterOverPEle[theProbe];
      double HoE      = hOverEEle[theProbe];
      
      dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
      EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );
      LATUnsplitEle           [iecal][iptbin] -> Fill ( lat );
      etaLATUnsplitEle        [iecal][iptbin] -> Fill ( etaLat );
      phiLATUnsplitEle        [iecal][iptbin] -> Fill ( phiLat );
      a20UnsplitEle           [iecal][iptbin] -> Fill ( a20 );
      a42UnsplitEle           [iecal][iptbin] -> Fill ( a42 );
      
      dPhiCaloClassEle        [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle         [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle            [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
      EoPoutClassEle          [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle             [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle            [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle           [iecal][iptbin][iclass] -> Fill ( s9s25 );
      LATClassEle             [iecal][iptbin][iclass] -> Fill ( lat );
      etaLATClassEle          [iecal][iptbin][iclass] -> Fill ( etaLat );
      phiLATClassEle          [iecal][iptbin][iclass] -> Fill ( phiLat );
      a20ClassEle             [iecal][iptbin][iclass] -> Fill ( a20 );
      a42ClassEle             [iecal][iptbin][iclass] -> Fill ( a42 );
      
      dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
      EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
      LATFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( lat );
      etaLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( etaLat );
      phiLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( phiLat );
      a20FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a20 );
      a42FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a42 );

      // fill the reduced tree
      reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,s1s9,sigmaIEtaIEta);    
      reducedTree.fillAttributesBackground(charge,eta,pt,theDeltaPhi,theInvMass,theMet);
      reducedTree.fillCategories(iecal,iptbin,iclass);
      reducedTree.store();
    } // no showering    

  } // loop on events

  // statistics
  cout << "statistics from QCD tag and probe: " << endl;
  cout << "allevents      = " << allevents      << endl;
  cout << "trigger        = " << trigger        << endl;
  cout << "one probe cand = " << oneele         << endl;
  cout << "one tag cand   = " << onejet         << endl;
  cout << "T and P found  = " << tagandprobe    << endl; 
  cout << "no crack       = " << nocrack        << endl;
  cout << "statistics from QCD tag and probe - electrons: " << endl;
  cout << "eleTot         = " << eleTot         << endl;
  cout << "deltaPhi ok    = " << deltaphi       << endl;
  cout << "inv mass       = " << invmass        << endl; 

  reducedTree.save();  
}


// make jet PDFs from W+jets 
void LHPdfsProducer::LoopWjets() {

  if(fChain == 0) return;

  bookHistos();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // trigger
    // Utils anaUtils;
    // bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
    // if(!passedHLT) continue;
    
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

    TVector3 mcEle(0,0,0);
    if(mceleindex>-1) mcEle = TVector3(pMc[mceleindex]*cos(phiMc[mceleindex])*sin(thetaMc[mceleindex]), 
                                       pMc[mceleindex]*sin(phiMc[mceleindex])*sin(thetaMc[mceleindex]),
                                       pMc[mceleindex]*cos(thetaMc[mceleindex]));

    // fill the PDFs for jets excluding the real electron in W+jets
    for(int iele=0;iele<nEle;iele++) {

      TLorentzVector eleP4(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);
      
      float deltaR = 1000;
      if(mceleindex>-1) deltaR = eleP4.Vect().DeltaR(mcEle);      
      if(deltaR<0.3) continue;
      
      if( m_selection->getSwitch("etaEleAcc") && 
          ! m_selection->passCut("etaEleAcc",etaEle[iele]) ) continue;
      
      if ( m_selection->getSwitch("applyIsolationOnProbe") &&
           ! m_selection->passCut("relSumPtTracks",dr04TkSumPtEle[iele]) ) continue;
      
      
      // define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = classificationEle[iele];
        
      int iclass = -1;
      int ifullclass = -1;
      if ( fullclassRaw == 0 || fullclassRaw == 100 ) { // golden
        iclass = 0;
        ifullclass = 0;
      }
      else if ( fullclassRaw == 10 || fullclassRaw == 110 ) { // bigbrem
        iclass = 0;
        ifullclass = 1;
      }
      else if ( fullclassRaw == 20 || fullclassRaw == 120 ) { // narrow
        iclass = 0;
        ifullclass = 2;
      }
      else if ( (fullclassRaw >= 30 && fullclassRaw <= 40) ||
                (fullclassRaw >= 130 && fullclassRaw <= 140) ) { // showering + cracks
        iclass = 1;
        ifullclass = 3;
      }

      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[iele]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiEle[iele]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiEle[iele]));
      float s1s9   = s1s9Ele[iele];
      float s9s25  = s9s25Ele[iele];
      float lat    = latEle[iele];
      float etaLat = etaLatEle[iele];
      float phiLat = phiLatEle[iele];
      float a20 = a20Ele[iele];
      float a42 = a42Ele[iele];

      double dPhiCalo = deltaPhiAtCaloEle[iele];
      double dPhiVtx  = deltaPhiAtVtxEle[iele];
      double dEta     = deltaEtaAtVtxEle[iele];
      double EoPout   = eSeedOverPoutEle[iele];
      double EoP      = eSuperClusterOverPEle[iele];
      double HoE      = hOverEEle[iele];
      double dxy      = eleTrackDxyEle[iele];
      double dxySig   = eleTrackDxyEle[iele]/eleTrackDxyErrorEle[iele];
      
      dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEta );
      EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );
      LATUnsplitEle           [iecal][iptbin] -> Fill ( lat );
      etaLATUnsplitEle        [iecal][iptbin] -> Fill ( etaLat );
      phiLATUnsplitEle        [iecal][iptbin] -> Fill ( phiLat );
      a20UnsplitEle           [iecal][iptbin] -> Fill ( a20 );
      a42UnsplitEle           [iecal][iptbin] -> Fill ( a42 );
      dxyUnsplitEle           [iecal][iptbin] -> Fill ( dxy );
      dxySigUnsplitEle        [iecal][iptbin] -> Fill ( dxySig );

      dPhiCaloClassEle      [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle       [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle          [iecal][iptbin][iclass] -> Fill ( dEta );
      EoPoutClassEle        [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle           [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle          [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle         [iecal][iptbin][iclass] -> Fill ( s9s25 );
      LATClassEle           [iecal][iptbin][iclass] -> Fill ( lat );
      etaLATClassEle        [iecal][iptbin][iclass] -> Fill ( etaLat );
      phiLATClassEle        [iecal][iptbin][iclass] -> Fill ( phiLat );
      a20ClassEle           [iecal][iptbin][iclass] -> Fill ( a20 );
      a42ClassEle           [iecal][iptbin][iclass] -> Fill ( a42 );
      dxyClassEle           [iecal][iptbin][iclass] -> Fill ( dxy );
      dxySigClassEle        [iecal][iptbin][iclass] -> Fill ( dxySig );

      dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEta );
      EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
      LATFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( lat );
      etaLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( etaLat );
      phiLATFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( phiLat );
      a20FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a20 );
      a42FullclassEle           [iecal][iptbin][ifullclass] -> Fill ( a42 );
      dxyFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( dxy );
      dxySigFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( dxySig );

    } // loop on electrons

  } // loop on events
  
}


// make jets PDF removing the 2 electrons from Z: TOCHECK!!! something strange!
void LHPdfsProducer::LoopZjets(const char *outname) {

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
  
  TH1F *AllFakesEta = new TH1F("AllFakesEta", "all reconstructed fakes", nbinsEta, minEta, maxEta);
  TH1F *UnmatchedFakesEta = new TH1F("UnmatchedFakesEta", "Fake electrons unmatched with true electrons", nbinsEta, minEta, maxEta);

  int nbinsPt = 20;
  float minPt = 10.0;
  float maxPt = 100.;

  TH1F *AllFakesPt = new TH1F("AllFakesPt", "all reconstructed fakes", nbinsPt, minPt, maxPt);
  TH1F *UnmatchedFakesPt = new TH1F("UnmatchedFakesPt", "Fake electrons unmatched with true electrons", nbinsPt, minPt, maxPt);

  if(fChain == 0) return;

  bookHistos();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // trigger
    // Utils anaUtils;
    // bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
    // if(!passedHLT) continue;

    bool tauPresence=false;
    bool spuriousElectronsMadgraph=false;
    int mcEleMinusIndex=-1, mcElePlusIndex=-1;
    for(int imc=0; imc<50; imc++) {
      // not only ele from Zee: there is enu emission (V_ud?) in madgraph: remove these spurious electrons
      if ( (fabs(idMc[imc])==11) && fabs(idMc[mothMc[imc]])!=23 && fabs(idMc[mothMc[imc]])!=11 ) {
        spuriousElectronsMadgraph=true;
        break;
      }
      if ( idMc[imc]==11 && fabs(idMc[mothMc[imc]])==23 ) mcEleMinusIndex = imc;
      if ( idMc[imc]==-11 && fabs(idMc[mothMc[imc]])==23 ) mcElePlusIndex = imc;
      // since the stable particle list is truncated, if there is a tau 
      // not possible to say what happens...
      if ( (fabs(idMc[imc])==15) ) {
        tauPresence=true;
        break;
      }
    }

    if(tauPresence || spuriousElectronsMadgraph) continue;

    TVector3 mcElePlus(0,0,0), mcEleMinus(0,0,0);
    if(mcEleMinusIndex>-1) mcEleMinus = TVector3(pMc[mcEleMinusIndex]*cos(phiMc[mcEleMinusIndex])*sin(thetaMc[mcEleMinusIndex]), 
                                                 pMc[mcEleMinusIndex]*sin(phiMc[mcEleMinusIndex])*sin(thetaMc[mcEleMinusIndex]),
                                                 pMc[mcEleMinusIndex]*cos(thetaMc[mcEleMinusIndex]));
    
    if(mcElePlusIndex>-1) mcEleMinus = TVector3(pMc[mcElePlusIndex]*cos(phiMc[mcElePlusIndex])*sin(thetaMc[mcElePlusIndex]), 
                                                pMc[mcElePlusIndex]*sin(phiMc[mcElePlusIndex])*sin(thetaMc[mcElePlusIndex]),
                                                pMc[mcElePlusIndex]*cos(thetaMc[mcElePlusIndex]));
    

    // fill the PDFs for jets excluding the real electron of Z->ee
    // the two electrons are the nes that give the |mee-mZ| closest to 0
    float minpullZee = 1000.;
    electrons[0]=-1; electrons[1]=-1;
    if(nEle>=2) {
      for(int iele1=0; iele1<nEle; iele1++) {
        TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);
        for(int iele2=iele1+1; iele2<nEle; iele2++) {
          TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);
          
          float mass = (electron1+electron2).M();
          m_Zmass->Fill(mass);
          float pull=fabs(mass-91.1876);
          
          if(pull < minpullZee) {
            minpullZee = pull;
            electrons[0] = iele1;
            electrons[1] = iele2;
          }
        }
      }
    }

    float minpullZmumu = 1000.;
    muons[0]=-1, muons[1]=-1;
    if(nMuon>=2) { 
      for(int imu1=0; imu1<nMuon; imu1++) {
        TLorentzVector muon1(pxMuon[imu1],pyMuon[imu1],pzMuon[imu1],momentumMuon[imu1]);
        for(int imu2=imu1+1; imu2<nMuon; imu2++) {
          TLorentzVector muon2(pxMuon[imu2],pyMuon[imu2],pzMuon[imu2],momentumMuon[imu2]);
          
          float mass = (muon1+muon2).M();
          float pull=fabs(mass-91.1876);
          
          if(pull < minpullZmumu) {
            minpullZmumu = pull;
            muons[0] = imu1;
            muons[1] = imu2;
          }
        }
      }
    }

    bool Ztag = false;
    // require the reconstruction of the Z
    // the two electrons in the acceptance, pTmin and Z mass
    // or the same with muons
    if (electrons[0]=!-1 && electrons[1]!=-1) {

      if( minpullZee < 9.0 &&
          fabs(etaEle[electrons[0]]) < 2.5 &&
          fabs(etaEle[electrons[1]]) < 2.5 &&
          etEle[electrons[0]] > 5.0 &&
          etEle[electrons[1]] > 5.0 ) Ztag = true;

    } else if (muons[0]!=-1 && muons[1]!=-1) {
      
      if( minpullZmumu < 9.0 &&
          fabs(etaMuon[muons[0]]) < 2.4 &&
          fabs(etaMuon[muons[1]]) < 2.4 &&
          etMuon[muons[0]] > 5.0 &&
          etMuon[muons[1]] > 5.0 ) Ztag = true;

    }


    if(!Ztag) continue;
    
    cout << "Looking for fakes" << endl;

    // find fake electrons
    for(int iele=0;iele<nEle;iele++) {

      cout << "iele = " << iele << "  electrons[0] = " << electrons[0] << "   electrons[1] = " << electrons[1] << endl;

      if( iele==electrons[0] || iele==electrons[1] ) continue;

      if( fabs(etaEle[iele]) > 2.5 || etEle[iele] < 5.0 ) continue;
      
      AllFakesEta->Fill(etaEle[iele]);
      AllFakesPt->Fill(etEle[iele]);

      TVector3 eleP3(pxEle[iele],pyEle[iele],pzEle[iele]);
      
      // evaluate the purity of the jets: p=#ele(unmatched)/#ele
      if(mcEleMinusIndex>-1 && mcElePlusIndex>-1) {
        float dRPlus = mcElePlus.DeltaR(eleP3);
        float dRMinus = mcEleMinus.DeltaR(eleP3);
        if(dRPlus>0.3 && dRMinus>0.3) {
          UnmatchedFakesEta->Fill(etaEle[iele]);
          UnmatchedFakesPt->Fill(etEle[iele]);
        }
      }
      
      /// define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (etEle[iele]<15.0) ? 0 : 1;
      
      int fullclassRaw = classificationEle[iele];
        
      int iclass = -1;
      int ifullclass = -1;
      if ( fullclassRaw == 0 || fullclassRaw == 100 ) { // golden
        iclass = 0;
        ifullclass = 0;
      }
      else if ( fullclassRaw == 10 || fullclassRaw == 110 ) { // bigbrem
        iclass = 0;
        ifullclass = 1;
      }
      else if ( fullclassRaw == 20 || fullclassRaw == 120 ) { // narrow
        iclass = 0;
        ifullclass = 2;
      }
      else if ( (fullclassRaw >= 30 && fullclassRaw <= 40) ||
                (fullclassRaw >= 130 && fullclassRaw <= 140) ) { // showering + cracks
        iclass = 1;
        ifullclass = 3;
      }

      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaEle[iele]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiEle[iele]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiEle[iele]));
      float s1s9 = s1s9Ele[iele];
      float s9s25 = s9s25Ele[iele];
      float lat = latEle[iele];
      float etaLat = etaLatEle[iele];
      float phiLat = phiLatEle[iele];
      float a20 = a20Ele[iele];
      float a42 = a42Ele[iele];

      double dPhiCalo = deltaPhiAtCaloEle[iele];
      double dPhiVtx = deltaPhiAtVtxEle[iele];
      double dEta = deltaEtaAtVtxEle[iele];
      double EoPout = eSeedOverPoutEle[iele];
      double EoP = eSuperClusterOverPEle[iele];
      double HoE = hOverEEle[iele];
      double dxy = eleTrackDxyEle[iele];
      double dxySig = eleTrackDxyEle[iele]/eleTrackDxyErrorEle[iele];

      dPhiCaloUnsplitEle    [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle     [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle        [iecal][iptbin] -> Fill ( dEta );
      EoPoutUnsplitEle      [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle         [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle        [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle       [iecal][iptbin] -> Fill ( s9s25 );
      LATUnsplitEle         [iecal][iptbin] -> Fill ( lat );
      etaLATUnsplitEle      [iecal][iptbin] -> Fill ( etaLat );
      phiLATUnsplitEle      [iecal][iptbin] -> Fill ( phiLat );
      a20UnsplitEle         [iecal][iptbin] -> Fill ( a20 );
      a42UnsplitEle         [iecal][iptbin] -> Fill ( a42 );
      dxyUnsplitEle         [iecal][iptbin] -> Fill ( dxy );
      dxySigUnsplitEle      [iecal][iptbin] -> Fill ( dxySig );

      dPhiCaloClassEle    [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle     [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle        [iecal][iptbin][iclass] -> Fill ( dEta );
      EoPoutClassEle      [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle         [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle        [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle       [iecal][iptbin][iclass] -> Fill ( s9s25 );
      LATClassEle         [iecal][iptbin][iclass] -> Fill ( lat );
      etaLATClassEle      [iecal][iptbin][iclass] -> Fill ( etaLat );
      phiLATClassEle      [iecal][iptbin][iclass] -> Fill ( phiLat );
      a20ClassEle         [iecal][iptbin][iclass] -> Fill ( a20 );
      a42ClassEle         [iecal][iptbin][iclass] -> Fill ( a42 );
      dxyClassEle         [iecal][iptbin][iclass] -> Fill ( dxy );
      dxySigClassEle      [iecal][iptbin][iclass] -> Fill ( dxySig );

      dPhiCaloFullclassEle    [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle     [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( dEta );
      EoPoutFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle        [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle       [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
      LATFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( lat );
      etaLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( etaLat );
      phiLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( phiLat );
      a20FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a20 );
      a42FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a42 );
      dxyFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( dxy );
      dxySigFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dxySig );

    } // loop on electrons

  } // loop on events

  char filename[200];
  sprintf(filename,"%s-JetPurityEta.root",outname);
  EfficiencyEvaluator PurityOfFakesEta(filename);
  PurityOfFakesEta.AddNumerator(AllFakesEta);
  PurityOfFakesEta.AddNumerator(UnmatchedFakesEta);
  PurityOfFakesEta.SetDenominator(AllFakesEta);
  PurityOfFakesEta.ComputeEfficiencies();
  PurityOfFakesEta.Write();

  sprintf(filename,"%s-JetPurityPt.root",outname);
  EfficiencyEvaluator PurityOfFakesPt(filename);
  PurityOfFakesPt.AddNumerator(AllFakesPt);
  PurityOfFakesPt.AddNumerator(UnmatchedFakesPt);
  PurityOfFakesPt.SetDenominator(AllFakesPt);
  PurityOfFakesPt.ComputeEfficiencies();
  PurityOfFakesPt.Write();
  
}


void LHPdfsProducer::bookHistos() {

  m_Zmass = new TH1F("Zmass", "Zmass", 260, 0., 130.);
  
  int nbins = 100;

  float dPhiCaloMin = -0.3;
  float dPhiCaloMax =  0.3;
  float dPhiVtxMin  = -0.1; // pixelMatchGsfElectron pre-selection: |dPhi| < 0.1
  float dPhiVtxMax  =  0.1;
  float dEtaMin     = -0.02;
  float dEtaMax     =  0.02;
  float EoPoutMin   =  0.0;
  float EoPoutMax   =  8.0;
  float HoEMin      =  0.0; // zero-suppression in HCAL
  float HoEMax      =  0.1; // ??
  float sigmaIEtaIEtaMin = 0.0;
  float sigmaIEtaIEtaMax = 0.05;
  float sigmaIEtaIPhiMin = 0.0;
  float sigmaIEtaIPhiMax = 0.05;
  float sigmaIPhiIPhiMin = 0.0;
  float sigmaIPhiIPhiMax = 0.05;
  float s1s9Min  = 0.0;
  float s1s9Max  = 1.0;
  float s9s25Min = 0.5;
  float s9s25Max = 1.0;
  float LATMin   = 0.0;
  float LATMax   = 1.0;
  float etaLATMin = 0.0;
  float etaLATMax = 1.0;
  float phiLATMin = 0.0;
  float phiLATMax = 1.0;
  float a20Min    = 0.0;
  float a20Max    = 1.0;
  float a42Min    = 0.0;
  float a42Max    = 1.0;
  float dxyMin    = -0.04;
  float dxyMax    = 0.04;
  float dxySigMin    = -10.0;
  float dxySigMax    = 10.0;

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {

    // iptbin = 0: < 15 GeV
    // iptbin = 1: > 15 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      char histo[200];
      
      sprintf(histo,"dPhiCaloUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiCaloUnsplitEle[iecal][iptbin]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
      sprintf(histo,"dPhiVtxUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiVtxUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
      sprintf(histo,"dEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      dEtaUnsplitEle[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      sprintf(histo,"EoPoutUnsplit_electrons_%d_%d",iecal,iptbin);
      EoPoutUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
      sprintf(histo,"HoEUnsplit_electrons_%d_%d",iecal,iptbin);
      HoEUnsplitEle[iecal][iptbin]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      sprintf(histo,"sigmaIEtaIEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
      sprintf(histo,"sigmaIEtaIPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaIEtaIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIPhiMin, sigmaIEtaIPhiMax);
      sprintf(histo,"sigmaIPhiIPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiMin, sigmaIPhiIPhiMax);
      sprintf(histo,"s1s9Unsplit_electrons_%d_%d",iecal,iptbin);
      s1s9UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
      sprintf(histo,"s9s25Unsplit_electrons_%d_%d",iecal,iptbin);
      s9s25UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
      sprintf(histo,"LATUnsplit_electrons_%d_%d",iecal,iptbin);
      LATUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, LATMin, LATMax);
      sprintf(histo,"etaLATUnsplit_electrons_%d_%d",iecal,iptbin);
      etaLATUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
      sprintf(histo,"phiLATUnsplit_electrons_%d_%d",iecal,iptbin);
      phiLATUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
      sprintf(histo,"a20Unsplit_electrons_%d_%d",iecal,iptbin);
      a20UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, a20Min, a20Max);
      sprintf(histo,"a42Unsplit_electrons_%d_%d",iecal,iptbin);
      a42UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, a42Min, a42Max);
      sprintf(histo,"dxyUnsplit_electrons_%d_%d",iecal,iptbin);
      dxyUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, dxyMin, dxyMax);
      sprintf(histo,"dxySigUnsplit_electrons_%d_%d",iecal,iptbin);
      dxySigUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, dxySigMin, dxySigMax);

      // iclass = 0: non-showering
      // iclass = 1: showering
      for(int iclass=0; iclass<2; iclass++) {
      
	sprintf(histo,"dPhiCaloClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiCaloClassEle[iecal][iptbin][iclass]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
	sprintf(histo,"dPhiVtxClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiVtxClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
	sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	EoPoutClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
	sprintf(histo,"sigmaIEtaIPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIEtaIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIPhiMin, sigmaIEtaIPhiMax);
	sprintf(histo,"sigmaIPhiIPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiMin, sigmaIPhiIPhiMax);
	sprintf(histo,"s1s9Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s1s9ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s9s25ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
	sprintf(histo,"LATClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	LATClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, LATMin, LATMax);
	sprintf(histo,"etaLATClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	etaLATClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
	sprintf(histo,"phiLATClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	phiLATClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
	sprintf(histo,"a20Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	a20ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, a20Min, a20Max);
	sprintf(histo,"a42Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	a42ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, a42Min, a42Max);
        sprintf(histo,"dxyClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
        dxyClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, dxyMin, dxyMax);
        sprintf(histo,"dxySigClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
        dxySigClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, dxySigMin, dxySigMax);

      }

      // iclass = 0: golden
      // iclass = 1: bigbrem
      // iclass = 2: narrow
      // iclass = 3: showering + cracks

      for(int ifullclass=0; ifullclass<4; ifullclass++) {
      
	sprintf(histo,"dPhiCaloFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dPhiCaloFullclassEle[iecal][iptbin][ifullclass]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
	sprintf(histo,"dPhiVtxFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dPhiVtxFullclassEle[iecal][iptbin][ifullclass]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
	sprintf(histo,"dEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dEtaFullclassEle[iecal][iptbin][ifullclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	EoPoutFullclassEle[iecal][iptbin][ifullclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"HoEFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	HoEFullclassEle[iecal][iptbin][ifullclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"sigmaIEtaIEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaIEtaIEtaFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
	sprintf(histo,"sigmaIEtaIPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaIEtaIPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaIEtaIPhiMin, sigmaIEtaIPhiMax);
	sprintf(histo,"sigmaIPhiIPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaIPhiIPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiMin, sigmaIPhiIPhiMax);
	sprintf(histo,"s1s9Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	s1s9FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	s9s25FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
	sprintf(histo,"LATFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	LATFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, LATMin, LATMax);
	sprintf(histo,"etaLATFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	etaLATFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
	sprintf(histo,"phiLATFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	phiLATFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
	sprintf(histo,"a20Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	a20FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, a20Min, a20Max);
	sprintf(histo,"a42Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	a42FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, a42Min, a42Max);
        sprintf(histo,"dxyFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
        dxyFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, dxyMin, dxyMax);
        sprintf(histo,"dxySigFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
        dxySigFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, dxySigMin, dxySigMax);

      }

    }

  }

}

void LHPdfsProducer::saveHistos(const char *filename) {

  char outfilename[200];
  sprintf(outfilename,"%s",filename);
  TFile *file = TFile::Open(outfilename,"recreate");
  file->mkdir("pdfsProducer","pdfs created from trees");
  file->cd("pdfsProducer");

  m_Zmass->Write();

  for (int iecal=0; iecal<2; iecal++) {

    for(int iptbin=0; iptbin<2; iptbin++) {
      
      dPhiCaloUnsplitEle[iecal][iptbin]->Write();
      dPhiVtxUnsplitEle[iecal][iptbin]->Write();
      dEtaUnsplitEle[iecal][iptbin]->Write();
      EoPoutUnsplitEle[iecal][iptbin]->Write();
      HoEUnsplitEle[iecal][iptbin]->Write();
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin]->Write();
      sigmaIEtaIPhiUnsplitEle[iecal][iptbin]->Write();
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin]->Write();
      s1s9UnsplitEle[iecal][iptbin]->Write();
      s9s25UnsplitEle[iecal][iptbin]->Write();
      LATUnsplitEle[iecal][iptbin]->Write();
      etaLATUnsplitEle[iecal][iptbin]->Write();
      phiLATUnsplitEle[iecal][iptbin]->Write();
      a20UnsplitEle[iecal][iptbin]->Write();
      a42UnsplitEle[iecal][iptbin]->Write();
      dxyUnsplitEle[iecal][iptbin]->Write();
      dxySigUnsplitEle[iecal][iptbin]->Write();

      for(int iclass=0; iclass<2; iclass++) {
      
	dPhiCaloClassEle[iecal][iptbin][iclass]->Write();
	dPhiVtxClassEle[iecal][iptbin][iclass]->Write();
	dEtaClassEle[iecal][iptbin][iclass]->Write();
	EoPoutClassEle[iecal][iptbin][iclass]->Write();
	HoEClassEle[iecal][iptbin][iclass]->Write();
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass]->Write();
	sigmaIEtaIPhiClassEle[iecal][iptbin][iclass]->Write();
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass]->Write();
	s1s9ClassEle[iecal][iptbin][iclass]->Write();
	s9s25ClassEle[iecal][iptbin][iclass]->Write();
	LATClassEle[iecal][iptbin][iclass]->Write();
	etaLATClassEle[iecal][iptbin][iclass]->Write();
	phiLATClassEle[iecal][iptbin][iclass]->Write();
	a20ClassEle[iecal][iptbin][iclass]->Write();
	a42ClassEle[iecal][iptbin][iclass]->Write();
	dxyClassEle[iecal][iptbin][iclass]->Write();
	dxySigClassEle[iecal][iptbin][iclass]->Write();

      }

      for(int ifullclass=0; ifullclass<4; ifullclass++) {
      
	dPhiCaloFullclassEle[iecal][iptbin][ifullclass]->Write();
	dPhiVtxFullclassEle[iecal][iptbin][ifullclass]->Write();
	dEtaFullclassEle[iecal][iptbin][ifullclass]->Write();
	EoPoutFullclassEle[iecal][iptbin][ifullclass]->Write();
	HoEFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaIEtaIEtaFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaIEtaIPhiFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaIPhiIPhiFullclassEle[iecal][iptbin][ifullclass]->Write();
	s1s9FullclassEle[iecal][iptbin][ifullclass]->Write();
	s9s25FullclassEle[iecal][iptbin][ifullclass]->Write();
	LATFullclassEle[iecal][iptbin][ifullclass]->Write();
	etaLATFullclassEle[iecal][iptbin][ifullclass]->Write();
	phiLATFullclassEle[iecal][iptbin][ifullclass]->Write();
	a20FullclassEle[iecal][iptbin][ifullclass]->Write();
	a42FullclassEle[iecal][iptbin][ifullclass]->Write();
	dxyFullclassEle[iecal][iptbin][ifullclass]->Write();
	dxySigFullclassEle[iecal][iptbin][ifullclass]->Write();

      }

    }

  }

  file->Close();
  
}

bool LHPdfsProducer::isEleID(int eleIndex) {

  TVector3 pTrkAtOuter(pxAtOuterEle[eleIndex],pyAtOuterEle[eleIndex],pzAtOuterEle[eleIndex]);
  EgammaCutBasedID.SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  EgammaCutBasedID.SetHOverE( hOverEEle[eleIndex] );
  EgammaCutBasedID.SetS9S25( s9s25Ele[eleIndex] );
  EgammaCutBasedID.SetDEta( deltaEtaAtVtxEle[eleIndex] );
  EgammaCutBasedID.SetDPhiIn( deltaPhiAtVtxEle[eleIndex] );
  EgammaCutBasedID.SetDPhiOut( deltaPhiAtCaloEle[eleIndex] );
  EgammaCutBasedID.SetInvEminusInvP( 1./ecalEle[eleIndex]-1./momentumEle[eleIndex] );
  EgammaCutBasedID.SetBremFraction( fabs(momentumEle[eleIndex]-pTrkAtOuter.Mag())/momentumEle[eleIndex] );
  EgammaCutBasedID.SetSigmaEtaEta( sqrt(covEtaEtaEle[eleIndex]) );
  EgammaCutBasedID.SetSigmaPhiPhi( sqrt(covPhiPhiEle[eleIndex]) );
  EgammaCutBasedID.SetEOverPout( eSeedOverPoutEle[eleIndex] );
  EgammaCutBasedID.SetEOverPin( eSuperClusterOverPEle[eleIndex] );
  EgammaCutBasedID.SetElectronClass ( classificationEle[eleIndex] );

  bool isIdentified = EgammaCutBasedID.output();

  return isIdentified;
}
