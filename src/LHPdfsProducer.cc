#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/LHPdfsProducer.hh"

LHPdfsProducer::LHPdfsProducer(TTree *tree)
  : EgammaBase(tree) {
  
  std::string fileCuts("/afs/cern.ch/user/e/emanuele/scratch0/Likelihood21X/OfflineAnalysis/EgammaAnalysisTools/config/LHPdfsProducer/cuts.txt");
  std::string fileSwitches("/afs/cern.ch/user/e/emanuele/scratch0/Likelihood21X/OfflineAnalysis/EgammaAnalysisTools/config/LHPdfsProducer/switches.txt");

  m_selection = new Selection(fileCuts,fileSwitches);
  m_selection->addSwitch("requireTrigger");
  m_selection->addSwitch("applyIsolationOnProbe");
  m_selection->addCut("meeWindow");
  m_selection->addCut("etaEleAcc");
  m_selection->addCut("ptEleAcc");
  m_selection->addCut("relSumPtTracks");
  m_selection->summary();

}

LHPdfsProducer::~LHPdfsProducer() { }

void LHPdfsProducer::LoopZTagAndProbe() {

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
    bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

    //    if(!passedHLT) continue;

    // best tag-probe pair = mee closest to Z mass
    float minpull = 1000.;

    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);
      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);

        if( m_selection->getSwitch("etaEleAcc") && 
            (!m_selection->passCut("etaEleAcc",etaEle[iele1]) ||
             !m_selection->passCut("etaEleAcc",etaEle[iele2]) ) ) continue;

        if( m_selection->getSwitch("ptEleAcc") && 
            (!m_selection->passCut("ptEleAcc",electron1.Pt()) ||
             !m_selection->passCut("ptEleAcc",electron2.Pt()) ) ) continue;

        float mass = (electron1+electron2).M();
        m_Zmass->Fill(mass);
        float pull=fabs(mass-91.1876);

        if(pull < minpull) {
          minpull = pull;
          electrons[0] = iele1;
          electrons[1] = iele2;
        }

      }

    }

    // start the tag & probe
    if( m_selection->passCut("meeWindow",minpull) ) {
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele=1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }

        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],pyEle[tag],pzEle[tag],energyEle[tag]);
        
        /// define the bins in which can be splitted the PDFs
        int iecal = (fabs( etaEle[probe])<1.479) ? 0 : 1;
        int iptbin = (probeP4.Pt()<15.0) ? 0 : 1;
      
        int fullclassRaw = eleClassEle[probe];
        
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
        
        // apply the electron ID loose on the tag electron        
        float tagIdentified = eleIdCutBasedEle[tag];
        
        // Tracker isolation
        // on the tag electron...
        
        bool tagIsolated = ( m_selection->passCut("relSumPtTracks",eleSumPt04Ele[tag]) );

        // on the probe electron...
        bool probeIsolated = true;
        if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
          probeIsolated = ( m_selection->passCut("relSumPtTracks",eleSumPt04Ele[probe]) );
        }
        
        float sigmaEtaEta = sqrt(fabs(covEtaEtaEle[probe]));
        float sigmaEtaPhi = sqrt(fabs(covEtaPhiEle[probe]));
        float sigmaPhiPhi = sqrt(fabs(covPhiPhiEle[probe]));
        float s1s9 = s1s9Ele[probe];
        float s9s25 = s9s25Ele[probe];
        float lat = latEle[probe];
        float etaLat = etaLatEle[probe];
        float phiLat = phiLatEle[probe];
        float a20 = a20Ele[probe];
        float a42 = a42Ele[probe];
        float fisher = -1000;

        if ( iecal==0 ) { // barrel
          if ( iptbin == 0 ) {  // low pt
            fisher = 0.693496 - 12.7018 * sigmaEtaEta + 1.23863 * s9s25 - 10.115 * etaLat;
          }
          else if ( iptbin == 1 ) {  // high pt
            fisher = 6.02184 - 49.2656 * sigmaEtaEta + 2.49634 * s9s25 - 30.1528 * etaLat;
          }
        }
        else if ( iecal == 1  ) { // endcap
          if ( iptbin == 0 ) {  // low pt
            fisher = -1.11814 - 5.3288 * sigmaEtaEta + 4.51575 * s9s25 - 6.47578 * etaLat;
          }
          else if ( iptbin == 1 ) {  // high pt
            fisher = 0.536351 - 11.7401 * sigmaEtaEta + 3.61809 * s9s25 - 9.3025 * etaLat;
          }
        }
        
        /// fill the electron ID pdfs only if:
        /// the tag is loose isolated and identified
        /// the probe is loose isolated
        if( tagIsolated && tagIdentified && probeIsolated ) {
          
          double dPhiCalo = eleDeltaPhiAtCaloEle[probe];
          double dPhiVtx = eleDeltaPhiAtVtxEle[probe];
          double dEta = eleDeltaEtaAtVtxEle[probe];
          double EoPout = eleCorrEoPoutEle[probe];
          double HoE = eleHoEEle[probe];

          dPhiCaloUnsplitEle    [iecal][iptbin] -> Fill ( dPhiCalo );
          dPhiVtxUnsplitEle     [iecal][iptbin] -> Fill ( dPhiVtx );
          dEtaUnsplitEle        [iecal][iptbin] -> Fill ( dEta );
          EoPoutUnsplitEle      [iecal][iptbin] -> Fill ( EoPout );
          HoEUnsplitEle         [iecal][iptbin] -> Fill ( HoE );
          shapeFisherUnsplitEle [iecal][iptbin] -> Fill ( fisher );
          sigmaEtaEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaEtaEta );
          sigmaEtaPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaEtaPhi );
          sigmaPhiPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaPhiPhi );
          s1s9UnsplitEle        [iecal][iptbin] -> Fill ( s1s9 );
          s9s25UnsplitEle       [iecal][iptbin] -> Fill ( s9s25 );
          LATUnsplitEle         [iecal][iptbin] -> Fill ( lat );
          etaLATUnsplitEle      [iecal][iptbin] -> Fill ( etaLat );
          phiLATUnsplitEle      [iecal][iptbin] -> Fill ( phiLat );
          a20UnsplitEle         [iecal][iptbin] -> Fill ( a20 );
          a42UnsplitEle         [iecal][iptbin] -> Fill ( a42 );


          dPhiCaloClassEle    [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
          dPhiVtxClassEle     [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
          dEtaClassEle        [iecal][iptbin][iclass] -> Fill ( dEta );
          EoPoutClassEle      [iecal][iptbin][iclass] -> Fill ( EoPout );
          HoEClassEle         [iecal][iptbin][iclass] -> Fill ( HoE );
          shapeFisherClassEle [iecal][iptbin][iclass] -> Fill ( fisher );
          sigmaEtaEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaEtaEta );
          sigmaEtaPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaEtaPhi );
          sigmaPhiPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaPhiPhi );
          s1s9ClassEle        [iecal][iptbin][iclass] -> Fill ( s1s9 );
          s9s25ClassEle       [iecal][iptbin][iclass] -> Fill ( s9s25 );
          LATClassEle         [iecal][iptbin][iclass] -> Fill ( lat );
          etaLATClassEle      [iecal][iptbin][iclass] -> Fill ( etaLat );
          phiLATClassEle      [iecal][iptbin][iclass] -> Fill ( phiLat );
          a20ClassEle         [iecal][iptbin][iclass] -> Fill ( a20 );
          a42ClassEle         [iecal][iptbin][iclass] -> Fill ( a42 );


          dPhiCaloFullclassEle    [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
          dPhiVtxFullclassEle     [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
          dEtaFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( dEta );
          EoPoutFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( EoPout );
          HoEFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( HoE );
          shapeFisherFullclassEle [iecal][iptbin][ifullclass] -> Fill ( fisher );
          sigmaEtaEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaEtaEta );
          sigmaEtaPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaEtaPhi );
          sigmaPhiPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaPhiPhi );
          s1s9FullclassEle        [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
          s9s25FullclassEle       [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
          LATFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( lat );
          etaLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( etaLat );
          phiLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( phiLat );
          a20FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a20 );
          a42FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a42 );

        } // fill histograms

      } // loop over the 2 Z electrons

    } // end tag and probe

  }

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
//     Utils anaUtils;
//     bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

//     if(!passedHLT) continue;

    // fill the PDFs for QCD with all the (isolated) reco'ed electrons
    for(int iele=0;iele<nEle;iele++) {

      if( m_selection->getSwitch("etaEleAcc") && 
          ! m_selection->passCut("etaEleAcc",etaEle[iele]) ) continue;
      
      if ( m_selection->getSwitch("applyIsolationOnProbe") &&
           ! m_selection->passCut("relSumPtTracks",eleSumPt04Ele[iele]) ) continue;
      
      TLorentzVector eleP4(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);

      /// define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = eleClassEle[iele];
        
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
        
      float sigmaEtaEta = sqrt(fabs(covEtaEtaEle[iele]));
      float sigmaEtaPhi = sqrt(fabs(covEtaPhiEle[iele]));
      float sigmaPhiPhi = sqrt(fabs(covPhiPhiEle[iele]));
      float s1s9 = s1s9Ele[iele];
      float s9s25 = s9s25Ele[iele];
      float lat = latEle[iele];
      float etaLat = etaLatEle[iele];
      float phiLat = phiLatEle[iele];
      float a20 = a20Ele[iele];
      float a42 = a42Ele[iele];
      float fisher = -1000;

      if ( iecal==0 ) { // barrel
        if ( iptbin == 0 ) {  // low pt
          fisher = 0.693496 - 12.7018 * sigmaEtaEta + 1.23863 * s9s25 - 10.115 * etaLat;
        }
        else if ( iptbin == 1 ) {  // high pt
          fisher = 6.02184 - 49.2656 * sigmaEtaEta + 2.49634 * s9s25 - 30.1528 * etaLat;
        }
      }
      else if ( iecal == 1  ) { // endcap
        if ( iptbin == 0 ) {  // low pt
          fisher = -1.11814 - 5.3288 * sigmaEtaEta + 4.51575 * s9s25 - 6.47578 * etaLat;
        }
        else if ( iptbin == 1 ) {  // high pt
          fisher = 0.536351 - 11.7401 * sigmaEtaEta + 3.61809 * s9s25 - 9.3025 * etaLat;
        }
      }

      double dPhiCalo = eleDeltaPhiAtCaloEle[iele];
      double dPhiVtx = eleDeltaPhiAtVtxEle[iele];
      double dEta = eleDeltaEtaAtVtxEle[iele];
      double EoPout = eleCorrEoPoutEle[iele];
      double HoE = eleHoEEle[iele];

      dPhiCaloUnsplitEle    [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle     [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle        [iecal][iptbin] -> Fill ( dEta );
      EoPoutUnsplitEle      [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle         [iecal][iptbin] -> Fill ( HoE );
      shapeFisherUnsplitEle [iecal][iptbin] -> Fill ( fisher );
      sigmaEtaEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaEtaEta );
      sigmaEtaPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaEtaPhi );
      sigmaPhiPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaPhiPhi );
      s1s9UnsplitEle        [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle       [iecal][iptbin] -> Fill ( s9s25 );
      LATUnsplitEle         [iecal][iptbin] -> Fill ( lat );
      etaLATUnsplitEle      [iecal][iptbin] -> Fill ( etaLat );
      phiLATUnsplitEle      [iecal][iptbin] -> Fill ( phiLat );
      a20UnsplitEle         [iecal][iptbin] -> Fill ( a20 );
      a42UnsplitEle         [iecal][iptbin] -> Fill ( a42 );


      dPhiCaloClassEle    [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle     [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle        [iecal][iptbin][iclass] -> Fill ( dEta );
      EoPoutClassEle      [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle         [iecal][iptbin][iclass] -> Fill ( HoE );
      shapeFisherClassEle [iecal][iptbin][iclass] -> Fill ( fisher );
      sigmaEtaEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaEtaEta );
      sigmaEtaPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaEtaPhi );
      sigmaPhiPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaPhiPhi );
      s1s9ClassEle        [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle       [iecal][iptbin][iclass] -> Fill ( s9s25 );
      LATClassEle         [iecal][iptbin][iclass] -> Fill ( lat );
      etaLATClassEle      [iecal][iptbin][iclass] -> Fill ( etaLat );
      phiLATClassEle      [iecal][iptbin][iclass] -> Fill ( phiLat );
      a20ClassEle         [iecal][iptbin][iclass] -> Fill ( a20 );
      a42ClassEle         [iecal][iptbin][iclass] -> Fill ( a42 );

      dPhiCaloFullclassEle    [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle     [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( dEta );
      EoPoutFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( HoE );
      shapeFisherFullclassEle [iecal][iptbin][ifullclass] -> Fill ( fisher );
      sigmaEtaEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaEtaEta );
      sigmaEtaPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaEtaPhi );
      sigmaPhiPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaPhiPhi );
      s1s9FullclassEle        [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle       [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
      LATFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( lat );
      etaLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( etaLat );
      phiLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( phiLat );
      a20FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a20 );
      a42FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a42 );

    } // loop on electrons

  } // loop on events
  
}

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
//     Utils anaUtils;
//     bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

//     if(!passedHLT) continue;

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
           ! m_selection->passCut("relSumPtTracks",eleSumPt04Ele[iele]) ) continue;
      
      
      /// define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = eleClassEle[iele];
        
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

      float sigmaEtaEta = sqrt(fabs(covEtaEtaEle[iele]));
      float sigmaEtaPhi = sqrt(fabs(covEtaPhiEle[iele]));
      float sigmaPhiPhi = sqrt(fabs(covPhiPhiEle[iele]));
      float s1s9 = s1s9Ele[iele];
      float s9s25 = s9s25Ele[iele];
      float lat = latEle[iele];
      float etaLat = etaLatEle[iele];
      float phiLat = phiLatEle[iele];
      float a20 = a20Ele[iele];
      float a42 = a42Ele[iele];
      float fisher = -1000;

      if ( iecal==0 ) { // barrel
        if ( iptbin == 0 ) {  // low pt
          fisher = 0.693496 - 12.7018 * sigmaEtaEta + 1.23863 * s9s25 - 10.115 * etaLat;
        }
        else if ( iptbin == 1 ) {  // high pt
          fisher = 6.02184 - 49.2656 * sigmaEtaEta + 2.49634 * s9s25 - 30.1528 * etaLat;
        }
      }
      else if ( iecal == 1  ) { // endcap
        if ( iptbin == 0 ) {  // low pt
          fisher = -1.11814 - 5.3288 * sigmaEtaEta + 4.51575 * s9s25 - 6.47578 * etaLat;
        }
        else if ( iptbin == 1 ) {  // high pt
          fisher = 0.536351 - 11.7401 * sigmaEtaEta + 3.61809 * s9s25 - 9.3025 * etaLat;
        }
      }

      double dPhiCalo = eleDeltaPhiAtCaloEle[iele];
      double dPhiVtx = eleDeltaPhiAtVtxEle[iele];
      double dEta = eleDeltaEtaAtVtxEle[iele];
      double EoPout = eleCorrEoPoutEle[iele];
      double HoE = eleHoEEle[iele];
      double dxy = eleTrackDxyEle[iele];
      double dxySig = eleTrackDxyEle[iele]/eleTrackDxyErrorEle[iele];

      dPhiCaloUnsplitEle    [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle     [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle        [iecal][iptbin] -> Fill ( dEta );
      EoPoutUnsplitEle      [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle         [iecal][iptbin] -> Fill ( HoE );
      shapeFisherUnsplitEle [iecal][iptbin] -> Fill ( fisher );
      sigmaEtaEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaEtaEta );
      sigmaEtaPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaEtaPhi );
      sigmaPhiPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaPhiPhi );
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
      shapeFisherClassEle [iecal][iptbin][iclass] -> Fill ( fisher );
      sigmaEtaEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaEtaEta );
      sigmaEtaPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaEtaPhi );
      sigmaPhiPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaPhiPhi );
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
      shapeFisherFullclassEle [iecal][iptbin][ifullclass] -> Fill ( fisher );
      sigmaEtaEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaEtaEta );
      sigmaEtaPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaEtaPhi );
      sigmaPhiPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaPhiPhi );
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
  float fisherMin   = -15.0;
  float fisherMax   =  15.0;
  float sigmaEtaEtaMin = 0.0;
  float sigmaEtaEtaMax = 0.05;
  float sigmaEtaPhiMin = 0.0;
  float sigmaEtaPhiMax = 0.05;
  float sigmaPhiPhiMin = 0.0;
  float sigmaPhiPhiMax = 0.05;
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
      sprintf(histo,"shapeFisherUnsplit_electrons_%d_%d",iecal,iptbin);
      shapeFisherUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
      sprintf(histo,"sigmaEtaEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaEtaEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
      sprintf(histo,"sigmaEtaPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaEtaPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
      sprintf(histo,"sigmaPhiPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaPhiPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
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
	sprintf(histo,"shapeFisherClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	shapeFisherClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
	sprintf(histo,"sigmaEtaEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaEtaEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
	sprintf(histo,"sigmaEtaPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaEtaPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
	sprintf(histo,"sigmaPhiPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaPhiPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
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
	sprintf(histo,"shapeFisherFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	shapeFisherFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
	sprintf(histo,"sigmaEtaEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaEtaEtaFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
	sprintf(histo,"sigmaEtaPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaEtaPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
	sprintf(histo,"sigmaPhiPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaPhiPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
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

  TFile *file = TFile::Open(filename,"recreate");
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
      shapeFisherUnsplitEle[iecal][iptbin]->Write();
      sigmaEtaEtaUnsplitEle[iecal][iptbin]->Write();
      sigmaEtaPhiUnsplitEle[iecal][iptbin]->Write();
      sigmaPhiPhiUnsplitEle[iecal][iptbin]->Write();
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
	shapeFisherClassEle[iecal][iptbin][iclass]->Write();
	sigmaEtaEtaClassEle[iecal][iptbin][iclass]->Write();
	sigmaEtaPhiClassEle[iecal][iptbin][iclass]->Write();
	sigmaPhiPhiClassEle[iecal][iptbin][iclass]->Write();
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
	shapeFisherFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaEtaEtaFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaEtaPhiFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaPhiPhiFullclassEle[iecal][iptbin][ifullclass]->Write();
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
