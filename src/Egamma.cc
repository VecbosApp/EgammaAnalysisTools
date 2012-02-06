#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <fstream>
#include <sstream>

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "EgammaAnalysisTools/include/Egamma.h"

using namespace bits;

Egamma::Egamma(TTree *tree) : EgammaBase(tree)
{
  jsonFile = "";
  lastFile = "";
}

Egamma::~Egamma()
{
  // By this time, the destructor of EgammaBase has not yet been called.
  // This means that the tree has not yet been deleted.
  // So, we do nothing here.
}

void Egamma::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggers=reqTriggers;
}

bool Egamma::hasPassedHLT() {
  Utils anaUtils;
  return anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
}

void Egamma::setJsonGoodRunList(const std::string& jsonFilePath)
{
  jsonFile=jsonFilePath;
}

void Egamma::fillRunLSMap()
{
  
  if (jsonFile == "")
    {
      std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
      return;
    }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open())
    {
      std::cout << "Unable to open file " << jsonFile << std::endl;
      return;
    }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
    {
      const json::Array& lsSegment = (*itRun).element;
      LSSegments thisRunSegments; 
      for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
	{
	  json::Array lsSegment=(*lsIterator);
	  json::Number lsStart=lsSegment[0];	   
	  json::Number lsEnd=lsSegment[1];
	  aLSSegment thisSegment;
	  thisSegment.first=lsStart.Value();
	  thisSegment.second=lsEnd.Value();
	  thisRunSegments.push_back(thisSegment);
	  //	   std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]); 
	}
      goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
    }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR)
    {
      std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
      for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
      std::cout << std::endl;
    }
}

bool Egamma::isGoodRunLS()
{
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //  std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //      std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
	return true;
    }
  return false;
}

bool Egamma::reloadTriggerMask(bool newVersion)
{
  if(newVersion) {
    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
      {   
        for(unsigned int i=0; i<nameHLT->size(); i++)
          {
            //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
            // nameHLT[i] has ..._vXXX
            if(nameHLT->at(i).find(*fIter) != std::string::npos)
              {
                triggerMask.push_back( indexHLT[i] ) ;
                break;
              }
          }
      }
    m_requiredTriggers = triggerMask;
  } else {
    TString fileName=((TChain*)fChain)->GetFile()->GetName();
    if ( TString(lastFile) != fileName )
      {

        std::cout << "[ReloadTriggerMask]::File has changed reloading trigger mask" << std::endl;
        lastFile = fileName;
        TTree *treeCond;
        std::cout << "[ReloadTriggerMask]::Opening " << fileName << std::endl;
        treeCond = (TTree*)((TChain*)fChain)->GetFile()->Get("Conditions");
        int           nHLT_;
        std::vector<std::string>  *nameHLT_;
        std::vector<unsigned int> *indexHLT_;

        //To get the pointers for the vectors
        nameHLT_=0;
        indexHLT_=0;

        treeCond->SetBranchAddress("nHLT", &nHLT_);
        treeCond->SetBranchAddress("nameHLT", &nameHLT_);
        treeCond->SetBranchAddress("indexHLT", &indexHLT_);
        treeCond->GetEntry(0);

        std::vector<int> triggerMask;
        for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
          {
            for(unsigned int i=0; i<nameHLT_->size(); i++) 
              {
                if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
                  {
                    triggerMask.push_back( indexHLT_->at(i) ) ;
                    break;
                  }
              }
          }
        m_requiredTriggers = triggerMask;
        for (int i=0;i<m_requiredTriggers.size();++i)
          std::cout << "[ReloadTriggerMask]::Requiring bit " << m_requiredTriggers[i] << " " << requiredTriggers[i] << std::endl;
      }
  }
}

float Egamma::mT3(TLorentzVector pl1, TLorentzVector pl2, TVector3 met) {
  float pTll = (pl1.Vect() + pl2.Vect()).Pt();
  float mll = (pl1 + pl2).M();
  float El = sqrt(pTll*pTll + mll*mll);
  float pTnu = met.Pt();
  float Enu = sqrt(pTnu*pTnu + mll*mll);
  float Ex = (pl1+pl2).X() + met.X();
  float Ey = (pl1+pl2).Y() + met.Y();
  float mnu = mll;

  return sqrt(mll*mll + mnu*mnu + 2*(El*Enu-Ex*Ex-Ey*Ey));
}

float Egamma::likelihoodRatio(int eleIndex, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  Utils anaUtils;
  bool inEB=anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB);
  measurements.pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
  if(inEB && fabs(etaEle[eleIndex])<1.0) measurements.subdet = 0;
  else if (inEB && fabs(etaEle[eleIndex])>=1.0) measurements.subdet = 1;
  else measurements.subdet = 2;
  measurements.deltaPhi = deltaPhiAtVtxEle[eleIndex];
  measurements.deltaEta = deltaEtaAtVtxEle[eleIndex];
  measurements.eSeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  measurements.eSuperClusterOverP = eSuperClusterOverPEle[eleIndex];
  measurements.hadronicOverEm = hOverEEle[eleIndex];
  measurements.sigmaIEtaIEta = SigmaiEiE(eleIndex);
  measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
  measurements.fBrem = fbremEle[eleIndex];
  measurements.nBremClusters = nbremsEle[eleIndex];
  int gsftrack = gsfTrackIndexEle[eleIndex];
  TVector3 pIn(pxGsfTrack[gsftrack],pyGsfTrack[gsftrack],pzGsfTrack[gsftrack]);
  measurements.OneOverEMinusOneOverP = 1./(eSuperClusterOverPEle[eleIndex]*pIn.Mag()) - 1./pIn.Mag();
  return lh.resultLog(measurements);
}

/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float Egamma::SigmaiEiE(int electron) {
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
float Egamma::SigmaiPiP(int electron) {
  float spp;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    spp = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      spp = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      spp = 999.;
    }
  }
  return spp;
}

bool Egamma::triggerMatch(float eta, float phi, float Dr){

  bool match=false;
  for( int i=0; i<m_requiredTriggers.size(); i++ ) {  // loop over require trigger paths
    
    int pathIndex=m_requiredTriggers[i];
    // std::cout << "testing trigger " << pathIndex << " with " << sizePassing[pathIndex] << " passing objects" << std::endl; 
    
    if( sizePassing[pathIndex]>  0 ) {  //some object has passed the required trigger 
      
      for(int np = 0; np < sizePassing[pathIndex]; np++ ){
        int iP = indexPassing[ indexPassingPerPath[pathIndex] +np];
        // std::cout << "passing object eta: " << triggerObsEta[iP] << " phi: " <<  triggerObsPhi[iP] << std::endl; 

        if(DeltaR(eta, phi,triggerObsEta[iP],  triggerObsPhi[iP] ) < Dr){
          match=true;
          //std::cout << "MATCH!" <<std::endl;	
          break;
        }
      }
    }
    if(match)  //it's enough if one path matches	
      break;
  }
  return match;
}

// dxy parameter with respect to PV for tracks
double Egamma::trackDxyPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  return ( - (trackVPos.X()-PVPos.X())*trackMom.Y() + (trackVPos.Y()-PVPos.Y())*trackMom.X() ) / trackMom.Pt(); 
}

/// dz parameter with respect to PV for tracks
double Egamma::trackDzPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  float trackPt = trackMom.Pt();
  return (trackVPos.Z()-PVPos.Z()) - ((trackVPos.X()-PVPos.X())*trackMom.X()+(trackVPos.Y()-PVPos.Y())*trackMom.Y())/trackPt *trackMom.Pz()/trackPt; 
}

/// dsz parameter with respect to PV for tracks
double Egamma::trackDszPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  float trackPt = trackMom.Pt();
  float trackP  = trackMom.Mag();
  return (trackVPos.Z()-PVPos.Z())*trackPt/trackP - ((trackVPos.X()-PVPos.X())*trackMom.X()+(trackVPos.Y()-PVPos.Y())*trackMom.Y())/trackPt *trackMom.Pz()/trackP; 
}

/// dxy, dz and dsz parameters with respect to PV for electrons
double Egamma::eleDxyPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDxyPV(PVPos,lepVPos,lepMom);
}

double Egamma::eleDzPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDzPV(PVPos,lepVPos,lepMom);
}

double Egamma::eleDszPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDszPV(PVPos,lepVPos,lepMom);
}

// using HWW BDT (paper 2011) 
float Egamma::eleBDT(ElectronIDMVA *mva, int eleIndex) {
  
  if(mva==0) {
    std::cout << "electron BDT not created/initialized. BIG PROBLEM. Returning false output -999!" << std::endl; 
    return -999.;
  }
  
  int gsfTrack = gsfTrackIndexEle[eleIndex]; 
  double gsfsign   = (-eleDxyPV(eleIndex,0) >=0 ) ? 1. : -1.;

  float ElePt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  float EleDEtaIn = deltaEtaAtVtxEle[eleIndex];
  float EleDPhiIn = deltaPhiAtVtxEle[eleIndex];
  float EleHoverE = hOverEEle[eleIndex];
  float EleD0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
  float EleFBrem = fbremEle[eleIndex];
  float EleEOverP = eSuperClusterOverPEle[eleIndex];
  float EleESeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  float EleNBrem = nbremsEle[eleIndex];
  TVector3 pInGsf(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);

  float EleIP3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
  float EleIP3dSig = EleIP3d/impactPar3DErrorGsfTrack[gsfTrack];

  // we have not pout and seed cluster energy in the trees. Gymnastyc...
  float Pout = pInGsf.Mag() - fbremEle[eleIndex] * pInGsf.Mag();
  float SCSeedEnergy = EleESeedClusterOverPout * Pout;
  float EleESeedClusterOverPIn = SCSeedEnergy/pInGsf.Mag();      

  float EleSigmaIEtaIEta, EleSigmaIPhiIPhi, EleOneOverEMinusOneOverP, EleSCEta;

  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    EleSigmaIEtaIEta = sqrt(covIEtaIEtaSC[sc]);
    EleSigmaIPhiIPhi = sqrt(covIPhiIPhiSC[sc]);
    EleOneOverEMinusOneOverP = 1./energySC[sc]  - 1./pInGsf.Mag();
    EleSCEta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      EleSigmaIEtaIEta = sqrt(covIEtaIEtaPFSC[sc]);
      EleSigmaIPhiIPhi = sqrt(covIPhiIPhiPFSC[sc]);
      EleOneOverEMinusOneOverP = 1./energyPFSC[sc]  - 1./pInGsf.Mag();
      EleSCEta = etaPFSC[sc];
    } else {
      EleSigmaIEtaIEta = 999.;
      EleSigmaIPhiIPhi = 999.;
      EleOneOverEMinusOneOverP = 999.;
      EleESeedClusterOverPIn = 999.;
      EleSCEta = 0.;
    }
  }

  return mva->MVAValue(ElePt, EleSCEta,
                       EleSigmaIEtaIEta,
                       EleDEtaIn,
                       EleDPhiIn,
                       EleHoverE,
                       EleD0,
                       EleFBrem,
                       EleEOverP,
                       EleESeedClusterOverPout,
                       EleSigmaIPhiIPhi,
                       EleNBrem,
                       EleOneOverEMinusOneOverP,
                       EleESeedClusterOverPIn,
                       EleIP3d,
                       EleIP3dSig );

}

// using HZZ BDT
float Egamma::eleBDT(ElectronIDMVAHZZ *mva, int eleIndex) {
  
  if(mva==0) {
    std::cout << "electron BDT not created/initialized. BIG PROBLEM. Returning false output -999!" << std::endl; 
    return -999.;
  }
  
  int gsfTrack = gsfTrackIndexEle[eleIndex]; 
  int kfTrack = trackIndexEle[eleIndex];
  float ElePt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  float EleFBrem = fbremEle[eleIndex];
  float EleDEtaIn = deltaEtaAtVtxEle[eleIndex];
  float EleDPhiIn = deltaPhiAtVtxEle[eleIndex];
  float EleHoverE = hOverEEle[eleIndex];
  float EleSuperClusterEOverP = eSuperClusterOverPEle[eleIndex];
  float EleEOverPout = eEleClusterOverPoutEle[eleIndex];
  float EleDEtaEleOut = deltaEtaEleClusterTrackAtCaloEle[eleIndex];
  float EleKFChi2 = (kfTrack>-1) ? trackNormalizedChi2Track[kfTrack] : 0.0;
  float EleKFHits = (kfTrack>-1) ? trackerLayersWithMeasurementTrack[kfTrack] : -1.0;
  float EleMissHits = expInnerLayersGsfTrack[gsfTrack];
  float EleDistConv = convDistEle[eleIndex];
  float EleDcotConv = convDcotEle[eleIndex];
  float NVtx = float(nPV);
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  float EleEcalSeeded = (ecaldriven) ? 1. : 0.;

  float EleSigmaIEtaIEta, EleE1x5E5x5, EleSCEta;

  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    EleSigmaIEtaIEta = sqrt(covIEtaIEtaSC[sc]);
    EleE1x5E5x5 = (e5x5SC[sc] - e1x5SC[sc])/e5x5SC[sc];
    EleSCEta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      EleSigmaIEtaIEta = sqrt(covIEtaIEtaPFSC[sc]);
      EleE1x5E5x5 = (e5x5PFSC[sc] - e1x5PFSC[sc])/e5x5PFSC[sc];
      EleSCEta = etaPFSC[sc];
    } else {
      EleSigmaIEtaIEta = 999.;
      EleE1x5E5x5 = 999.;
      EleSCEta = 0.;
    }
  }

  return mva->MVAValue(ElePt, EleSCEta,
                       EleFBrem,
                       EleDEtaIn,
                       EleDPhiIn,
                       EleDEtaEleOut,
                       EleSigmaIEtaIEta,
                       EleHoverE,
                       EleSuperClusterEOverP,
                       EleE1x5E5x5,
                       EleEOverPout,
                       EleKFChi2,
                       EleKFHits,
                       EleMissHits,
                       EleDistConv,
                       EleDcotConv,
                       NVtx,
                       EleEcalSeeded);

}

bool Egamma::passEleBDT(float pt, float eta, float bdt) {

  if(pt < 20 && fabs(eta) < 1.0) return (bdt > 0.139);
  if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdt > 0.525);
  if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdt > 0.543);
  if(pt >= 20 && fabs(eta) < 1.0) return (bdt > 0.947);
  if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdt > 0.950);
  if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdt > 0.884);

  // here we are cutting the events with |SC eta|>2.5. If the acceptance is done with |ele eta|<2.5 then will cut some event more. Fine. Synch with this.
  return false;

}
