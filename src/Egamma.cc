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
            if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
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
  measurements.pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
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
