#define prepareQCDhistos_cxx
#include "prepareQCDhistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <iostream>
#include <math.h>

using namespace std;

void prepareQCDhistos::Loop() {

  if (fChain == 0) return;
  
  // histos booking
  bookFullHistos();

  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *LHdir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  LH = new ElectronLikelihood(&(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);

  Long64_t nentries = fChain->GetEntriesFast();
  int firstEvent = 0;
  int lastEvent  = nentries;
  std::cout << firstEvent << " " << lastEvent << std::endl;

  Long64_t nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    if (jentry%10000 == 0) std::cout << "Processing event # " << jentry << " / " << nentries << std::endl;

    // remove residual spikes (if the cleaning is not done)    
    if(see<1e-4) continue; 
    
    // eta, pt and class 
    int jecal, jptbin, jclass;
    if (fabs(eta)<1.479) jecal = 0;
    else if (fabs(eta)>=1.479) jecal = 1;
    else continue;

    if(pt>=10 && pt<=20) jptbin = 0;
    else if(pt>20) jptbin = 1;
    else continue;

    if(nbrem==0) jclass = 0;
    else jclass = 1;
    
    // filling histos: unsplit
    dPhiUnsplitEle  [jecal][jptbin] -> Fill ( dphi );
    dEtaUnsplitEle  [jecal][jptbin] -> Fill ( deta );
    EoPUnsplitEle   [jecal][jptbin] -> Fill ( EoP );
    HoEUnsplitEle   [jecal][jptbin] -> Fill ( HoE );
    sigmaIEtaIEtaUnsplitEle [jecal][jptbin] -> Fill ( see );
    fBremUnsplitEle [jecal][jptbin] -> Fill ( fbrem );
    sigmaIPhiIPhiUnsplitEle [jecal][jptbin] -> Fill ( spp );
    
    // splitted
    dPhiClassEle  [jecal][jptbin][jclass] -> Fill ( dphi );
    dEtaClassEle  [jecal][jptbin][jclass] -> Fill ( deta );
    EoPClassEle   [jecal][jptbin][jclass] -> Fill ( EoP );
    HoEClassEle   [jecal][jptbin][jclass] -> Fill ( HoE );
    sigmaIEtaIEtaClassEle [jecal][jptbin][jclass] -> Fill ( see );

    // where not used, flat distribution
    float minSPP   = sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinLowEdge(0)-1.0;
    float maxSPP   = sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinUpEdge(sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetNbinsX()+1)+1.0;
    float minFBrem = fBremClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinLowEdge(0)-1.0;
    float maxFBrem = fBremClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinUpEdge(fBremClassEle [jecal][jptbin][jclass]->GetNbinsX()+1)+1.0;
    TRandom3 rnd(jentry);
    // to create PDFs 
    if (jclass==0) fBremClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minFBrem,maxFBrem) );
    // for performance plots
    // if (jclass==0) fBremClassEle [jecal][jptbin][jclass] -> Fill ( fbrem );
    if (jclass==1) fBremClassEle [jecal][jptbin][jclass] -> Fill ( fbrem );    
    if( jclass==0) sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( spp );
    // to create PDFs     
    if( jclass==1) sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minSPP,maxSPP) );
    // for performance plots 
    // if( jclass==1) sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( spp );

    // the likelihood output
    float f_lh = likelihoodRatio(*LH);
    lhUnsplitEle [jecal][jptbin] -> Fill ( f_lh );
    lhClassEle [jecal][jptbin][jclass] -> Fill ( f_lh );
    
  } // end loop over entries


  // saving likelihood pdfs
  char buffer[200];
  char bufferLikelihood[200];
  sprintf(buffer,"pdfsQCD_histograms.root");
  sprintf(bufferLikelihood,"pdfsQCD.root");
  
  TFile *fileOut2 = TFile::Open(bufferLikelihood, "recreate");
  
  for (int iecal=0; iecal<2; iecal++) {
    
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      dPhiUnsplitEle[iecal][iptbin] -> Write();
      dEtaUnsplitEle[iecal][iptbin] -> Write();
      EoPUnsplitEle[iecal][iptbin]  -> Write();
      HoEUnsplitEle[iecal][iptbin]  -> Write();
      fBremUnsplitEle[iecal][iptbin]-> Write();
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin] -> Write();
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin] -> Write();
      lhUnsplitEle[iecal][iptbin]->Write();

      for(int iclass=0; iclass<2; iclass++) {
      	dPhiClassEle[iecal][iptbin][iclass] -> Write();
	dEtaClassEle[iecal][iptbin][iclass] -> Write();
	EoPClassEle[iecal][iptbin][iclass]  -> Write();
	HoEClassEle[iecal][iptbin][iclass]  -> Write();
        fBremClassEle[iecal][iptbin][iclass]-> Write();
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] -> Write();
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] -> Write();
	lhClassEle[iecal][iptbin][iclass]->Write();
      }
    }
  }
  
  fileOut2->Close();
}

void prepareQCDhistos::bookFullHistos() {

  int nbins = 200;
  
  float dPhiMin  = -0.1; // pixelMatchGsfElectron pre-selection: |dPhi| < 0.1
  float dPhiMax  =  0.1;
  float dEtaMin  = -0.015;
  float dEtaMax  =  0.015;
  float EoPMin   =  0.0;
  float EoPMax   =  10.0;
  float HoEMin   =  0.0;   // zero-suppression in HCAL
  float HoEMax   =  0.15;  // pixelMatchGsfElectron pre-selection: H/E<0.15
  float sigmaIEtaIEtaEBMin = 0.0;
  float sigmaIEtaIEtaEBMax = 0.03;
  float sigmaIEtaIEtaEEMin = 0.01;
  float sigmaIEtaIEtaEEMax = 0.05;
  float sigmaIPhiIPhiEBMin = 0.0;
  float sigmaIPhiIPhiEBMax = 0.03;
  float sigmaIPhiIPhiEEMin = 0.01;
  float sigmaIPhiIPhiEEMax = 0.09;
  float fBremMin = 0.0;
  float fBremMax = 1.0;
  float lhMin = 0.0;
  float lhMax = 1.0;

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {

    // iptbin = 0: 10 pt < 20 GeV    
    // iptbin = 1: pt > 20 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      char histo[200];
      
      sprintf(histo,"dPhiUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
      dPhiUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
      sprintf(histo,"dEtaUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
      dEtaUnsplitEle[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      sprintf(histo,"EoPUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
      EoPUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
      sprintf(histo,"HoEUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
      HoEUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      if(iecal==0) {
        sprintf(histo,"sigmaIEtaIEtaUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
        sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEBMin, sigmaIEtaIEtaEBMax);
        sprintf(histo,"sigmaIPhiIPhiUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
        sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEBMin, sigmaIPhiIPhiEBMax);
      } else {
        sprintf(histo,"sigmaIEtaIEtaUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
        sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEEMin, sigmaIEtaIEtaEEMax);
        sprintf(histo,"sigmaIPhiIPhiUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
        sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEEMin, sigmaIPhiIPhiEEMax);
      }
      sprintf(histo,"fBremUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
      fBremUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, fBremMin, fBremMax);

      sprintf(histo,"lhUnsplit_hadrons_subdet%d_ptbin%d",iecal,iptbin);
      lhUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, lhMin, lhMax);

      // iclass = 0: 0 - brem clusters
      // iclass = 1: >=1 - brem clusters
      for(int iclass=0; iclass<2; iclass++) {
	sprintf(histo,"dPhiClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
	dPhiClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
	sprintf(histo,"dEtaClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
	EoPClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
	sprintf(histo,"HoEClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
        if(iecal==0) {
          sprintf(histo,"sigmaIEtaIEtaClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
          sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEBMin, sigmaIEtaIEtaEBMax);
          sprintf(histo,"sigmaIPhiIPhiClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
          sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEBMin, sigmaIPhiIPhiEBMax);
        } else {
          sprintf(histo,"sigmaIEtaIEtaClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
          sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEEMin, sigmaIEtaIEtaEEMax);
          sprintf(histo,"sigmaIPhiIPhiClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
          sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEEMin, sigmaIPhiIPhiEEMax);
        }
        sprintf(histo,"fBremClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
        fBremClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, fBremMin, fBremMax);

	sprintf(histo,"lhClass_hadrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
	lhClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, lhMin, lhMax);
      }
    }
  }
}

float prepareQCDhistos::likelihoodRatio(ElectronLikelihood &lh) {

  LikelihoodMeasurements measurements;

  measurements.pt = pt;
  measurements.subdet = (fabs(eta)<1.479) ? 0 : 1;
  measurements.deltaPhi = dphi;
  measurements.deltaEta = deta;
  measurements.eSeedClusterOverPout = -1;
  measurements.eSuperClusterOverP = EoP;
  measurements.hadronicOverEm = HoE;
  measurements.sigmaIEtaIEta = see;
  measurements.sigmaIPhiIPhi = spp;
  measurements.fBrem = fbrem;
  measurements.nBremClusters = nbrem;
  return lh.result(measurements);
}
