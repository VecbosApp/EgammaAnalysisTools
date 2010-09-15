#define prepareQCDhistos_cxx
#include "prepareQCDhistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void prepareQCDhistos::Loop() {

  if (fChain == 0) return;
  
  // histos booking
  bookFullHistos();


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
    
    // filling histos
    dPhiUnsplitEle  [jecal][jptbin] -> Fill ( dphi );
    dEtaUnsplitEle  [jecal][jptbin] -> Fill ( deta );
    EoPUnsplitEle   [jecal][jptbin] -> Fill ( EoP );
    HoEUnsplitEle   [jecal][jptbin] -> Fill ( HoE );
    fBremUnsplitEle [jecal][jptbin] -> Fill ( fbrem );
    sigmaIEtaIEtaUnsplitEle [jecal][jptbin] -> Fill ( see );
    sigmaIPhiIPhiUnsplitEle [jecal][jptbin] -> Fill ( spp );
    
    dPhiClassEle  [jecal][jptbin][jclass] -> Fill ( dphi );
    dEtaClassEle  [jecal][jptbin][jclass] -> Fill ( deta );
    EoPClassEle   [jecal][jptbin][jclass] -> Fill ( EoP );
    HoEClassEle   [jecal][jptbin][jclass] -> Fill ( HoE );
    fBremClassEle [jecal][jptbin][jclass] -> Fill ( fbrem );
    sigmaIEtaIEtaClassEle [jecal][jptbin][jclass] -> Fill ( see );
    sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( spp );
    
  } // end loop over entries


  // saving likelihood pdfs
  char buffer[200];
  char bufferLikelihood[200];
  sprintf(buffer,"pdfsQCD_histograms_data.root");
  sprintf(bufferLikelihood,"pdfsQCD_data.root");
  
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

      for(int iclass=0; iclass<2; iclass++) {
      	dPhiClassEle[iecal][iptbin][iclass] -> Write();
	dEtaClassEle[iecal][iptbin][iclass] -> Write();
	EoPClassEle[iecal][iptbin][iclass]  -> Write();
	HoEClassEle[iecal][iptbin][iclass]  -> Write();
        fBremClassEle[iecal][iptbin][iclass]-> Write();
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] -> Write();
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] -> Write();
      }
    }
  }
  
  fileOut2->Close();
}

void prepareQCDhistos::bookFullHistos() {

  int nbins = 50;
  
  float dPhiMin  = -0.05; // pixelMatchGsfElectron pre-selection: |dPhi| < 0.1
  float dPhiMax  =  0.05;
  float dEtaMin  = -0.008;
  float dEtaMax  =  0.008;
  float EoPMin   =  0.0;
  float EoPMax   =  5.0;
  float HoEMin   =  0.0;   // zero-suppression in HCAL
  float HoEMax   =  0.15;  // pixelMatchGsfElectron pre-selection: H/E<0.15
  float sigmaIEtaIEtaEBMin = 0.0;
  float sigmaIEtaIEtaEBMax = 0.03;
  float sigmaIEtaIEtaEEMin = 0.01;
  float sigmaIEtaIEtaEEMax = 0.05;
  float sigmaIPhiIPhiEBMin = 0.0;
  float sigmaIPhiIPhiEBMax = 0.03;
  float sigmaIPhiIPhiEEMin = 0.01;
  float sigmaIPhiIPhiEEMax = 0.05;
  float fBremMin = 0.0;
  float fBremMax = 1.0;

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {

    // iptbin = 0: 10 pt < 20 GeV    
    // iptbin = 1: pt > 20 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      char histo[200];
      
      sprintf(histo,"dPhiUnsplit_hadrons_%d_%d",iecal,iptbin);
      dPhiUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
      sprintf(histo,"dEtaUnsplit_hadrons_%d_%d",iecal,iptbin);
      dEtaUnsplitEle[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      sprintf(histo,"EoPUnsplit_hadrons_%d_%d",iecal,iptbin);
      EoPUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
      sprintf(histo,"HoEUnsplit_hadrons_%d_%d",iecal,iptbin);
      HoEUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      if(iecal==0) {
        sprintf(histo,"sigmaIEtaIEtaUnsplit_hadrons_%d_%d",iecal,iptbin);
        sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEBMin, sigmaIEtaIEtaEBMax);
        sprintf(histo,"sigmaIPhiIPhiUnsplit_hadrons_%d_%d",iecal,iptbin);
        sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEBMin, sigmaIPhiIPhiEBMax);
      } else {
        sprintf(histo,"sigmaIEtaIEtaUnsplit_hadrons_%d_%d",iecal,iptbin);
        sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEEMin, sigmaIEtaIEtaEEMax);
        sprintf(histo,"sigmaIPhiIPhiUnsplit_hadrons_%d_%d",iecal,iptbin);
        sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEEMin, sigmaIPhiIPhiEEMax);
      }
      sprintf(histo,"fBremUnsplit_hadrons_%d_%d",iecal,iptbin);
      fBremUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, fBremMin, fBremMax);

      // iclass = 0: 0 - brem clusters
      // iclass = 1: >=1 - brem clusters
      for(int iclass=0; iclass<2; iclass++) {
	sprintf(histo,"dPhiClass_hadrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
	sprintf(histo,"dEtaClass_hadrons_%d_%d_%d",iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPClass_hadrons_%d_%d_%d",iecal,iptbin,iclass);
	EoPClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
	sprintf(histo,"HoEClass_hadrons_%d_%d_%d",iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
        if(iecal==0) {
          sprintf(histo,"sigmaIEtaIEtaClass_hadrons%d_%d_%d",iecal,iptbin,iclass);
          sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEBMin, sigmaIEtaIEtaEBMax);
          sprintf(histo,"sigmaIPhiIPhiClass_hadrons%d_%d_%d",iecal,iptbin,iclass);
          sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEBMin, sigmaIPhiIPhiEBMax);
        } else {
          sprintf(histo,"sigmaIEtaIEtaClass_hadrons%d_%d_%d",iecal,iptbin,iclass);
          sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEEMin, sigmaIEtaIEtaEEMax);
          sprintf(histo,"sigmaIPhiIPhiClass_hadrons%d_%d_%d",iecal,iptbin,iclass);
          sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEEMin, sigmaIPhiIPhiEEMax);
        }
        sprintf(histo,"fBremClass_hadrons_%d_%d_%d",iecal,iptbin,iclass);
        fBremClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, fBremMin, fBremMax);
      }
    }
  }
}
