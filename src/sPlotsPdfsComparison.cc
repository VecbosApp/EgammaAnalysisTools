#define sPlotsPdfsComparison_cxx
#include "include/sPlotsPdfsComparison.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

using namespace std;

void sPlotsPdfsComparison::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L sPlotsPdfsComparison.C
//      Root > sPlotsPdfsComparison t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   bookHistos();

   Long64_t nentries = fChain->GetEntries();

   int firstEvent = 0;
   int lastEvent = 0;
   if( m_isMC ) {
     lastEvent = nentries/2 - 1;
   } else {
     firstEvent = nentries/2;
     lastEvent = nentries;
   }

   cout << firstEvent << " " << lastEvent << endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=firstEvent; jentry<lastEvent;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      int jecal = (int)iecal;
      int jptbin = (int)iptbin;
      int jclass = (int)iclass;

      float wgt = (m_isMC) ? 1.0 : N_sig_sw;

      dPhiClassEle          [jecal][jptbin][jclass] -> Fill ( deltaPhi, wgt );
      dEtaClassEle          [jecal][jptbin][jclass] -> Fill ( deltaEta, wgt );
      EoPoutClassEle        [jecal][jptbin][jclass] -> Fill ( EoPout, wgt );
      EoPClassEle           [jecal][jptbin][jclass] -> Fill ( EoP, wgt );
      HoEClassEle           [jecal][jptbin][jclass] -> Fill ( HoE, wgt );
      sigmaIEtaIEtaClassEle [jecal][jptbin][jclass] -> Fill ( sigmaIEtaIEta, wgt );
      s1s9ClassEle          [jecal][jptbin][jclass] -> Fill ( s1s9, wgt );
      s9s25ClassEle         [jecal][jptbin][jclass] -> Fill ( s9s25, wgt );

   }

   char buffer[200];
   if( m_isMC ) sprintf(buffer,"pdfs_histograms_MC.root");
   else sprintf(buffer,"pdfs_histograms_data.root");
   TFile *fileOut = TFile::Open(buffer,"recreate");

   for (int jecal=0; jecal<2; jecal++) {
     for(int jptbin=0; jptbin<2; jptbin++) {
       for(int jclass=0; jclass<2; jclass++) {
        dPhiClassEle[jecal][jptbin][jclass]->Write();
        dEtaClassEle[jecal][jptbin][jclass]->Write();
        EoPoutClassEle[jecal][jptbin][jclass]->Write();
        EoPClassEle[jecal][jptbin][jclass]->Write();           
        HoEClassEle[jecal][jptbin][jclass]->Write();
        sigmaIEtaIEtaClassEle[jecal][jptbin][jclass]->Write();
        s1s9ClassEle[jecal][jptbin][jclass]->Write();
        s9s25ClassEle[jecal][jptbin][jclass]->Write();
       }
     }
   }
   
   fileOut->Close();


}

void sPlotsPdfsComparison::bookHistos() {

  int nbins = 100;

  float dPhiMin = -0.15;
  float dPhiMax =  0.15;
  float dEtaMin     = -0.02;
  float dEtaMax     =  0.02;
  float EoPoutMin   =  0.0;
  float EoPoutMax   =  20.0;
  float EoPMin   =  0.0;
  float EoPMax   =  20.0;
  float HoEMin      =  0.0; // zero-suppression in HCAL
  float HoEMax      =  0.1; // ??
  float sigmaIEtaIEtaMin = 0.0;
  float sigmaIEtaIEtaMax = 0.05;
  float s1s9Min  = 0.0;
  float s1s9Max  = 1.0;
  float s9s25Min = 0.5;
  float s9s25Max = 1.0;

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {

    // iptbin = 0: < 15 GeV
    // iptbin = 1: > 15 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      // iclass = 0: non-showering
      // iclass = 1: showering
      for(int iclass=0; iclass<2; iclass++) {

        char histo[200];
      
	sprintf(histo,"dPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiClassEle[iecal][iptbin][iclass]    = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
	sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	EoPoutClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"EoPClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	EoPClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
	sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
	sprintf(histo,"s1s9Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s1s9ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s9s25ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);

      }

    }

  }

}
