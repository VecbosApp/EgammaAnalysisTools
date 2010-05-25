#define sPlotsPdfsComparison_cxx
#include "include/sPlotsPdfsComparison.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <math.h>

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
  
  if (fChain == 0) return;

  bookHistosFixedBinning();
  InitCuts();
  float sel, norm;
  sel = norm = 0;
  
  Long64_t nentries = fChain->GetEntries();
  
  int firstEvent = 0;
  int lastEvent  = nentries;
  cout << firstEvent << " " << lastEvent << endl;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=firstEvent; jentry<lastEvent;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    int jecal;
    if(m_isMC) {
      if (fabs(f_eta)<1.46)  jecal = 0;
      else if (fabs(f_eta)>=1.46) jecal = 1;
      else continue;
    } else {
      if (fabs(eta)<1.46)  jecal = 0;
      else if (fabs(eta)>=1.46) jecal = 1;
      else continue;
    }

    // da cambiare (solo x dati, selezionare sgn o fondo)
    float wgt = (m_isMC) ? 1.0 : N_sig_sw;
    // float wgt = (m_isMC) ? 1.0 : N_qcd_sw;

    vector<float> WP_inf, WP_sup;
    if(jecal == 0) {
      WP_inf = WP70_EB_inf;
      WP_sup = WP70_EB_sup;
    } else {
      WP_inf = WP70_EE_inf;
      WP_sup = WP70_EE_sup;
    }

    if(m_isMC) {
      dPhiClassEle          [jecal] -> Fill ( f_dphi, wgt );
      dEtaClassEle          [jecal] -> Fill ( f_deta, wgt );
      EoPClassEle           [jecal] -> Fill ( f_eop, wgt );
      HoEClassEle           [jecal] -> Fill ( f_hoe, wgt );
      sigmaIEtaIEtaClassEle [jecal] -> Fill ( f_see, wgt );
      etaClassEle                   -> Fill ( f_eta, wgt );
      
      norm += wgt;
      if(f_see >= WP_inf[0] && f_see <= WP_sup[0])
        //         f_dphi >= WP_inf[1] && f_dphi <= WP_sup[1] &&
        //         f_deta >= WP_inf[2] && f_deta <= WP_sup[2] &&
        //         f_hoe >= WP_inf[3] && f_hoe <= WP_sup[3]) 
         sel += wgt;         

    } else {
      dPhiClassEle          [jecal] -> Fill ( dphi, wgt );
      dEtaClassEle          [jecal] -> Fill ( deta, wgt );
      EoPClassEle           [jecal] -> Fill ( eop, wgt );
      HoEClassEle           [jecal] -> Fill ( hoe, wgt );
      sigmaIEtaIEtaClassEle [jecal] -> Fill ( see, wgt );
      etaClassEle                   -> Fill ( eta, wgt );

      norm += wgt;
      if(see >= WP_inf[0] && see <= WP_sup[0])
        //         dphi >= WP_inf[1] && dphi <= WP_sup[1] &&
        //         deta >= WP_inf[2] && deta <= WP_sup[2] &&
        //         hoe >= WP_inf[3] && hoe <= WP_sup[3]) 
        sel += wgt;         

    }
  } // loop over events

  float eff = sel/norm;
  float erreff = sqrt(eff*(1-eff)/fabs(norm));

  std::cout << "EFFICIENCY = " << eff << " +/- " << erreff << std::endl; 
  std::cout << "(selected events = " << sel << " over total of " << norm << " events " << std::endl;

  char buffer[200];
  if( m_isMC ) sprintf(buffer,"pdfs_histograms_MC_x%d.root",m_factor);
  else sprintf(buffer,"pdfs_histograms_data_x%d.root",m_factor);
  TFile *fileOut = TFile::Open(buffer,"recreate");

  etaClassEle -> Write();  
  for (int jecal=0; jecal<2; jecal++) {
    dPhiClassEle[jecal]->Write();
    dEtaClassEle[jecal]->Write();
    EoPClassEle[jecal]->Write();           
    HoEClassEle[jecal]->Write();
    sigmaIEtaIEtaClassEle[jecal]->Write();
  }
  
  fileOut->Close();
}

void sPlotsPdfsComparison::bookHistosVariableBinning() {
  
  float dPhiBins[12] = {-0.1, -0.0348 , -0.0272 , -0.0196 , -0.012 , -0.0044 , 0.0032 , 0.0108 , 0.0184 , 
                        0.026 , 0.0336 , 0.1};

  float dEtaBins[12] = {-0.01 , -0.00424  , -0.00272 , -0.00196 , -0.0012 , -0.00044 , 0.00032 , 0.00108 , 
                        0.00184 , 0.0026 , 0.00412 , 0.01};

  float EoPBinsEB[12] = {0, 0.5, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 2, 4.0, 7., 10.};
  float EoPBinsEE[12] = {0, 0.5, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 2, 4.0, 7., 10.};

  float seeBinsEB[14] = {0., 0.005 , 0.0055 , 0.006 , 0.0065 , 0.007 , 0.0075 , 0.008 , 0.0085 , 0.009 , 0.0095, 0.01,
                         0.015, 0.020 };

  float seeBinsEE[14] = {0., 0.015 , 0.02, 0.021 , 0.022 , 0.023 , 0.024 , 0.025 , 0.026 , 0.027 , 0.028 , 0.029, 
                         0.03, 0.04};

  float HoEBins[6] = {0., 0.005, 0.01, 0.03, 0.06, 0.1 };

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  char histo[200];
  for (int iecal=0; iecal<2; iecal++) {
      
    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiClassEle[iecal] = new TH1F(histo, histo, 15, dPhiBins);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle[iecal] = new TH1F(histo, histo, 15, dEtaBins);
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEClassEle[iecal] = new TH1F(histo, histo, 5, HoEBins);

  }

  
  int iecal = 0;
  sprintf(histo,"EoPClass_electrons_%d",iecal);
  EoPClassEle[iecal] = new TH1F(histo, histo, 11, EoPBinsEB);
  sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
  sigmaIEtaIEtaClassEle[iecal] = new TH1F(histo, histo, 13, seeBinsEB);

  iecal = 1;
  sprintf(histo,"EoPClass_electrons_%d",iecal);
  EoPClassEle[iecal] = new TH1F(histo, histo, 11, EoPBinsEE);
  sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
  sigmaIEtaIEtaClassEle[iecal] = new TH1F(histo, histo, 13, seeBinsEE);

}



void sPlotsPdfsComparison::bookHistosFixedBinning() {
  
  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  char histo[200];

  sprintf(histo,"etaClass_electrons");
  etaClassEle = new TH1F(histo, histo, 10, -2.5, 2.5);
  etaClassEle ->Sumw2();

  for (int iecal=0; iecal<2; iecal++) {
      
    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiClassEle[iecal] = new TH1F(histo, histo, 15, -0.15, 0.15);
    dPhiClassEle[iecal]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle[iecal] = new TH1F(histo, histo, 15, -0.02, 0.02);
    dEtaClassEle[iecal]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEClassEle[iecal] = new TH1F(histo, histo, 15, 0., 0.15);
    HoEClassEle[iecal]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",iecal);
    EoPClassEle[iecal] = new TH1F(histo, histo, 15, 0., 10.);
    EoPClassEle[iecal]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
    sigmaIEtaIEtaClassEle[iecal] = new TH1F(histo, histo, 15, 0., 0.05);
    sigmaIEtaIEtaClassEle[iecal]->Sumw2();
  }
  
}

void sPlotsPdfsComparison::InitCuts() {

  // see, dphi, deta, H/E
  WP70_EB_inf.push_back(0.0);  
  WP70_EB_inf.push_back(-0.02);
  WP70_EB_inf.push_back(-0.006);
  WP70_EB_inf.push_back(0.00);

  WP70_EB_sup.push_back(0.01);
  WP70_EB_sup.push_back(0.02);
  WP70_EB_sup.push_back(0.006);
  WP70_EB_sup.push_back(0.02);

  WP70_EE_inf.push_back(0.0);
  WP70_EE_inf.push_back(-0.02);
  WP70_EE_inf.push_back(-0.003);
  WP70_EE_inf.push_back(0.00);

  WP70_EE_sup.push_back(0.03);
  WP70_EE_sup.push_back(0.02);
  WP70_EE_sup.push_back(0.003);
  WP70_EE_sup.push_back(0.0025);

}
