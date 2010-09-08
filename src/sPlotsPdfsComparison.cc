#define sPlotsPdfsComparison_cxx
#include "include/sPlotsPdfsComparison.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

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
  float selEB, selEE, normEB, normEE;
  selEB = selEE = normEB = normEE = 0;
  
  Long64_t nentries = fChain->GetEntries();
  
  int firstEvent = 0;
  int lastEvent  = nentries;
  cout << firstEvent << " " << lastEvent << endl;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=firstEvent; jentry<lastEvent;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // da cambiare (solo x dati, selezionare sgn o fondo)
    float weight = (m_doSignal) ? N_sig_sw : N_qcd_sw;
    float wgt = (m_isMC) ? f_weight : weight;
    // float wgt = (m_isMC) ? 1.0 : N_qcd_sw;
    
    int jecal;
    if(m_isMC) {
      if(wgt>500) continue;
      if (fabs(f_eta)<1.479)  jecal = 0;
      else if (fabs(f_eta)>=1.479) jecal = 1;
      else continue;
      if(m_doSignal) {
        if(f_nJets>0 || f_pfmet/f_pt<0.3) continue;
      }
    } else {
      if(see<1e-4) continue; // remove residual spikes (if the cleaning is not done)
      if (fabs(eta)<1.479)  jecal = 0;
      else if (fabs(eta)>=1.479) jecal = 1;
      else continue;
      if(m_doSignal) {
        if(nJets>0 || pfmet/pt<0.3) continue;
      }
    }

    vector<float> WP_inf, WP_sup;
    if(jecal == 0) {
      WP_inf = WP70_EB_inf;
      WP_sup = WP70_EB_sup;
      normEB += wgt;
    } else {
      WP_inf = WP70_EE_inf;
      WP_sup = WP70_EE_sup;
      normEE += wgt;
    }

    if(m_isMC) {
      dPhiClassEle          [jecal] -> Fill ( f_dphi, wgt );
      dEtaClassEle          [jecal] -> Fill ( f_deta, wgt );
      EoPClassEle           [jecal] -> Fill ( f_eop, wgt );
      HoEClassEle           [jecal] -> Fill ( f_hoe, wgt );
      sigmaIEtaIEtaClassEle [jecal] -> Fill ( f_see, wgt );
      fbremClassEle         [jecal] -> Fill ( f_fbrem, wgt );
      etaClassEle                   -> Fill ( f_eta, wgt );
      phiClassEle           [jecal] -> Fill ( f_phi, wgt );
      chargeClassEle        [jecal] -> Fill ( f_charge, wgt );
      
      if(f_see >= WP_inf[0] && f_see <= WP_sup[0] &&
         f_dphi >= WP_inf[1] && f_dphi <= WP_sup[1] &&
         f_deta >= WP_inf[2] && f_deta <= WP_sup[2] &&
         f_hoe >= WP_inf[3] && f_hoe <= WP_sup[3])
        if(jecal == 0) selEB += wgt;         
        else selEE += wgt;

    } else {
      dPhiClassEle          [jecal] -> Fill ( dphi, wgt );
      dEtaClassEle          [jecal] -> Fill ( deta, wgt );
      EoPClassEle           [jecal] -> Fill ( eop, wgt );
      HoEClassEle           [jecal] -> Fill ( hoe, wgt );
      sigmaIEtaIEtaClassEle [jecal] -> Fill ( see, wgt );
      fbremClassEle         [jecal] -> Fill ( fbrem, wgt );
      etaClassEle                   -> Fill ( eta, wgt );
      phiClassEle           [jecal] -> Fill ( phi, wgt );
      chargeClassEle        [jecal] -> Fill ( charge, wgt );

      if(see >= WP_inf[0] && see <= WP_sup[0] &&
         dphi >= WP_inf[1] && dphi <= WP_sup[1] &&
         deta >= WP_inf[2] && deta <= WP_sup[2] &&
         hoe >= WP_inf[3] && hoe <= WP_sup[3])
        if(jecal == 0) selEB += wgt;         
        else selEE += wgt;
    }
  } // loop over events

  float effEB = selEB/normEB;
  float effEE = selEE/normEE;
  float erreffEB = sqrt(effEB*(1-effEB)/fabs(normEB));
  float erreffEE = sqrt(effEE*(1-effEE)/fabs(normEE));

  std::cout << "EB EFFICIENCY = " << effEB << " +/- " << erreffEB << std::endl; 
  std::cout << "(selected events = " << selEB << " over total of " << normEB << " events " << std::endl;

  std::cout << "EE EFFICIENCY = " << effEE << " +/- " << erreffEE << std::endl; 
  std::cout << "(selected events = " << selEE << " over total of " << normEE << " events " << std::endl;

  char buffer[200];
  if( m_isMC ) sprintf(buffer,"pdfs_histograms_MC.root");
  else sprintf(buffer,"pdfs_histograms_data.root");
  TFile *fileOut = TFile::Open(buffer,"recreate");

  etaClassEle->Write();
  
  for (int jecal=0; jecal<2; jecal++) {
    dPhiClassEle[jecal]->Write();
    dEtaClassEle[jecal]->Write();
    EoPClassEle[jecal]->Write();           
    HoEClassEle[jecal]->Write();
    sigmaIEtaIEtaClassEle[jecal]->Write();
    fbremClassEle[jecal]->Write();
    phiClassEle[jecal]->Write();    
    chargeClassEle[jecal]->Write();    
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
    dPhiClassEle[iecal] = new TH1F(histo, histo, 10, dPhiBins);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle[iecal] = new TH1F(histo, histo, 10, dEtaBins);
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
  //  for (int iecal=0; iecal<2; iecal++) {

  sprintf(histo,"etaClass_electrons");
  etaClassEle = new TH1F(histo, histo, 30, -3.0, 3.0);
  etaClassEle->Sumw2();

  if(m_doSignal) {
    // EB histograms
    sprintf(histo,"dPhiClass_electrons_%d",0);
    dPhiClassEle[0] = new TH1F(histo, histo, 50, -0.05, 0.05);
    dPhiClassEle[0]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",0);
    dEtaClassEle[0] = new TH1F(histo, histo, 50, -0.008, 0.008);
    dEtaClassEle[0]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",0);
    HoEClassEle[0] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEClassEle[0]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",0);
    EoPClassEle[0] = new TH1F(histo, histo, 50, 0., 3.);
    EoPClassEle[0]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",0);
    sigmaIEtaIEtaClassEle[0] = new TH1F(histo, histo, 50, 0., 0.03);
    sigmaIEtaIEtaClassEle[0]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",0);
    fbremClassEle[0] = new TH1F(histo, histo, 40, -1.0, 1.0);
    fbremClassEle[0]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",0);
    phiClassEle[0] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiClassEle[0]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",0);
    chargeClassEle[0] = new TH1F(histo, histo, 2, -2., 2.);
    chargeClassEle[0]->Sumw2();

    // EE histograms
    sprintf(histo,"dPhiClass_electrons_%d",1);
    dPhiClassEle[1] = new TH1F(histo, histo, 50, -0.1, 0.1);
    dPhiClassEle[1]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",1);
    dEtaClassEle[1] = new TH1F(histo, histo, 50, -0.03, 0.03);
    dEtaClassEle[1]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",1);
    HoEClassEle[1] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEClassEle[1]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",1);
    EoPClassEle[1] = new TH1F(histo, histo, 50, 0., 5.);
    EoPClassEle[1]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",1);
    sigmaIEtaIEtaClassEle[1] = new TH1F(histo, histo, 50, 0.01, 0.05);
    sigmaIEtaIEtaClassEle[1]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",1);
    fbremClassEle[1] = new TH1F(histo, histo, 30, -1.0, 1.0);
    fbremClassEle[1]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",1);
    phiClassEle[1] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiClassEle[1]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",1);
    chargeClassEle[1] = new TH1F(histo, histo, 2, -2., 2.);
    chargeClassEle[1]->Sumw2();
  } else {
    // EB histograms
    sprintf(histo,"dPhiClass_electrons_%d",0);
    dPhiClassEle[0] = new TH1F(histo, histo, 50, -0.15, 0.15);
    dPhiClassEle[0]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",0);
    dEtaClassEle[0] = new TH1F(histo, histo, 50, -0.02, 0.02);
    dEtaClassEle[0]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",0);
    HoEClassEle[0] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEClassEle[0]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",0);
    EoPClassEle[0] = new TH1F(histo, histo, 50, 0., 10.);
    EoPClassEle[0]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",0);
    sigmaIEtaIEtaClassEle[0] = new TH1F(histo, histo, 50, 0., 0.03);
    sigmaIEtaIEtaClassEle[0]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",0);
    fbremClassEle[0] = new TH1F(histo, histo, 50, -1.0, 1.0);
    fbremClassEle[0]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",0);
    phiClassEle[0] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiClassEle[0]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",0);
    chargeClassEle[0] = new TH1F(histo, histo, 2, -2., 2.);
    chargeClassEle[0]->Sumw2();

    // EE histograms
    sprintf(histo,"dPhiClass_electrons_%d",1);
    dPhiClassEle[1] = new TH1F(histo, histo, 50, -0.3, 0.3);
    dPhiClassEle[1]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",1);
    dEtaClassEle[1] = new TH1F(histo, histo, 50, -0.02, 0.02);
    dEtaClassEle[1]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",1);
    HoEClassEle[1] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEClassEle[1]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",1);
    EoPClassEle[1] = new TH1F(histo, histo, 50, 0., 10.);
    EoPClassEle[1]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",1);
    sigmaIEtaIEtaClassEle[1] = new TH1F(histo, histo, 50, 0.01, 0.08);
    sigmaIEtaIEtaClassEle[1]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",1);
    fbremClassEle[1] = new TH1F(histo, histo, 50, -1.0, 1.0);
    fbremClassEle[1]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",1);
    phiClassEle[1] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiClassEle[1]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",1);
    chargeClassEle[1] = new TH1F(histo, histo, 2, -2., 2.);
    chargeClassEle[1]->Sumw2();
  }
  
}

void sPlotsPdfsComparison::InitCuts() {

  // see, dphi, deta, H/E
  WP70_EB_inf.push_back(0.0);  
  WP70_EB_inf.push_back(-0.06);
  WP70_EB_inf.push_back(-0.004);
  WP70_EB_inf.push_back(0.00);

  WP70_EB_sup.push_back(0.01);
  WP70_EB_sup.push_back(0.06);
  WP70_EB_sup.push_back(0.004);
  WP70_EB_sup.push_back(0.04);

  WP70_EE_inf.push_back(0.0);
  WP70_EE_inf.push_back(-0.03);
  WP70_EE_inf.push_back(-0.007);
  WP70_EE_inf.push_back(0.00);

  WP70_EE_sup.push_back(0.03);
  WP70_EE_sup.push_back(0.03);
  WP70_EE_sup.push_back(0.007);
  WP70_EE_sup.push_back(0.025);

}
