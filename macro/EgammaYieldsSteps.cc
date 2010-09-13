#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

int UseWjetsCuts[12];
string WjetsCuts[12];
float Wj_Wsel[12];
float Zj_Wsel[12];
float ttj_Wsel[12];
float QCD_em_Wsel[12];
float QCD_bctoe_Wsel[12];
float SingleTop_Wsel[12];
float PhotonJet_Wsel[12];

float Wj_eff_Wsel[12];
float Zj_eff_Wsel[12];
float ttj_eff_Wsel[12];
float QCD_em_eff_Wsel[12];
float QCD_bctoe_eff_Wsel[12];
float SingleTop_eff_Wsel[12];
float PhotonJet_eff_Wsel[12];

float Wj_finaleff_Wsel;
float Zj_finaleff_Wsel;
float ttj_finaleff_Wsel;
float QCD_em_finaleff_Wsel;
float QCD_bctoe_finaleff_Wsel;
float SingleTop_finaleff_Wsel;
float PhotonJet_finaleff_Wsel;

// xsec
float Wj_MADGRAPH_xsec  = 31314.;  
float Zj_MADGRAPH_xsec  = 3048.;
float ttj_MADGRAPH_xsec = 157.5;

// here xsec = x-sec * filter_eff (pb)
float QCD_em20to30_xsec = 235.5E+06 * 0.0073;
float QCD_em30to80_xsec = 59.3E+06 * 0.059;
float QCD_em80to170_xsec = 0.906E+06 * 0.148;
float QCD_bctoe20to30_xsec = 235.5E+06 * 0.00046;
float QCD_bctoe30to80_xsec = 59.3E+06 * 0.00234;
float QCD_bctoe80to170_xsec = 0.906E+06 * 0.0104;

// xsec (pb)
float SingleTopS_xsec = 1.49;
float SingleTopT_xsec = 22.0;
float SingleTopTW_xsec = 10.6;

// xsec (pb)
float PhotonJet_Pt0to15_xsec = 84.46E+06;
float PhotonJet_Pt15to20_xsec = 114700;
float PhotonJet_Pt20to30_xsec = 57180;
float PhotonJet_Pt30to50_xsec = 16520;
float PhotonJet_Pt50to80_xsec = 2723;
float PhotonJet_Pt80to120_xsec = 446.2;
float PhotonJet_Pt120to170_xsec = 84.43;
float PhotonJet_Pt170to300_xsec = 22.55;
float PhotonJet_Pt300to500_xsec = 1.545;
float PhotonJet_Pt500toInf_xsec = 0.0923;

void computeYields(float lumi=100.) {

  char nametree[200];
  sprintf(nametree,"EVENT_COUNTER");
  cout << "nametree = " << nametree << endl;

  TChain *chains[22];
  for(int isample=0; isample<22; isample++) {
    chains[isample] = new TChain(nametree);
  }

  // W/Z+jets with MC truth request
  chains[0]->Add("results/WJetsMADGRAPH/WJets-madgraph/4/*counters.root");
  chains[1]->Add("results/ZJetsMADGRAPH/ZJets-madgraph/4/*counters.root");
  chains[2]->Add("results/TTbar/TTbarJets-madgraph/4/*counters.root");
  chains[3]->Add("results/QCD/QCD_EMEnriched_Pt20to30/4/*counters.root");
  chains[4]->Add("results/QCD/QCD_EMEnriched_Pt30to80/4/*counters.root");
  chains[5]->Add("results/QCD/QCD_EMEnriched_Pt80to170/4/*counters.root");
  chains[6]->Add("results/QCD/QCD_BCtoE_Pt20to30/4/*counters.root");
  chains[7]->Add("results/QCD/QCD_BCtoE_Pt30to80/4/*counters.root");
  chains[8]->Add("results/QCD/QCD_BCtoE_Pt80to170/4/*counters.root");
  chains[9]->Add("results/SingleTop/SingleTop_sChannel-madgraph/4/*counters.root");
  chains[10]->Add("results/SingleTop/SingleTop_tChannel-madgraph/4/*counters.root");
  chains[11]->Add("results/SingleTop/SingleTop_tWChannel-madgraph/4/*counters.root");
  chains[12]->Add("results/PhotonJet/PhotonJet_Pt0to15/4/*counters.root");
  chains[13]->Add("results/PhotonJet/PhotonJet_Pt15to20/4/*counters.root");
  chains[14]->Add("results/PhotonJet/PhotonJet_Pt20to30/4/*counters.root");
  chains[15]->Add("results/PhotonJet/PhotonJet_Pt30to50/4/*counters.root");
  chains[16]->Add("results/PhotonJet/PhotonJet_Pt50to80/4/*counters.root");
  chains[17]->Add("results/PhotonJet/PhotonJet_Pt80to120/4/*counters.root");
  chains[18]->Add("results/PhotonJet/PhotonJet_Pt120to170/4/*counters.root");
  chains[19]->Add("results/PhotonJet/PhotonJet_Pt170to300/4/*counters.root");
  chains[20]->Add("results/PhotonJet/PhotonJet_Pt300to500/4/*counters.root");
  chains[21]->Add("results/PhotonJet/PhotonJet_Pt500toInf/4/*counters.root");

  float nSelTot[12][22];
  for(int icut=0; icut<12; icut++) {
    for(int isample=0; isample<22; isample++) {
      nSelTot[icut][isample] = 0.0;
    }
  }

  int nCutsAna = 12;

  for(int isample=0; isample<22; isample++) {

    cout << "\tProcessing sample # " << isample << "..." << endl;

    Int_t           nCuts;
    Float_t         nSel[12];   //[nCuts]
    
    // List of branches
    TBranch        *b_nCuts;   //!
    TBranch        *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;

      // loop over cuts
      for(int icut=0; icut<nCuts; icut++) nSelTot[icut][isample] += nSel[icut];
    }
  }

  for(int icut=0; icut<nCutsAna; icut++) {

    // Nj = L * x-sec * eff (eff = Nj_fin / Nincl_ini)
    if(nSelTot[0][0]>0) Wj_Wsel[icut] = lumi * Wj_MADGRAPH_xsec * nSelTot[icut][0]/nSelTot[0][0];
    if(nSelTot[0][1]>0) Zj_Wsel[icut] = lumi * Zj_MADGRAPH_xsec * nSelTot[icut][1]/nSelTot[0][1];
    if(nSelTot[0][2]>0) ttj_Wsel[icut] = lumi * ttj_MADGRAPH_xsec * nSelTot[icut][2]/nSelTot[0][2];

    float qcdem_tmp=0.;
    if(nSelTot[0][3]>0) qcdem_tmp += lumi * QCD_em20to30_xsec * nSelTot[icut][3]/nSelTot[0][3];
    if(nSelTot[0][4]>0) qcdem_tmp += lumi * QCD_em30to80_xsec * nSelTot[icut][4]/nSelTot[0][4];
    if(nSelTot[0][5]>0) qcdem_tmp += lumi * QCD_em80to170_xsec * nSelTot[icut][5]/nSelTot[0][5];
    QCD_em_Wsel[icut] = qcdem_tmp;
    
    float qcdbctoe_tmp=0.;
    if(nSelTot[0][6]>0) qcdbctoe_tmp += lumi * QCD_bctoe20to30_xsec * nSelTot[icut][6]/nSelTot[0][6];
    if(nSelTot[0][7]>0) qcdbctoe_tmp += lumi * QCD_bctoe30to80_xsec * nSelTot[icut][7]/nSelTot[0][7];
    if(nSelTot[0][8]>0) qcdbctoe_tmp += lumi * QCD_bctoe80to170_xsec * nSelTot[icut][8]/nSelTot[0][8];
    QCD_bctoe_Wsel[icut] = qcdbctoe_tmp;
    
    float singletop_tmp=0.;
    if(nSelTot[0][9]>0) singletop_tmp += lumi * SingleTopS_xsec * nSelTot[icut][9]/nSelTot[0][9];
    if(nSelTot[0][10]>0) singletop_tmp += lumi * SingleTopT_xsec * nSelTot[icut][10]/nSelTot[0][10];
    if(nSelTot[0][11]>0) singletop_tmp += lumi * SingleTopTW_xsec * nSelTot[icut][11]/nSelTot[0][11];
    SingleTop_Wsel[icut] = singletop_tmp;
    
    float photonjet_tmp=0.;
    if(nSelTot[0][12]>0) photonjet_tmp += lumi * PhotonJet_Pt0to15_xsec * nSelTot[icut][12]/nSelTot[0][12];
    if(nSelTot[0][13]>0) photonjet_tmp += lumi * PhotonJet_Pt15to20_xsec * nSelTot[icut][13]/nSelTot[0][13];
    if(nSelTot[0][14]>0) photonjet_tmp += lumi * PhotonJet_Pt20to30_xsec * nSelTot[icut][14]/nSelTot[0][14];
    if(nSelTot[0][15]>0) photonjet_tmp += lumi * PhotonJet_Pt30to50_xsec * nSelTot[icut][15]/nSelTot[0][15];
    if(nSelTot[0][16]>0) photonjet_tmp += lumi * PhotonJet_Pt50to80_xsec * nSelTot[icut][16]/nSelTot[0][16];
    if(nSelTot[0][17]>0) photonjet_tmp += lumi * PhotonJet_Pt80to120_xsec * nSelTot[icut][17]/nSelTot[0][17];
    if(nSelTot[0][18]>0) photonjet_tmp += lumi * PhotonJet_Pt120to170_xsec * nSelTot[icut][18]/nSelTot[0][18];
    if(nSelTot[0][19]>0) photonjet_tmp += lumi * PhotonJet_Pt170to300_xsec * nSelTot[icut][19]/nSelTot[0][19];
    if(nSelTot[0][20]>0) photonjet_tmp += lumi * PhotonJet_Pt300to500_xsec * nSelTot[icut][20]/nSelTot[0][20];
    if(nSelTot[0][21]>0) photonjet_tmp += lumi * PhotonJet_Pt500toInf_xsec * nSelTot[icut][21]/nSelTot[0][21];
    PhotonJet_Wsel[icut] = photonjet_tmp;
    
    // EFFICIENCY
    int previousCut = icut-1;
    
    if (icut!=7) {
      if(icut>0 && nSelTot[previousCut][0]>0) Wj_eff_Wsel[icut] = nSelTot[icut][0] / nSelTot[previousCut][0];
      else Wj_eff_Wsel[icut] = 0.0;
      if(icut>0 && nSelTot[previousCut][1]>0) Zj_eff_Wsel[icut] = nSelTot[icut][1] / nSelTot[previousCut][1];
      else Zj_eff_Wsel[icut] = 0.0;
      if(icut>0 && nSelTot[previousCut][2]>0) ttj_eff_Wsel[icut] = nSelTot[icut][2] / nSelTot[previousCut][2];
      else ttj_eff_Wsel[icut] = 0.0;
      if(icut>0 && QCD_em_Wsel[previousCut]>0) QCD_em_eff_Wsel[icut] = QCD_em_Wsel[icut] / QCD_em_Wsel[previousCut];
      else QCD_em_eff_Wsel[icut] = 0.0;
      if(icut>0 && QCD_bctoe_Wsel[previousCut]>0) QCD_bctoe_eff_Wsel[icut] = QCD_bctoe_Wsel[icut] / QCD_bctoe_Wsel[previousCut];
      else QCD_bctoe_eff_Wsel[icut] = 0.0;
      if(icut>0 && SingleTop_Wsel[previousCut]>0) SingleTop_eff_Wsel[icut] = SingleTop_Wsel[icut] / SingleTop_Wsel[previousCut];
      else SingleTop_eff_Wsel[icut] = 0.0;
      if(icut>0 && PhotonJet_Wsel[previousCut]>0) PhotonJet_eff_Wsel[icut] = PhotonJet_Wsel[icut] / PhotonJet_Wsel[previousCut];
      else PhotonJet_eff_Wsel[icut] = 0.0;
    } else {
      if(nSelTot[2][0]>0) Wj_eff_Wsel[icut] = nSelTot[icut][0] / nSelTot[2][0];
      else Wj_eff_Wsel[icut] = 0.0;
      if(nSelTot[2][1]>0) Zj_eff_Wsel[icut] = nSelTot[icut][1] / nSelTot[2][1];
      else Zj_eff_Wsel[icut] = 0.0;
      if(nSelTot[2][2]>0) ttj_eff_Wsel[icut] = nSelTot[icut][2] / nSelTot[2][2];
      else ttj_eff_Wsel[icut] = 0.0;
      if(QCD_em_Wsel[2]>0) QCD_em_eff_Wsel[icut] = QCD_em_Wsel[icut] / QCD_em_Wsel[2];
      else QCD_em_eff_Wsel[icut] = 0.0;
      if(QCD_bctoe_Wsel[2]>0) QCD_bctoe_eff_Wsel[icut] = QCD_bctoe_Wsel[icut] / QCD_bctoe_Wsel[2];
      else QCD_bctoe_eff_Wsel[icut] = 0.0;
      if(SingleTop_Wsel[2]>0) SingleTop_eff_Wsel[icut] = SingleTop_Wsel[icut] / SingleTop_Wsel[2];
      else SingleTop_eff_Wsel[icut] = 0.0;
      if(PhotonJet_Wsel[2]>0) PhotonJet_eff_Wsel[icut] = PhotonJet_Wsel[icut] / PhotonJet_Wsel[2];
      else PhotonJet_eff_Wsel[icut] = 0.0;
    }
  } 
  if(nSelTot[1][0]>0) Wj_finaleff_Wsel = nSelTot[nCutsAna-1][0] / nSelTot[1][0]; // normalize to MC truth events
  else Wj_finaleff_Wsel = 0.0;
  if(nSelTot[0][1]>0) Zj_finaleff_Wsel = nSelTot[nCutsAna-1][1] / nSelTot[0][1]; 
  else Zj_finaleff_Wsel = 0.0;
  if(nSelTot[0][2]>0) ttj_finaleff_Wsel = nSelTot[nCutsAna-1][2] / nSelTot[0][2]; 
  else ttj_finaleff_Wsel = 0.0;
  float qcdem_ini = nSelTot[0][3]+nSelTot[0][4]+nSelTot[0][5];
  float qcdem_fin = nSelTot[nCutsAna-1][3]+nSelTot[nCutsAna-1][4]+nSelTot[nCutsAna-1][5];
  if(qcdem_ini>0) QCD_em_finaleff_Wsel = qcdem_fin / qcdem_ini;
  else QCD_em_finaleff_Wsel = 0.0;
  float qcdbctoe_ini = nSelTot[0][6]+nSelTot[0][7]+nSelTot[0][8];
  float qcdbctoe_fin = nSelTot[nCutsAna-1][6]+nSelTot[nCutsAna-1][7]+nSelTot[nCutsAna-1][8];
  if(qcdbctoe_ini>0) QCD_bctoe_finaleff_Wsel = qcdbctoe_fin / qcdbctoe_ini;
  else QCD_bctoe_finaleff_Wsel = 0.0;
  float singletop_ini = nSelTot[0][9]+nSelTot[0][10]+nSelTot[0][11];
  float singletop_fin = nSelTot[nCutsAna-1][9]+nSelTot[nCutsAna-1][10]+nSelTot[nCutsAna-1][11];
  if(singletop_ini>0) SingleTop_finaleff_Wsel = singletop_fin / singletop_ini;
  else SingleTop_finaleff_Wsel = 0.0;
  float photonjet_ini = nSelTot[0][12]+nSelTot[0][13]+nSelTot[0][14]+nSelTot[0][15]+nSelTot[0][16]+nSelTot[0][17]+nSelTot[0][18]+nSelTot[0][19]+nSelTot[0][20]+nSelTot[0][21];
  float photonjet_fin = nSelTot[nCutsAna-1][12]+nSelTot[nCutsAna-1][13]+nSelTot[nCutsAna-1][14]+nSelTot[nCutsAna-1][15]+
    nSelTot[nCutsAna-1][16]+nSelTot[nCutsAna-1][17]+nSelTot[nCutsAna-1][18]+nSelTot[nCutsAna-1][19]+nSelTot[nCutsAna-1][20]+nSelTot[nCutsAna-1][21];
  if(photonjet_ini>0) PhotonJet_finaleff_Wsel = photonjet_fin / photonjet_ini;
  else PhotonJet_finaleff_Wsel = 0.0;
}

void setupCuts() {
  
  for(int i=0; i<12; i++) {
    UseWjetsCuts[i] = 1;
  }
  
  WjetsCuts[0]="event";
  WjetsCuts[1]="pthat";
  WjetsCuts[2]="HLT";
  WjetsCuts[3]="probe found (nio)";
  WjetsCuts[4]="tag found (nio)";
  WjetsCuts[5]="eletot (nio)";
  WjetsCuts[6]="deltaPhi (nio)";
  WjetsCuts[7]="tag and probe found";
  WjetsCuts[8]="Z mass veto";
  WjetsCuts[9]="tracker not isol";
  WjetsCuts[10]="ecal not isol";
  WjetsCuts[11]="full selection";
}

void printLatex(float lumi) {

  setupCuts();

  computeYields(lumi);
  
  char namefile[200];
  sprintf(namefile,"yieldsByCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::trunc);
  textfile.precision(0);
  
  textfile << "\\documentclass{article}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl;
  
  // W+jets events
  textfile << "\\begin{sidewaystable}[p]" << endl
	   << "\\begin{small}" << endl
	   << "\\begin{center}" << endl
	   << "\\begin{tabular}{|c|c|c|c|c|c|c|c}" << endl
	   << "\\hline" << endl
	   << "selection & W$(e \\nu)$+jets & Z+jets & $t\\bar{t}$ & single $t$ & QCD (uds) & QCD (bc) & $\\gamma$+jets \t\\\\" << endl
	   << "\\hline" << endl; 
  
  for(int icut=0; icut<12; icut++) {
    
    if(!UseWjetsCuts[icut]) continue;
    
    textfile << WjetsCuts[icut] << "\t&\t";

    if (icut!=3 && icut!=4 && icut!=5 && icut!=6) {
      textfile << fixed
	       << Wj_Wsel[icut]        << " (" << 100. * Wj_eff_Wsel[icut] << "\\%)" << "\t&\t"
	       << Zj_Wsel[icut]        << " (" << 100. * Zj_eff_Wsel[icut] << "\\%)" << "\t&\t"
	       << ttj_Wsel[icut]       << " (" << 100. * ttj_eff_Wsel[icut] << "\\%)" << "\t&\t"
	       << SingleTop_Wsel[icut] << " (" << 100. * SingleTop_eff_Wsel[icut] << "\\%)" << "\t&\t"
	       << QCD_em_Wsel[icut]    << " (" << 100. * QCD_em_eff_Wsel[icut] << "\\%)" << "\t&\t"
	       << QCD_bctoe_Wsel[icut] << " (" << 100. * QCD_bctoe_eff_Wsel[icut] << "\\%)" << "\t&\t"
	       << PhotonJet_Wsel[icut] << " (" << 100. * PhotonJet_eff_Wsel[icut] << "\\%)" << "\t\\\\" << endl;
    } else {
      textfile << fixed
	       << Wj_Wsel[icut]        << "\t&\t"
	       << Zj_Wsel[icut]        << "\t&\t"
	       << ttj_Wsel[icut]       << "\t&\t"
	       << SingleTop_Wsel[icut] << "\t&\t"
	       << QCD_em_Wsel[icut]    << "\t&\t"
	       << QCD_bctoe_Wsel[icut] << "\t&\t"
	       << PhotonJet_Wsel[icut] << "\t\\\\" << endl;
    }
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total " << "\t&\t"
	   << Wj_Wsel[11] << " (" << 100. * Wj_finaleff_Wsel << "\\%)" << "\t&\t"
	   << Zj_Wsel[11] << " (" << 100. * Zj_finaleff_Wsel << "\\%)" << "\t&\t"
	   << ttj_Wsel[11] << " (" << 100. * ttj_finaleff_Wsel << "\\%)" << "\t&\t"
	   << SingleTop_Wsel[11] << " (" << 100. * SingleTop_finaleff_Wsel << "\\%)" << "\t&\t"
	   << QCD_em_Wsel[11] << " (" << 100. * QCD_em_finaleff_Wsel << "\\%)" << "\t&\t"
	   << QCD_bctoe_Wsel[11] << " (" << 100. * QCD_bctoe_finaleff_Wsel << "\\%)" << "\t&\t"
	   << PhotonJet_Wsel[11] << " (" << 100. * PhotonJet_finaleff_Wsel << "\\%)"
	   << "\t\\\\" << endl;
  
  textfile << "\\hline" << endl
	   << "\\end{tabular}" << endl
	   << "\\end{center}" << endl
	   << "\\end{small}" << endl
	   << "\\end{sidewaystable}" << endl;
    
  textfile << "\\end{document}" << endl;
}
