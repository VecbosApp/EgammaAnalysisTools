#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

void countEvents() {

  char nametree[200];
  sprintf(nametree,"EVENT_COUNTER");

  cout << "nametree = " << nametree << endl;
  
  TChain *chains[22];
  for(int isample=0; isample<22; isample++) {
    chains[isample] = new TChain(nametree);
  }

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

  cout << "chains added. " << endl;

  std::vector<std::string> sampleName;

  sampleName.push_back("WJetsMADGRAPH");

  sampleName.push_back("ZJetsMADGRAPH");

  sampleName.push_back("TTbar_mcatnlo");

  sampleName.push_back("QCD_QCD_EMEnriched_Pt20to30");
  sampleName.push_back("QCD_QCD_EMEnriched_Pt30to80");
  sampleName.push_back("QCD_QCD_EMEnriched_Pt80to170");

  sampleName.push_back("QCD_QCD_BCtoE_Pt20to30");
  sampleName.push_back("QCD_QCD_BCtoE_Pt30to80");
  sampleName.push_back("QCD_QCD_BCtoE_Pt80to170");

  sampleName.push_back("SingleTop_SingleTop_sChannel_madgraph");
  sampleName.push_back("SingleTop_SingleTop_tChannel_madgraph");
  sampleName.push_back("SingleTop_SingleTop_tWChannel_madgraph");

  sampleName.push_back("PhotonJet_Pt0to15");
  sampleName.push_back("PhotonJet_Pt15to20");
  sampleName.push_back("PhotonJet_Pt20to30");
  sampleName.push_back("PhotonJet_Pt30to50");
  sampleName.push_back("PhotonJet_Pt50to80");
  sampleName.push_back("PhotonJet_Pt80to120");
  sampleName.push_back("PhotonJet_Pt120to170");
  sampleName.push_back("PhotonJet_Pt170to300");
  sampleName.push_back("PhotonJet_Pt300to500");
  sampleName.push_back("PhotonJet_Pt500toInf");



  float nEv[22];

  for(int isample=0; isample<22; isample++) {
    nEv[isample] = 0.0;
  }

  for(int isample=0; isample<22; isample++) {

    cout << "\tProcessing sample # " << isample << "..." << endl;

    Int_t           nCuts;
    Float_t         nSel[20];   //[nCuts]
    
    // List of branches
    TBranch        *b_nCuts;   //!
    TBranch        *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;

      nEv[isample] += nSel[0];
    }
  }

  for(int isample=0; isample<22; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
  }
  
}

