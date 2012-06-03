#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TEventList.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

using namespace std;
void makeList(const char* cut) {
  TFile *f1 = new TFile("../results_data/electrons_zeemc_hzzidbitsFriend.root");
  TTree *ntuple = (TTree*) f1->Get("eleIDdir/T1");
  ntuple->AddFriend("eleIDdir/T1 = eleIDdir/T1","../results_data/electrons_zeemc.root");
  ntuple->Draw(">>elist",cut);
  TEventList *elist = (TEventList*)gDirectory->Get("elist");
  TFile ef("elist.root","recreate");
  elist->Write();
}

void makeSmall() {
  TFile *f = new TFile("elist.root");
  TEventList *elist = (TEventList*)f->Get("elist");
  TFile *f1 = new TFile("../results_data/electrons_zeemc_hzzidbitsFriend.root");
  TTree *ntuple = (TTree*) f1->Get("eleIDdir/T1");
  ntuple->SetEventList(elist);
  TFile *f2 = new TFile("../results_data/electrons_zeemc_hzzidbitsFriend_promptrate.root","recreate");
  f2->mkdir("eleIDdir");
  TTree *small = ntuple->CopyTree("");
  f2->cd("eleIDdir");
  small->Write();
  small->Print();
}
