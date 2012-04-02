#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TEventList.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

using namespace std;
void makeList() {
  TFile *f1 = new TFile("../results_data/electrons_hzzidbitsFriend.root");
  TTree *ntuple = (TTree*) f1->Get("eleIDdir/T1");
  ntuple->Draw(">>elist","denom==1");
  TEventList *elist = (TEventList*)gDirectory->Get("elist");
  TFile ef("elist.root","recreate");
  elist->Write();
}

void makeSmall() {
  TFile *f = new TFile("elist.root");
  TEventList *elist = (TEventList*)f->Get("elist");
  TFile *f1 = new TFile("../results_data/electrons_hzzidbitsFriend.root");
  TTree *ntuple = (TTree*) f1->Get("eleIDdir/T1");
  ntuple->SetEventList(elist);
  TFile *f2 = new TFile("../results_data/electrons_hzzidbitsFriend_promptrate.root","recreate");
  f2->mkdir("eleIDdir");
  TTree *small = ntuple->CopyTree("");
  f2->cd("eleIDdir");
  small->Write();
  small->Print();
}
