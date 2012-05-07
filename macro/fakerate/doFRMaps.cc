#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include <iostream>
#include <vector>

enum {
  kEleLoose = 0,
  kEleTight
};

enum {
  kMuPfIso = 0,
  kMuMvaPfIso
};
  
using namespace std;

void doEleFRMapsHzz4l(int wp) {

  if(wp>1) {
    cout << "wp must be kLoose / kTight" << endl;
    return;
  }
  TString wpstr = (wp==kEleLoose) ? TString("Loose") : TString("Tight");

  TFile *file = TFile::Open("fakerates_zee1fake.root");

  TH1F *barrel = (TH1F*)file->Get(TString("NoTrgElehzzMva")+wpstr+TString("PtBarrel_Eff"));
  TH1F *endcap = (TH1F*)file->Get(TString("NoTrgElehzzMva")+wpstr+TString("PtEndcap_Eff"));

  std::vector<TH1F*> fakerates1D;
  fakerates1D.push_back(barrel);
  fakerates1D.push_back(endcap);

  Float_t LowerPt[11] = {0.0,5.0,7.0,10.0,14.0,20.0,25.0,30.0,40.0,50.0,80.0};
  Float_t LowerEta[3] = {0.0,1.479,2.5};


  TString namefile = TString("hzz4lfr")+wpstr+TString(".root");
  TFile *targetf = TFile::Open(namefile,"recreate");

  TH2F *map = new TH2F("frmap","",10,LowerPt,2,LowerEta);
  for(int ieta=0; ieta<2; ieta++) {
    for(int ipt=0;ipt<11;ipt++) {
      float fr=fakerates1D[ieta]->GetBinContent(ipt);
      float fre=fakerates1D[ieta]->GetBinError(ipt);
      map->SetBinContent(ipt,ieta+1,fr);
      map->SetBinError(ipt,ieta+1,fre);
    }
  }

  map->Write();
  targetf->Close();

}

void doMuFRMapsHzz4l(int wp) {

  if(wp>1) {
    cout << "wp must be kPfIso / kMvaPfIso" << endl;
    return;
  }
  TString wpstr = (wp==kMuPfIso) ? TString("PfIso") : TString("MvaPfIso");

  TFile *file = TFile::Open("fakerates_zll1m.root");

  TH1F *barrel = (TH1F*)file->Get(TString("NoTrgMuonhzz")+wpstr+TString("PtBarrel_Eff"));
  TH1F *endcap = (TH1F*)file->Get(TString("NoTrgMuonhzz")+wpstr+TString("PtEndcap_Eff"));

  std::vector<TH1F*> fakerates1D;
  fakerates1D.push_back(barrel);
  fakerates1D.push_back(endcap);

  Float_t LowerPt[11] = {0.0,5.0,7.0,10.0,14.0,20.0,25.0,30.0,40.0,50.0,80.0};
  Float_t LowerEta[3] = {0.0,1.479,2.5};


  TString namefile = TString("hzz4lmufr")+wpstr+TString(".root");
  TFile *targetf = TFile::Open(namefile,"recreate");

  TH2F *map = new TH2F("frmap","",10,LowerPt,2,LowerEta);
  for(int ieta=0; ieta<2; ieta++) {
    for(int ipt=0;ipt<11;ipt++) {
      float fr=fakerates1D[ieta]->GetBinContent(ipt);
      float fre=fakerates1D[ieta]->GetBinError(ipt);
      map->SetBinContent(ipt,ieta+1,fr);
      map->SetBinError(ipt,ieta+1,fre);
    }
  }

  map->Write();
  targetf->Close();

}
