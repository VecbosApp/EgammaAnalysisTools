#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"

void packTH2F(const char *filename) {

  TFile *file = TFile::Open(filename);
  file->cd("tpTreeElEl/PASSING_all/fit_eff_plots");

  float ptbins[8] = { 7, 10, 15, 20, 30, 40, 50, 100 };
  float etabins[6] = { 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5 };

  TH2F *map = new TH2F("heff","",7,ptbins,5,etabins);

  for(int ieta=0; ieta<5; ++ieta) {
    char plot[1000];
    sprintf(plot,"pt_PLOT_abseta_bin%d_&_tag_HLT_Ele20_Ele4_TnP_Ele20Leg_pass",ieta);
    TCanvas *c = (TCanvas*)gDirectory->Get(plot);
    for(int ipt=0; ipt<7; ++ipt) {
      TGraphAsymmErrors* fit = (TGraphAsymmErrors*)c->FindObject("hxy_fit_eff");
      double x = fit->GetX()[ipt], y = fit->GetY()[ipt], eyl = fit->GetErrorYlow(ipt), eyh = fit->GetErrorYhigh(ipt);
      map->SetBinContent(ipt+1,ieta+1,y);
      map->SetBinError(ipt+1,ieta+1,(eyl+eyh)/2.0);
    }
  }

  TFile *filemap = TFile::Open("mapeff.root","recreate");
  map->Write();
  filemap->Close();

}
