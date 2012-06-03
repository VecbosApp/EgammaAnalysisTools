#include "TH2F.h"
#include "TFile.h"
void mod() {
  
  Float_t xbins[8] = {7,10,15,20,30,40,50,1000};
  Float_t ybins[6] = {0.,0.8,1.4442,1.566,2.0,2.5};

  TH2F *hnew = new TH2F("heff","heff",7,xbins,5,ybins);

  Float_t eff[7][5];

  eff[0][0] = 0.980;
  eff[0][1] = 0.739;
  eff[0][2] = 1.000; // too few stat
  eff[0][3] = 0.526;
  eff[0][4] = 0.584;

  eff[1][0] = 1.008;
  eff[1][1] = 0.900;
  eff[1][2] = 1.000; // too few stat
  eff[1][3] = 1.050;
  eff[1][4] = 1.032;

  eff[2][0] = 1.045;
  eff[2][1] = 1.024;
  eff[2][2] = 1.000; // too few stat
  eff[2][3] = 1.070;
  eff[2][4] = 1.090;
  
  eff[3][0] = 1.016;
  eff[3][1] = 1.000;
  eff[3][2] = 1.040;
  eff[3][3] = 1.004;
  eff[3][4] = 1.009;

  eff[4][0] = 1.004;
  eff[4][1] = 0.998;
  eff[4][2] = 0.987;
  eff[4][3] = 0.982;
  eff[4][4] = 0.979;

  eff[5][0] = 1.000;
  eff[5][1] = 0.994;
  eff[5][2] = 0.992;
  eff[5][3] = 0.980;
  eff[5][4] = 0.974;

  eff[6][0] = 1.000;
  eff[6][1] = 0.991;
  eff[6][2] = 0.992;
  eff[6][3] = 0.984;
  eff[6][4] = 0.977;


  Float_t efferr[7][5];

  efferr[0][0] = 0.097;
  efferr[0][1] = 0.082;
  efferr[0][2] = 0.073; // too few stat
  efferr[0][3] = 0.080;
  efferr[0][4] = 0.106;

  efferr[1][0] = 0.077;
  efferr[1][1] = 0.078;
  efferr[1][2] = 0.095; // too few stat
  efferr[1][3] = 0.072;
  efferr[1][4] = 0.078;

  efferr[2][0] = 0.066;
  efferr[2][1] = 0.064;
  efferr[2][2] = 0.071; // too few stat
  efferr[2][3] = 0.060;
  efferr[2][4] = 0.067;
  
  efferr[3][0] = 0.060;
  efferr[3][1] = 0.060;
  efferr[3][2] = 0.067;
  efferr[3][3] = 0.061;
  efferr[3][4] = 0.063;

  efferr[4][0] = 0.057;
  efferr[4][1] = 0.058;
  efferr[4][2] = 0.062;
  efferr[4][3] = 0.057;
  efferr[4][4] = 0.061;

  efferr[5][0] = 0.056;
  efferr[5][1] = 0.057;
  efferr[5][2] = 0.059;
  efferr[5][3] = 0.058;
  efferr[5][4] = 0.061;

  efferr[6][0] = 0.055;
  efferr[6][1] = 0.056;
  efferr[6][2] = 0.058;
  efferr[6][3] = 0.057;
  efferr[6][4] = 0.058;

  for(int ix=0;ix<7;ix++) {
    for(int iy=0;iy<5;iy++) {
      hnew->SetBinContent(ix+1,iy+1,eff[ix][iy]);
      hnew->SetBinError(ix+1,iy+1,efferr[ix][iy]);
    }
  }

  TFile *filenew = TFile::Open("hel_sf_2011.root","recreate");
  hnew->Write();
  filenew->Close();

}
