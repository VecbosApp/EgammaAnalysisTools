#include <iostream>
#include <string>

using namespace std;

double weight(string sample, double ngen, double xsec, double filtereff, double lumi = 10, double prescale = 1);

int main(int argc, char* argv[]) {
  
  double wantLumi = 1.0; //pb-1  
  double w;

  std::cout << "sample" << "\t" << "xsec" << "\t" << "weight" << std:: endl;
  
  w = weight("WjetsMADGRAPH", 9.7689e+06,  31314., 1.0, wantLumi);    // https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  w = weight("Zjets",         1.03492e+06,  3048., 1.0, wantLumi);    //       " "
  w = weight("TTbar",         1.2834e+06,  157.5,  1.0, wantLumi);    //       " " 

  w = weight("QCD_EMenriched_Pt20to30",  3.03098e+07, 235.5E+06, 0.0073, wantLumi);
  w = weight("QCD_EMenriched_Pt30to80",  3.36773e+07, 59.3E+06,  0.059,  wantLumi);
  w = weight("QCD_EMenriched_Pt80to170", 5.02093e+06, 0.906E+06, 0.148,  wantLumi);

  w = weight("QCD_BCtoE_Pt20to30", 2.26102e+06, 235.5E+06, 0.00046, wantLumi);
  w = weight("QCD_BCtoE_Pt30to80",      875597,  59.3E+06, 0.00234, wantLumi);
  w = weight("QCD_BCtoE_Pt80to170",       8674, 0.906E+06, 0.0104,  wantLumi);

  w = weight("QCD_Pt-20_TuneD6T", 2.35361e+07,  2.988E+08 , 1.0, wantLumi);   // lumi da cfg in DBS

  w = weight("PhotonJet_Pt0to15",    115265, 84.46E+06, 1.0, wantLumi);
  w = weight("PhotonJet_Pt15to20",   108560, 114700,    1.0, wantLumi);
  w = weight("PhotonJet_Pt20to30",    60000, 57180,     1.0, wantLumi);
  w = weight("PhotonJet_Pt30to50",   110000, 16520,     1.0, wantLumi);
  w = weight("PhotonJet_Pt50to80",   109730,  2723,     1.0, wantLumi);
  w = weight("PhotonJet_Pt80to120",  110827, 446.2,     1.0, wantLumi);
  w = weight("PhotonJet_Pt120to170", 122281, 84.43,     1.0, wantLumi);
  w = weight("PhotonJet_Pt170to300", 125128, 22.55,     1.0, wantLumi);
  w = weight("PhotonJet_Pt300to500", 107606, 1.545,     1.0, wantLumi);
  w = weight("PhotonJet_Pt500toInf",  56895, 0.0923,    1.0, wantLumi);

  w = weight("SingleTop_sChannel",  312055, 1.49,  1.0, wantLumi);
  w = weight("SingleTop_tChannel",  478593, 22.00, 1.0, wantLumi);   // https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf 
  w = weight("SingleTop_tWChannel", 16437,  10.6,  1.0, wantLumi);

  return 0;

}

double weight(string sample, double ngen, double xsec, double filtereff, 
	   double lumi, double prescale) {

  double W = xsec * filtereff * lumi * prescale / ngen;

  std::cout << sample << "\t" << xsec << "\t" << W << std:: endl;

  return W;

}
